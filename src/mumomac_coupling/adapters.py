"""Score adapters: reuse the official scorer, never infer algorithms from scales."""
from dataclasses import dataclass,field,asdict
import json,hashlib
import numpy as np,pandas as pd
from .registry import ORDER,registry,digest
UNKNOWN={'','unknown','nan','na','none','null','missing','unverified'}
def known(v):return v is not None and str(v).strip().lower() not in UNKNOWN
def confirmed(v):return str(v).lower() in {'true','1','confirmed','verified'}
@dataclass
class ScoreSpec:
 namespace:str
 signature_version:str
 species:str
 scoring_method:str
 scoring_version:str
 program_hashes:dict
 reference_version:str='1.0.0'
 parameters:dict=field(default_factory=dict)
 mapping_version:str='mouse_native'
 mapping_hash:str=''
 def validate(self):
  if self.reference_version!='1.0.0' or self.signature_version not in ['original','expanded'] or self.species not in ['mouse','human']:raise ValueError('Unsupported frozen reference/version/species.')
  expected={f['feature_id']:f['signatures'][self.signature_version]['gene_set_hash'] for f in registry()['features']}
  if self.program_hashes!=expected:raise ValueError('Program gene-set hashes do not match frozen definitions.')
  if self.scoring_method=='UCell':
   if self.namespace!='CURRENT_UCELL_V1_'+self.signature_version.upper():raise ValueError('UCell namespace/version mismatch.')
   if self.parameters.get('max_rank')!=1500 or self.parameters.get('ties_method')!='average' or '0.8.0' not in self.scoring_version:raise ValueError('Current frozen UCell parameters/version not verified.')
  elif self.scoring_method=='LegacyMeanQ20':
   if not self.namespace.startswith('LEGACY_V4_') or self.signature_version!='original':raise ValueError('Legacy/current values cannot share a namespace.')
  else:raise NotImplementedError('AUCell/singscore/AddModuleScore adapters are interfaces only; no unverified implementation is supplied.')
  if self.species=='human' and (not self.mapping_hash or not self.mapping_version):raise ValueError('Human input requires the audited mapping version and hash.')
  return True
 def to_dict(self):return asdict(self)
@dataclass
class ScoreBundle:
 scores:pd.DataFrame
 metadata:pd.DataFrame
 spec:ScoreSpec
 coverage:pd.DataFrame
 detected_members:pd.DataFrame|None=None
 def __post_init__(self):
  self.spec.validate();s=self.scores.copy();m=self.metadata.copy()
  if 'cell_id' in s:s=s.set_index('cell_id')
  if 'cell_id' in m:m=m.set_index('cell_id')
  s.index=s.index.astype(str);m.index=m.index.astype(str)
  if s.index.has_duplicates or m.index.has_duplicates:raise ValueError('Duplicate cell IDs; observations cannot be counted twice.')
  if len(s)==0 or set(s.index)!=set(m.index):raise ValueError('Cell IDs in scores and metadata must match exactly.')
  self.missing_programs=[x for x in ORDER if x not in s]
  unknown=[c for c in s if c.startswith(('S0','S1','L0','A0')) and c not in ORDER]
  if unknown:raise ValueError('Unknown program identifiers: '+str(unknown))
  s=s.reindex(columns=ORDER).astype(np.float32);m=m.loc[s.index].copy()
  if self.spec.scoring_method=='UCell':
   vals=s.to_numpy();finite=vals[np.isfinite(vals)]
   if len(finite) and (finite.min() < -1e-6 or finite.max()>1+1e-6):raise ValueError('Declared UCell scores outside [0,1].')
  for c in ['sample_id','study_id','source_block_id','biological_unit_id','coarse_identity','tissue','condition','modality','batch']:
   if c not in m:m[c]='Unknown'
   m[c]=m[c].map(lambda v:str(v) if known(v) else 'Unknown').astype(str)
  if 'species' not in m:m['species']=self.spec.species
  if not m.species.eq(self.spec.species).all():raise ValueError('Mixed or inconsistent species in one score namespace.')
  if 'biological_unit_confirmed' not in m:m['biological_unit_confirmed']=False
  m['biological_unit_confirmed']=m.biological_unit_confirmed.map(confirmed).astype(bool)&m.biological_unit_id.map(known).astype(bool)
  for c in ['nUMI','nGenes']:
   if c not in m:m[c]=np.nan
   m[c]=pd.to_numeric(m[c],errors='coerce')
   if (m[c].dropna()<0).any():raise ValueError('Negative count QC covariate: '+c)
  for key,g in m.groupby(['study_id','source_block_id','sample_id'],dropna=False,observed=True):
   for col in ['biological_unit_id','condition','species','modality']:
    if g[col].nunique()>1:raise ValueError('One sample record has inconsistent '+col+'; provide a correct disaggregated sample mapping.')
  cov=self.coverage.copy()
  if 'feature_id' in cov:cov=cov.set_index('feature_id')
  if cov.index.has_duplicates:raise ValueError('Coverage must have one row per program in the selected version.')
  cov=cov.reindex(ORDER)
  if 'coverage' not in cov:cov['coverage']=np.nan
  cov['coverage']=pd.to_numeric(cov.coverage,errors='coerce')
  if ((cov.coverage.dropna()<0)|(cov.coverage.dropna()>1)).any():raise ValueError('Gene coverage must be a fraction in[0,1].')
  if 'scoring_eligible' not in cov:cov['scoring_eligible']=False
  cov['scoring_eligible']=cov.scoring_eligible.map(confirmed).astype(bool)
  if self.detected_members is not None:self.detected_members=self.detected_members.reindex(index=s.index,columns=ORDER)
  self.scores=s;self.metadata=m;self.coverage=cov
 def fingerprint(self):
  return digest(dict(cell_ids=list(self.scores.index),score_bytes=hashlib.sha256(self.scores.to_numpy().tobytes()).hexdigest(),spec=self.spec.to_dict()))

def from_official(adata,metadata,version='expanded'):
 """Read genuine official AnnData scores and their stored hashes/configuration."""
 import mumomac
 mm=adata.uns.get('mumomac_reference')
 if mm is None:raise ValueError('Official score provenance is missing.')
 conf={k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in mm['configuration'].items()}
 if version not in conf['versions']:raise ValueError('Requested score version was not computed.')
 values=np.asarray(adata.obsm['X_mumomac_reference_'+version])
 if hashlib.sha256(values.tobytes()).hexdigest()!=mm['score_hashes'][version]:raise ValueError('Official score matrix hash changed.')
 fs=list(mm['feature_order']);assert fs==list(ORDER)
 ref=mumomac.load_reference();species=conf['species'];mapping_hash='';mapping_version='mouse_native'
 if species=='human':
  from importlib import resources
  files=resources.files('mumomac.resources');mapping_hash=digest({n:hashlib.sha256(files.joinpath(n).read_bytes()).hexdigest() for n in ['human_la_one_to_one_v1.tsv','human_one_to_one.tsv','reference_version_manifest.json','signatures_manifest.json']});mapping_version='Frozen_v1_LA_global_one_to_one_plus_legacy_S16_Ensembl_mapping; legacy_S_database_release_missing'
 # Selecting Expanded from a dual-version run has the same algorithm parameters as an Expanded-only run.
 semantic_conf={**conf,'signature_version':version,'versions':[version]}
 spec=ScoreSpec('CURRENT_UCELL_V1_'+version.upper(),version,species,'UCell','Toolkit0.3.0/pyucell0.8.0', {f['feature_id']:f['signatures'][version]['gene_set_hash'] for f in ref['features']},parameters=semantic_conf,mapping_version=mapping_version,mapping_hash=mapping_hash)
 coverage=mm['gene_coverage'];coverage=coverage[coverage.signature_version==version]
 return ScoreBundle(pd.DataFrame(values,index=adata.obs_names,columns=ORDER),metadata,spec,coverage,pd.DataFrame(adata.obsm['X_mumomac_detected_'+version],index=adata.obs_names,columns=ORDER))

def score_provided(adata,metadata,*,species='mouse',version='expanded',chunk_size=500):
 """Only invokes unchanged formal scoring. No Atlas mapping or classifier dependency."""
 import mumomac
 scored=mumomac.score_signature_versions(adata,species=species,signature_version=version,layer='counts',expression_type='counts',coverage_policy='reference',max_rank=1500,chunk_size=chunk_size,copy=True)
 return from_official(scored,metadata,version)
