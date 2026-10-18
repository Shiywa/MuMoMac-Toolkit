"""Versioned 22-feature registry. Inherited S16 definitions are untouched."""
from importlib import resources
import copy, hashlib, json
import pandas as pd
from .signatures import load_signatures, _bundle

REFERENCE_ORDER=tuple([f'S{i:02}' for i in range(1,17)]+[f'L{i:02}' for i in range(1,6)]+['A01'])

def load_reference():
    root=resources.files('mumomac.resources')
    manifest=json.loads(root.joinpath('reference_version_manifest.json').read_text())
    for name,digest in manifest['resource_sha256'].items():
        if hashlib.sha256(root.joinpath(name).read_bytes()).hexdigest()!=digest:
            raise RuntimeError(f'Reference resource hash mismatch: {name}')
    registry=json.loads(root.joinpath('reference_registry_v1.json').read_text())
    if tuple(f['feature_id'] for f in registry['features'])!=REFERENCE_ORDER:
        raise RuntimeError('Invalid reference feature order')
    for f in registry['features']:
        for s in f['signatures'].values():
            genes=s['genes']
            digest=hashlib.sha256(json.dumps(sorted(genes),separators=(',',':')).encode()).hexdigest()
            if len(genes)!=len(set(genes)) or not genes or digest!=s['gene_set_hash']:
                raise RuntimeError(f'Invalid frozen signature {f["feature_id"]}')
        if not set(f['signatures']['original']['genes'])<=set(f['signatures']['expanded']['genes']):
            raise RuntimeError('Expanded lacks frozen Seed members')
    return copy.deepcopy(registry)

def resolve_feature(value):
    for f in load_reference()['features']:
        if value in [f['feature_id'],f['display_name'],f['canonical_name']]+f['aliases']:
            return f
    raise KeyError(f'Unknown reference feature {value!r}')

def reference_signatures(version='original',species='mouse'):
    if version not in ('original','expanded'):raise ValueError('Choose original or expanded')
    if species not in ('mouse','human'):raise ValueError('Species must be mouse or human')
    root=resources.files('mumomac.resources')
    with root.joinpath('human_la_one_to_one_v1.tsv').open() as handle: h=pd.read_csv(handle,sep='\t')
    if h.mouse_gene.duplicated().any() or h.human_gene.duplicated().any():raise RuntimeError('Ambiguous one-to-one mapping')
    la=dict(zip(h.mouse_gene,h.human_gene));smap=_bundle()[1];out={}
    for f in load_reference()['features']:
        s=copy.deepcopy(f['signatures'][version]);pid=f['feature_id'];mapper=smap if pid.startswith('S') else la
        s['mouse_genes']=list(s['genes'])
        s['unmapped_genes']=[g for g in s['genes'] if g not in mapper] if species=='human' else []
        if species=='human':
            for key in ['genes','original_genes','added_genes']:s[key]=[mapper[g] for g in s[key] if g in mapper]
        s.update(display_name=f['display_name'],feature_class=f['feature_class'],species=species,
                 construct_validity=f['construct_validity'],default_scope=f['default_scope'],
                 biological_limits=f['biological_limitations'],reference_version='1.0.0')
        out[pid]=s
    return out

def validate_reference():
    r=load_reference()
    for version in ('original','expanded'):
        old=load_signatures(version)
        for f in r['features'][:16]:
            if tuple(f['signatures'][version]['genes'])!=old[f['feature_id']].mouse_genes:
                raise RuntimeError('Inherited S16 member order changed')
    return {'reference_version':r['reference_version'],'n_features':22,'n_signatures':44,'hashes':'PASS','legacy_s16':'PASS'}
