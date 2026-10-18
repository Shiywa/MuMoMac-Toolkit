"""Frozen definitions. No signature selection or refitting at runtime."""
from dataclasses import dataclass
from functools import lru_cache
from importlib import resources
import hashlib,json
import pandas as pd

PROGRAM_ORDER=tuple(f'S{i:02d}' for i in range(1,17))
@dataclass(frozen=True)
class Program:
    program_id: str
    name: str
    genes: tuple[str,...]
    original_genes: tuple[str,...]
    added_genes: tuple[str,...]
    mouse_genes: tuple[str,...]
    unmapped_genes: tuple[str,...]
    species: str
    signature_version: str
    definition_version: str
    description: str
    evidence_status: str
    limitations: str

@lru_cache(maxsize=1)
def _bundle():
    root=resources.files('mumomac.resources')
    m=json.loads(root.joinpath('signatures_manifest.json').read_text())
    for name,digest in m['resource_sha256'].items():
        if hashlib.sha256(root.joinpath(name).read_bytes()).hexdigest()!=digest:
            raise RuntimeError(f'Frozen MuMoMac resource integrity failure: {name}')
    if tuple(m['program_order'])!=PROGRAM_ORDER:
        raise RuntimeError('Frozen program order is invalid.')
    with root.joinpath('human_one_to_one.tsv').open('r') as f:
        h=pd.read_csv(f,sep='\t')
    if h.groupby('mouse_gene').human_gene.nunique().max()>1 or h.groupby('human_gene').mouse_gene.nunique().max()>1:
        raise RuntimeError('Ambiguous frozen human mapping.')
    return m,dict(zip(h.mouse_gene,h.human_gene))

def signature_manifest():
    """Return an isolated copy of versions, provenance, evidence and resource hashes."""
    return json.loads(json.dumps(_bundle()[0]))

def load_signatures(signature_version='original',species='mouse'):
    """Return 16 immutable Programs; human is a declared experimental conversion."""
    if signature_version not in ('original','expanded'):
        raise ValueError("signature_version must be 'original' or 'expanded'.")
    if species not in ('mouse','human'):
        raise ValueError("species must be 'mouse' or 'human'.")
    m,h=_bundle();out={}
    for pid in PROGRAM_ORDER:
        x=m['programs'][pid];original=tuple(x['original_genes']);full=tuple(x['expanded_genes']);selected=original if signature_version=='original' else full
        convert=lambda gs:tuple(h[g] for g in gs if g in h) if species=='human' else tuple(gs)
        additions=tuple(g for g in full if g not in original) if signature_version=='expanded' else ()
        out[pid]=Program(pid,x['program_name'],convert(selected),convert(original),convert(additions),selected,tuple(g for g in selected if g not in h) if species=='human' else (),species,signature_version,m['definition_versions'][signature_version],x['biological_description'],x['evidence_status'],x['known_limitations']+' '+x['external_evidence'])
    return out
