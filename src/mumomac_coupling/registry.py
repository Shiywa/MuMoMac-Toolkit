import json,hashlib,itertools,math
from pathlib import Path
from importlib import resources
ORDER=tuple([f'S{i:02}' for i in range(1,17)]+[f'L{i:02}' for i in range(1,6)]+['A01'])
VERSION='0.1.0-experimental'
def resource_json(name):
 p=resources.files('mumomac_coupling.resources').joinpath(name)
 return json.loads(p.read_text())
def registry():return resource_json('FROZEN_REFERENCE_V1.json')
def signatures(species,version):return resource_json(f'MAPPED_SIGNATURES_{species}_{version}.json')
def canonicalize(x):
 if isinstance(x,dict):return {str(k):canonicalize(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [canonicalize(v) for v in x]
 if hasattr(x,'tolist'):return canonicalize(x.tolist())
 if isinstance(x,float) and not math.isfinite(x):return None
 if x is None or isinstance(x,(str,int,float,bool)):return x
 return str(x)
def digest(x):return hashlib.sha256(json.dumps(canonicalize(x),sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def pairs(optional_lifecycle=(),only=None):
 r=registry();hist={tuple(v):k for k,v in r['combinations'].items()};out=[]
 allowed=set(ORDER[:16])|set(optional_lifecycle)
 if 'A01' in optional_lifecycle or any(x not in ORDER[16:21] for x in optional_lifecycle):raise ValueError('Only explicitly applicable L01–L05 can be added; A01 is continuous.')
 candidates=list(itertools.combinations(ORDER[:16],2))
 # Optional lifecycle pairs are explicit historical/tutorial pairs, not an automatic22-feature binary screen.
 for ab in r['combinations'].values():
  if set(ab)<=allowed and tuple(ab) not in candidates:candidates.append(tuple(ab))
 for a,b in candidates:
  pid=hist.get((a,b),a+'__'+b)
  if only is None or pid in only:out.append(dict(pair_id=pid,program_A=a,program_B=b,historical_candidate=pid in r['combinations']))
 if only is not None and set(only)-{x['pair_id'] for x in out}:raise ValueError('Requested pair absent or lifecycle applicability not declared.')
 return out
