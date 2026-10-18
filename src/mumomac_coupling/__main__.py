import argparse,json,hashlib
from pathlib import Path
import pandas as pd
from .adapters import ScoreSpec,ScoreBundle
from .engine import run_engine
def main():
 p=argparse.ArgumentParser(prog='python -m mumomac_coupling',description='Experimental, provenance-gated sample-stratified program coupling; no Atlas dependency.')
 for name in ['scores','metadata','score-spec','coverage','config','output']:p.add_argument('--'+name,required=True)
 p.add_argument('--frozen-thresholds');a=p.parse_args()
 def read(path):return pd.read_parquet(path) if str(path).endswith('.parquet') else pd.read_csv(path,sep='\t',keep_default_na=False)
 s=read(a.scores);m=read(a.metadata);spec=ScoreSpec(**json.loads(Path(a.score_spec).read_text()));bundle=ScoreBundle(s,m,spec,read(a.coverage));cfg=json.loads(Path(a.config).read_text());model=json.loads(Path(a.frozen_thresholds).read_text()) if a.frozen_thresholds else None
 cfg['input_manifest']=[dict(path=str(Path(p).resolve()),sha256=hashlib.sha256(Path(p).read_bytes()).hexdigest()) for p in [a.scores,a.metadata,a.score_spec,a.coverage,a.config]+([a.frozen_thresholds] if a.frozen_thresholds else [])]
 r=run_engine(bundle,cfg,a.output,model);print(json.dumps(dict(output=a.output,n_samples=len(r['samples']),threshold_mode=r['thresholds']['mode'],resource=r['resource']),indent=2))
if __name__=='__main__':main()
