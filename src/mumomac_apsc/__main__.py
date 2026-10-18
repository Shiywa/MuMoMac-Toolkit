import argparse,json
import pandas as pd
from mumomac_coupling.adapters import ScoreBundle,ScoreSpec
from .engine import run
p=argparse.ArgumentParser(description='APSC state extension for frozen MuMoMac scores; requires frozen C1 on PYTHONPATH.')
for k in ['scores','metadata','coverage','spec','config','output']:p.add_argument('--'+k,required=True)
p.add_argument('--detected-members');p.add_argument('--frozen-calibration')
a=p.parse_args();s=pd.read_parquet(a.scores);m=pd.read_csv(a.metadata,sep='\t');cov=pd.read_csv(a.coverage,sep='\t')
b=ScoreBundle(s,m,ScoreSpec(**json.load(open(a.spec))),cov,pd.read_parquet(a.detected_members).set_index('cell_id') if a.detected_members else None)
r=run(b,a.config,a.output,model=json.load(open(a.frozen_calibration)) if a.frozen_calibration else None)
print(json.dumps({'n_cells':len(b.scores),'parameter_hash':r['model']['parameter_hash'],'output':a.output}))
