import argparse,json
import pandas as pd
from mumomac_coupling.adapters import ScoreBundle,ScoreSpec
from .engine import run
p=argparse.ArgumentParser(description='MuMoMac experimental continuous-first representation from frozen scores only.')
for k in ['scores','metadata','coverage','spec','config','output']:p.add_argument('--'+k,required=True)
p.add_argument('--detected-members');p.add_argument('--frozen-model');p.add_argument('--mode',choices=['auto','force_gmm','force_tail','continuous_only','relative_high']);p.add_argument('--relative-quantile',type=float);p.add_argument('--joint-index',action='store_true');p.add_argument('--tail-profile',action='store_true');a=p.parse_args()
b=ScoreBundle(pd.read_parquet(a.scores),pd.read_csv(a.metadata,sep='\t'),ScoreSpec(**json.load(open(a.spec))),pd.read_csv(a.coverage,sep='\t'),pd.read_parquet(a.detected_members).set_index('cell_id') if a.detected_members else None)
r=run(b,a.config,a.output,model=json.load(open(a.frozen_model)) if a.frozen_model else None,mode=a.mode,relative={'quantile':a.relative_quantile} if a.relative_quantile is not None else None,include_joint=a.joint_index,include_tail=a.tail_profile)
print(json.dumps({'output':a.output,'n_cells':len(b.scores),'model_hash':r['model']['model_hash'],'J_soft_requested':a.joint_index,'tail_profile_requested':a.tail_profile,'parameter_hash':r['calls']['parameter_hash']}))
