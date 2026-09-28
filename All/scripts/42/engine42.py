"""Fixed-margin and structural-zero constrained site-interaction analysis.
Python standard library + numpy. Four independent constrained row-trade chains.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import sys,json,csv,math,time
from pathlib import Path
import numpy as np

OUT=Path(sys.argv[1]);DATA=json.loads((OUT/'42_inputs.json').read_text(encoding='utf-8'))
CHAINS=4;DRAWS=250;BURN_SWEEPS=100;THIN_SWEEPS=5
def arr(x):return x if isinstance(x,list) else [x]
def write(path,rows):
 if not rows:return
 with open(path,'w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def ploss(k,n,m):
 return np.array([0. if v>m else math.exp(math.lgamma(m+1)-math.lgamma(m-int(v)+1)-math.lgamma(n+1)+math.lgamma(n-int(v)+1)) for v in k])
def pair_metrics(M):
 r=M.sum(1);i,j=np.triu_indices(len(r),1);a=(M.astype(np.float32)@M.T.astype(np.float32))[i,j].astype(float)
 lo=np.minimum(r[i],r[j]);hi=np.maximum(r[i],r[j]);valid=lo>0;unequal=valid&(lo<hi)
 containment=np.divide(a,lo,out=np.full_like(a,np.nan),where=valid)
 sor=np.divide(r[i]+r[j]-2*a,r[i]+r[j],out=np.full_like(a,np.nan),where=(r[i]+r[j])>0)
 turn=1-containment;nest=sor-turn
 vals=np.array([np.nanmean(containment[unequal]),np.nanmean(turn[valid]),np.nanmean(nest[valid])])
 return vals,(i,j,lo,hi,a,containment,sor,turn,nest,valid,unequal)
def bh(p):
 p=np.asarray(p);o=np.argsort(p);z=p[o]*len(p)/np.arange(1,len(p)+1);z=np.minimum.accumulate(z[::-1])[::-1];ans=np.empty(len(p));ans[o]=np.minimum(1,z);return ans
def diagnostics(x):
 # Split R-hat and positive-sequence ESS for scalar chain traces; diagnostics
 # cannot certify irreducibility of a structural-zero state space.
 c,n=x.shape;half=n//2;y=np.concatenate([x[:,:half],x[:,-half:]],axis=0)
 W=np.mean(np.var(y,axis=1,ddof=1));B=half*np.var(y.mean(1),ddof=1)
 if W<1e-20:return float('nan'),0.
 rhat=math.sqrt(((half-1)/half*W+B/half)/W)
 centered=x-x.mean(1,keepdims=True);v=np.mean(np.sum(centered**2,axis=1))
 ac=[]
 for lag in range(1,min(100,n-1)):
  rho=np.mean(np.sum(centered[:,:-lag]*centered[:,lag:],axis=1))/v
  if rho<=0:break
  ac.append(rho)
 return rhat,min(c*n,c*n/(1+2*sum(ac)))
def trades(rows,allowed,nsteps,rng):
 moved=0;n=len(rows)
 for _ in range(nsteps):
  i=int(rng.integers(n));j=int(rng.integers(n-1));j+=j>=i
  ai=rows[i]-rows[j];aj=rows[j]-rows[i]
  mi=ai if allowed is None else ai&allowed[j]
  mj=aj if allowed is None else aj&allowed[i]
  if not mi or not mj:continue
  pool=sorted(mi|mj);rng.shuffle(pool);ni=set(pool[:len(mi)]);nj=set(pool[len(mi):])
  if ni!=mi:
   rows[i]=(rows[i]-mi)|ni;rows[j]=(rows[j]-mj)|nj;moved+=1
 return moved
def to_matrix(rows,L):
 M=np.zeros((len(rows),L),dtype=np.uint8)
 for i,row in enumerate(rows):M[i,list(row)]=1
 return M

# Known cases: strict nestedness, disjointness and fixed-support mean retention.
checks=[]
def check(name,ok):
 checks.append(dict(check=name,pass_check=bool(ok)))
 if not ok:raise AssertionError(name)
v,_=pair_metrics(np.array([[1,0,0],[1,1,0],[1,1,1]],dtype=np.uint8));check('nested containment = 1',abs(v[0]-1)<1e-12);check('nested turnover = 0',abs(v[1])<1e-12)
v,_=pair_metrics(np.array([[1,0,0],[0,1,1]],dtype=np.uint8));check('disjoint containment = 0',abs(v[0])<1e-12);check('disjoint turnover = 1',abs(v[1]-1)<1e-12)
for n in range(3,7):
 import itertools
 for m in range(1,n):
  for k in range(1,n+1):
   exact=np.mean([set(range(k)).issubset(z) for z in itertools.combinations(range(n),m)])
   check(f'exact_loss_{n}_{m}_{k}',abs(exact-ploss([k],n,m)[0])<1e-10)
summary=[];global_rows=[];all_sites=[];diag_rows=[];null_records=[];removal=[]
for di,(name,d) in enumerate(DATA.items()):
 start=time.time();folder=OUT/name;folder.mkdir(exist_ok=True)
 n=len(d['sites']);L=len(d['observed']);M=np.zeros((n,L),dtype=np.uint8);A=M.copy()
 for j,ss in enumerate(d['observed']):M[arr(ss),j]=1
 for j,ss in enumerate(d['allowed']):A[arr(ss),j]=1
 check(name+' valid observed opportunities',np.all(M<=A))
 K=M.sum(0).astype(int);rich=M.sum(1).astype(int);obs,pair=pair_metrics(M)
 i,j,lo,hi,a,cont,sor,turn,nest,valid,uneq=pair
 write(folder/'42_pairwise_site_metrics.csv',[dict(site1=d['sites'][x],site2=d['sites'][y],richness1=int(rich[x]),richness2=int(rich[y]),shared=int(a[t]),containment=cont[t],sorensen=sor[t],turnover=turn[t],nestedness_resultant=nest[t],included_containment=bool(uneq[t])) for t,(x,y) in enumerate(zip(i,j))])
 summary.append(dict(dataset=name,sites=n,links=L,occurrences=int(M.sum()),nonempty_pairs=int(valid.sum()),unequal_nonempty_pairs=int(uneq.sum()),containment=obs[0],turnover=obs[1],nestedness_resultant=obs[2],sorensen=float(np.nanmean(sor[valid]))))
 # Exact focal-site excess link loss at approximately half removal.
 m=max(1,min(n-1,round(n*.5)));base=ploss(K,n,m).mean();p0=ploss(K,n-1,m-1);delta=ploss(K-1,n-1,m-1)-p0
 importance=M@delta/L+p0.mean()-base;unique=M@(K==1).astype(float)
 site_samples={};snapshots={};models=('Margins','Opportunity')
 for model in models:
  allowed=None if model=='Margins' else [set(np.flatnonzero(A[r])) for r in range(n)]
  traces=np.empty((CHAINS,DRAWS,3));scores=[];snaps=[];movement=[]
  for chain in range(CHAINS):
   rng=np.random.default_rng(42000+di*100+chain*2+(model=='Opportunity'))
   rows=[set(np.flatnonzero(M[r])) for r in range(n)]
   trades(rows,allowed,BURN_SWEEPS*n,rng)
   changes=0;prev=M.copy();changed_fr=[]
   for draw in range(DRAWS):
    changes+=trades(rows,allowed,THIN_SWEEPS*n,rng);B=to_matrix(rows,L)
    check_ok=np.array_equal(B.sum(0),K) and np.array_equal(B.sum(1),rich) and (model=='Margins' or np.all(B<=A))
    if not check_ok:raise AssertionError('Margins/opportunity changed')
    vals,_=pair_metrics(B);traces[chain,draw]=vals;scores.append(B@delta/L+p0.mean()-base)
    changed_fr.append(float(np.sum((B!=M)&(M==1))/M.sum()))
    if draw%25==0:snaps.append(B.copy())
    null_records.append(dict(dataset=name,model=model,chain=chain,draw=draw,containment=vals[0],turnover=vals[1],nestedness_resultant=vals[2]))
   movement.append(changes/(DRAWS*THIN_SWEEPS*n))
   diag_rows.append(dict(dataset=name,model=model,chain=chain,changed_trade_fraction=movement[-1],mean_observed_edges_relocated=np.mean(changed_fr)))
  site_samples[model]=np.array(scores);snapshots[model]=snaps
  for t,metric in enumerate(('containment','turnover','nestedness_resultant')):
   x=traces[:,:,t];flat=x.ravel();rh,ess=diagnostics(x);sd=flat.std(ddof=1)
   degenerate=sd<1e-12
   reliable=not degenerate and np.isfinite(rh) and rh<1.05 and ess>=100 and min(movement)>0
   p=min(1.,2*min((1+np.sum(flat>=obs[t]-1e-12))/(len(flat)+1),(1+np.sum(flat<=obs[t]+1e-12))/(len(flat)+1)))
   global_rows.append(dict(dataset=name,model=model,metric=metric,observed=obs[t],null_mean=flat.mean(),null_low=np.quantile(flat,.025),null_high=np.quantile(flat,.975),ses=(obs[t]-flat.mean())/sd if not degenerate else float('nan'),p_approx=p if reliable else float('nan'),rhat=rh,ess=ess,diagnostics_ok=reliable,degenerate=degenerate))
  print(name,model,'draws',CHAINS*DRAWS,'trade fraction',round(np.mean(movement),3),flush=True)
 check(name+' all null draws preserve both margins and allowed cells',True)
 check(name+' focal-site importance averages to zero',abs(importance.mean())<1e-10)
 for s in range(n):
  for model in models:
   x=site_samples[model][:,s];sd=x.std(ddof=1);rh,ess=diagnostics(x.reshape(CHAINS,DRAWS));reliable=sd>1e-12 and np.isfinite(rh) and rh<1.05 and ess>=100
   p=min(1.,2*min((1+sum(x>=importance[s]-1e-12))/(len(x)+1),(1+sum(x<=importance[s]+1e-12))/(len(x)+1)))
   all_sites.append(dict(dataset=name,site=d['sites'][s],model=model,richness=int(rich[s]),unique_links=int(unique[s]),removal_fraction=m/n,importance=importance[s],null_mean=x.mean(),null_low=np.quantile(x,.025),null_high=np.quantile(x,.975),adjusted_importance=importance[s]-x.mean(),p_approx=p if reliable else float('nan'),rhat=rh,ess=ess,diagnostics_ok=reliable))
 # Secondary removal diagnostics are computed for every dataset, but the
 # README limits substantive interpretation to supported, well-mixed results.
 fractions=[.2,.4,.6,.8];rng=np.random.default_rng(8000+di)
 matrices=[('Observed',M)]+[(mod,B) for mod in models for B in snapshots[mod]]
 for mod,B in matrices:
  # Repeated random tie-breaking for equal site richness; same row margins.
  for direction in ('Rich first','Poor first','Random'):
   curves=[]
   for rep in range(100 if mod=='Observed' else 10):
    perm=rng.permutation(n)
    if direction!='Random':perm=perm[np.argsort((-rich if direction=='Rich first' else rich)[perm],kind='stable')]
    cumulative=np.cumsum(B[perm],axis=0)
    curves.append([float(np.mean(K-cumulative[max(1,min(n-1,round(n*f)))-1]>0)) for f in fractions])
   for t,f in enumerate(fractions):
    x=np.array(curves)[:,t]
    removal.append(dict(dataset=name,model=mod,order=direction,removed=max(1,min(n-1,round(n*f))),fraction=max(1,min(n-1,round(n*f)))/n,retention_mean=x.mean(),retention_variance=x.var(ddof=1),exact_random_mean=float(np.mean(1-ploss(K,n,max(1,min(n-1,round(n*f))))))))
 # Similarity heatmap order is supplied by the R plotting stage.
 np.savez_compressed(folder/'42_matrices.npz',observed=M,allowed=A)
 write(folder/'42_incidence.csv',[dict(site=d['sites'][s],consumer=arr(d['consumer'])[j],resource=arr(d['resource'])[j],support=int(K[j])) for s,j in zip(*np.where(M))])
 print(name,'complete',round(time.time()-start,1),'seconds',flush=True)

# Global BH correction by metric across 10 datasets x 2 models, reliable only.
for metric in ('containment','turnover','nestedness_resultant'):
 ids=[i for i,r in enumerate(global_rows) if r['metric']==metric and np.isfinite(r['p_approx'])]
 for i,q in zip(ids,bh([global_rows[i]['p_approx'] for i in ids])):global_rows[i]['q_bh']=q
for r in global_rows:r.setdefault('q_bh',float('nan'))
# Site BH families are dataset x model. No pairwise pseudo-replication tests.
for name in DATA:
 for model in ('Margins','Opportunity'):
  ids=[i for i,r in enumerate(all_sites) if r['dataset']==name and r['model']==model and np.isfinite(r['p_approx'])]
  for i,q in zip(ids,bh([all_sites[i]['p_approx'] for i in ids])):all_sites[i]['q_bh']=q
for r in all_sites:r.setdefault('q_bh',float('nan'))
for name in DATA:
 write(OUT/name/'42_null_comparisons.csv',[r for r in global_rows if r['dataset']==name])
 write(OUT/name/'42_site_importance.csv',[r for r in all_sites if r['dataset']==name])
write(OUT/'42_dataset_summary.csv',summary);write(OUT/'42_null_comparisons.csv',global_rows)
write(OUT/'42_site_importance.csv',all_sites);write(OUT/'42_chain_movement.csv',diag_rows)
write(OUT/'42_null_traces.csv',null_records);write(OUT/'42_removal_diagnostics.csv',removal)
write(OUT/'42_validation_checks.csv',checks)
(OUT/'42_settings.json').write_text(json.dumps(dict(chains=CHAINS,draws_per_chain=DRAWS,burn_sweeps=BURN_SWEEPS,thinning_sweeps=THIN_SWEEPS,seed_base=42000)),encoding='utf-8')
print('All datasets completed',flush=True)
