#!/usr/bin/env python3
"""Cross-validated error model for ComPair.sh calls on the read-labelled set (tests/build_training_set.py).
Logistic regression (class-balanced, L2) on the ComPair metrics + call type, held out by dvl scaffold (5 folds).
Usage: fit_error_model.py TRAIN.tsv   (needs numpy)"""
import sys
import csv, numpy as np, zlib, collections
rows=list(csv.DictReader(open(sys.argv[1]),delimiter='\t'))
M=['LF','FL','LFcp','RFlength','FR','RFcp','OneTwo','OneSINE','TwoSINE','SL','SO','ST']
def num(x):
    try: return float(x)
    except: return np.nan
X=np.array([[num(r[m]) for m in M]+[r['status']=='PM', r['status']=='MP', r['status']=='SINE'] for r in rows],float)
X=np.where(np.isnan(X), np.nanmedian(X,0), X)
y=np.array([1-int(r['label']) for r in rows])          # 1 = call contradicted by reads
fold=np.array([zlib.crc32(r['dvl_scaffold'].encode())%5 for r in rows])
def fit(X,y,l2=1.0,it=300):
    mu,sd=X.mean(0),X.std(0)+1e-9; Z=(X-mu)/sd; Z=np.c_[np.ones(len(Z)),Z]; w=np.zeros(Z.shape[1])
    pos=y.mean(); sw=np.where(y==1,0.5/pos,0.5/(1-pos))   # class-balanced weights
    for _ in range(it):
        p=1/(1+np.exp(-Z@w)); g=Z.T@(sw*(p-y))/len(y)+l2*np.r_[0,w[1:]]/len(y)
        H=(Z.T*(sw*p*(1-p)))@Z/len(y)+np.diag(np.r_[0,np.ones(len(w)-1)])*l2/len(y)
        w-=np.linalg.solve(H,g)
    return mu,sd,w
def pred(m,X):
    mu,sd,w=m; Z=np.c_[np.ones(len(X)),(X-mu)/sd]; return 1/(1+np.exp(-Z@w))
def auc(s,y):
    o=np.argsort(s); r=np.empty(len(s)); r[o]=np.arange(1,len(s)+1)
    n1=y.sum(); n0=len(y)-n1; return (r[y==1].sum()-n1*(n1+1)/2)/(n1*n0)
scores=np.zeros(len(y))
for f in range(5):
    m=fit(X[fold!=f],y[fold!=f]); scores[fold==f]=pred(m,X[fold==f])
print(f'labelled {len(y)}, contradicted {y.sum()} ({100*y.mean():.2f}%)')
print(f'held-out-by-scaffold AUC (all calls) = {auc(scores,y):.3f}')
for st in ('SINE','PM','MP'):
    k=np.array([r['status']==st for r in rows]); print(f'  {st}: AUC {auc(scores[k],y[k]):.3f}  n={k.sum()} contradicted={y[k].sum()}')
# single-feature AUCs (lower value = more suspicious -> use -x)
for j,mn in enumerate(M):
    a=auc(-X[:,j],y); print(f'  feature {mn:9s} AUC(low = suspicious) {a:.3f}')
# precision gain: contradicted fraction in top 1%/5% risk
for q in (0.01,0.05,0.2):
    t=np.quantile(scores,1-q); k=scores>=t; print(f'  top {int(q*100)}% risk: contradicted {100*y[k].mean():.2f}% vs {100*y.mean():.2f}% overall')
m=fit(X,y); print('coefficients (standardized):', dict(zip(['icpt']+M+['isPM','isMP','isSINE'], np.round(m[2],2))))
