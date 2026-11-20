import json, numpy as np
from collections import defaultdict
R=json.load(open('results/test3.json'))
print(len(R), 'all finite', all(r['finite'] for r in R))
def pool(rs):
    n=sum(r['runs'] for r in rs); v=sum(r['viol'] for r in rs)
    return dict(n=n,v=v,pct=100*v/n,est=sum(r['est'] for r in rs),sat=float(np.median([r['sat'] for r in rs])),rms=float(np.mean([r['rms'] for r in rs])),e2=float(np.mean([r['e2'] for r in rs])),utv=float(np.mean([r['utv'] for r in rs])),exc=max(r['excess'] for r in rs),p2p=float(np.mean([r['p2p'] for r in rs])))
D=["base","low-pass","dead zone 1","dead zone 2","linear"]
G=sorted({r['lam'] for r in R if r['part']=='S2'})
for part in ('S1','S2'):
    print('==',part)
    for d in D:
        print(d)
        for lam in G:
            p=pool([r for r in R if r['part']==part and r['design']==d and r['lam']==lam])
            print(f"  {lam:5.0f} viol {p['v']:3d}/{p['n']} {p['pct']:5.1f}% est {p['est']:3d} sat {100*p['sat']:5.1f}% rms {p['rms']:.4f} e2 {p['e2']:.4f} utv {p['utv']:7.0f} exc {p['exc']:.3f}")
print('== R1')
for dz in ("dead zone 1","dead zone 2"):
    for lam in (40.,57.):
        p=pool([r for r in R if r['part']=='S2' and r['design']==dz and r['lam']==lam])
        best=min(pool([r for r in R if r['part']=='S2' and r['design']==d and r['lam']==lam])['rms'] for d in D)
        print(dz,lam,round(p['pct'],2),round(100*p['sat'],1),round(p['rms'],4),'best',round(best,4))
print('== R2 selection on seed 7, eval 8,9')
sel={}
for d in D:
    c=[]
    for lam in G:
        r=[x for x in R if x['part']=='S2' and x['design']==d and x['lam']==lam and x['seed']==7][0]
        if r['viol']/r['runs']<=0.01 and r['sat']<0.10: c.append((r['e2'],lam))
    sel[d]=min(c)[1] if c else None
    if c:
        p=pool([x for x in R if x['part']=='S2' and x['design']==d and x['lam']==sel[d] and x['seed'] in (8,9)])
        print(d,'gain',sel[d],'eval',p)
    else: print(d,'no eligible gain')
print('== S3')
for d in D:
    for lam in (28.,40.,57.):
        p=pool([r for r in R if r['part']=='S3' and r['design']==d and r['lam']==lam])
        print(d,lam,'viol',p['v'],'/',p['n'],'sat',round(100*p['sat'],1),'rms',round(p['rms'],4),'e2',round(p['e2'],4))
print('== small (count=',2*np.pi/4096,')')
for sc in ('slow','hold'):
    for d in D:
        for lam in (28.,40.,57.):
            p=pool([r for r in R if r['part']=='small' and r['scen']==sc and r['design']==d and r['lam']==lam])
            print(sc,d,lam,'rms',round(p['rms'],5),'p2p counts',round(p['p2p']/(2*np.pi/4096),2),'e2',round(p['e2'],4),'utv',round(p['utv']),'sat',round(100*p['sat'],1),'viol',p['v'])
print('== R4')
cs={}
for b in (10,12,14,16):
    row=[]
    for lam in G:
        rs=[r for r in R if r['design']=='base' and r['lam']==lam and r['bits']==b and r['sensor']=='S2' and r['scen']=='std']
        row.append((lam,round(100*float(np.median([x['sat'] for x in rs])))))
    print(b,row)
    c=[l for l,s in row if s>10]; cs[b]=c[0] if c else None
print(cs)
x=[np.log(2*np.pi/2**b) for b,c in cs.items() if c]; y=[np.log(c) for c in cs.values() if c]
print('slope',np.polyfit(x,y,1)[0] if len(x)>=2 else None,len(x))
