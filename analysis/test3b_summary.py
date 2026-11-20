import json, numpy as np
R=json.load(open('results/test3b.json'))
def pool(rs):
    n=sum(r['runs'] for r in rs); v=sum(r['viol'] for r in rs)
    return n,v,100*v/n,float(np.median([r['sat'] for r in rs])),float(np.mean([r['rms'] for r in rs])),float(np.mean([r['e2'] for r in rs])),float(np.mean([r['utv'] for r in rs])),max(r['excess'] for r in rs)
print('T torque ceilings; predicted 14.3, 20.2, 28.6')
G=sorted({r['lam'] for r in R if r['part']=='T'})
cs={}
for um in (4.,8.,16.):
    row=[(l,round(100*pool([r for r in R if r['part']=='T' and r['umax']==um and r['lam']==l])[3])) for l in G]
    print(um,row); c=[l for l,s in row if s>10]; cs[um]=c[0] if c else None
print(cs,'slope vs umax',np.polyfit(np.log(list(cs)),np.log(list(cs.values())),1)[0],'predicted 0.5')
print('N no cancel')
for l in sorted({r['lam'] for r in R if r['part']=='N'}):
    p=pool([r for r in R if r['part']=='N' and r['lam']==l]); print(l,p[1],'/',p[0],round(p[2],1),'sat',round(100*p[3],1),'rms',round(p[4],4),'e2',round(p[5],4),'utv',round(p[6]),'exc',round(p[7],3))
print('C clearance 5%')
for d in sorted({r['design'] for r in R if r['part']=='C'}):
    for l in (14.,20.,28.,40.,57.):
        p=pool([r for r in R if r['part']=='C' and r['design']==d and r['lam']==l]); print(d,l,p[1],'/',p[0],round(p[2],1),'sat',round(100*p[3],1),'rms',round(p[4],4),'e2',round(p[5],4),'utv',round(p[6]),'exc',round(p[7],3))
