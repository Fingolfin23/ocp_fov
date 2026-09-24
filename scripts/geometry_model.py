"""Smooth illustrative road; all events and hits are computed, not drawn by hand.

The centerline is a radius-10 circular arc parameterized by arc length s.
Positive, Gaussian-modulated widths preserve an injective radial Frenet map.
This is a mathematical illustration, not a run of the original cache pipeline.
"""
import numpy as np
from scipy.optimize import brentq

R=10.; H=18.
bumpsL=[(1.7,3.8,.85),(1.3,9.,1.)]
bumpsR=[(1.2,6.5,1.2),(1.6,15.,1.3)]

def wd(s,side):
    s=np.asarray(s,dtype=float)
    w=np.ones_like(s)*.8; dw=np.zeros_like(s)
    for amp,mu,sigma in (bumpsL if side=='L' else bumpsR):
        b=amp*np.exp(-.5*((s-mu)/sigma)**2)
        w+=b;dw+=b*(-(s-mu)/sigma**2)
    return w,dw

def B(s,side,derivative=False):
    s=np.asarray(s,dtype=float); w,dw=wd(s,side); sign=-1 if side=='L' else 1
    r=R+sign*w; rp=sign*dw;a=s/R
    if derivative:
        return np.stack([rp*np.sin(a)+r/R*np.cos(a),-rp*np.cos(a)+r/R*np.sin(a)],axis=-1)
    return np.stack([r*np.sin(a),R-r*np.cos(a)],axis=-1)

def cross(a,b):return a[...,0]*b[...,1]-a[...,1]*b[...,0]

def roots(f,lo,hi,n=10001):
    s=np.linspace(lo,hi,n);vals=f(s);out=[]
    for i in np.where(vals[:-1]*vals[1:]<0)[0]:
        out.append(brentq(f,s[i],s[i+1],xtol=1e-13))
    return out

def events():
    out=[]
    for side in ['L','R']:
        f=lambda s:cross(B(s,side),B(s,side,True))
        for st in roots(f,.001,H):
            delta=(f(st+.0001)-f(st-.0001))/.0002
            if (side=='L' and delta<0) or (side=='R' and delta>0):continue
            T=B(st,side);d=T/np.linalg.norm(T)
            if np.dot(d,[np.cos(st/R),np.sin(st/R)])<=0:continue
            hits={}
            for hit_side in ['L','R']:
                hs=roots(lambda s:cross(B(s,hit_side),T),st+.0001,H)
                hs=[s for s in hs if np.dot(B(s,hit_side)-T,d)>1e-7]
                hits[hit_side]=hs
            first=[(v[0],k) for k,v in hits.items() if v]
            first.sort()
            out.append(dict(side=side,st=st,T=T,hits=hits,first=first[0] if first else None))
    return sorted(out,key=lambda d:d['st'])

if __name__=='__main__':
    for ev in events():print(ev)
