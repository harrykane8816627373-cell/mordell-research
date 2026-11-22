from fractions import Fraction as F
from math import factorial
from pathlib import Path
from functools import lru_cache
import sys

S=2**20
def need(b,msg):
    if not b: raise RuntimeError(msg)

def atan_pair(b):
    x=F(1,b)
    part=sum((-1)**k*x**(2*k+1)/(2*k+1) for k in range(12))
    return part,part+x**25/25

def pi_bounds():
    a,b=atan_pair(5);c,d=atan_pair(239)
    p,q=F(31415926,10000000),F(31415927,10000000)
    need(p<16*a-4*d<16*b-4*c<q,'Machin pi enclosure')
    return p,q

def make_phases(n):
    p,q=pi_bounds()
    need(n%16==0,'phase subdivision compatible with quarters')
    @lru_cache(None)
    def sine(j):
        j=j%n
        if j>n//2:
            l,h=sine(n-j);return -h,-l
        if j>n//4:return sine(n//2-j)
        if j==0:return 0,0
        if j==n//4:return S,S
        x=2*p*j/n
        v=sum((-1)**k*x**(2*k+1)/factorial(2*k+1) for k in range(8))
        hi=v+x**17/factorial(17)+2*(q-p)*j/n
        return int(v*S//1),int(-((-hi*S)//1))
    def complex_enclosure(i,k):
        cells=range(i*k,(i+1)*k+1)
        cosines=[sine(t+n//4) for t in cells]
        sines=[sine(t) for t in cells]
        return [min(l for l,h in cosines),max(h for l,h in cosines),
                -max(h for l,h in sines),-min(l for l,h in sines)]
    lower=max(i for i in range(n//4) if F(sine(i+n//4)[0],S)>F(78,100))
    upper=min(i for i in range(n//4) if F(sine(i+n//4)[1],S)<F(59,100))
    return '\n'.join([f'{n} {lower} {upper}']+[' '.join(str(v) for k in (1,2,4) for v in complex_enclosure(i,k)) for i in range(n)])
if __name__=='__main__':print(make_phases(int(sys.argv[1])))
