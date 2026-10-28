"""Independent exact spot audit of the separate rational-ascent certificate."""
from fractions import Fraction as F
from pathlib import Path
import json
import rational_ascent as a


def require(v,name):
    if not v:raise ArithmeticError(name)


def run():
    cases=0
    for d in (1,3):
        for n in range(3,10):
            points=[(i,(i*i+3*i)%7) for i in range(n)]
            row=a.ascent_step(points,d)
            if row['stationary']:continue
            v=row['points'];r=row['residual'];Q=row['Q'];R=row['R'];N=row['N'];t=row['step']
            require(sum(z[0] for z in r)==0 and sum(z[1] for z in r)==0,'zero-sum residual')
            pairing=(F(0),F(0))
            for z,w in zip(v,r):pairing=a.add(pairing,a.mul(a.conjugate(z),w,d))
            require(pairing==(0,0),'complex tangent orthogonality')
            require(sum(a.norm(z,d) for z in row['new_points'])==Q+t*t*R,'exact new norm')
            require(row['C']**2<=(d+1)*row['A'],'d-dependent rational modulus cap')
            ratio=row['after']/row['before'];y=(ratio-1)/(ratio+1)
            lower=2*sum((y**(2*j+1)/F(2*j+1) for j in range(20)),F(0))
            require(lower>=row['log_gain_lower'],'independent positive atanh log lower')
            require(row['log_gain_lower']==R*R/(row['H']+2*row['C']*R),'rational gain formula')
            cases+=1
    # Residual zero need not mean global maximality.
    line=a.ascent_step([(-1,0),(0,0),(1,0)])
    require(line['stationary'] and line['after']==F(27,2),'suboptimal stationary line')
    square=a.ascent_step([(0,0),(1,0),(0,1),(-1,0),(0,-1)])
    require(square['stationary'] and square['after']<5**5,'suboptimal stationary center-square')
    triangle=a.ascent_step([(0,0),(1,0),(F(1,2),F(1,2))],d=3)
    require(triangle['stationary'] and triangle['after']==27,'equilateral d3 control')
    return {'passed':True,'independent_nonstationary_steps':cases,
            'stationary_controls':3,'arithmetic':'Fraction; positive atanh lower series',
            'scope':'Independent spot and identity audit of the per-step certificate; no global-optimality inference.'}


if __name__=='__main__':
    report=run();(Path(__file__).parent/'ascent_phase_independent_checks.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))
