"""Shrinking-bulk interpolation gives an effective whole-background cover."""
from fractions import Fraction as F
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'joint-bounds'))
sys.path.insert(0,str(ROOT/'distribution-geometry'))
from subset_moment import log_interval
from microscopic_covering import sqrt_interval


def boundary_certificate(n,scale=None):
    if type(n)is not int or n<64:
        raise ValueError('integer n>=64 required')
    if scale is None:scale=max(10**12,n*n)
    if type(scale)is not int or scale<100:
        raise ValueError('integer n>=64 and root scale>=100 required')
    llo,lhi=log_interval(8*n)
    eps=sqrt_interval(4*lhi/n,scale)[1]
    if eps>=1:raise ValueError('rounded shrinking bulk empty; increase root scale')
    bulk_squared=2*(1-eps)
    lower_bulk_root=sqrt_interval(bulk_squared,scale)[0]
    upper_sqrt2=sqrt_interval(2,scale)[1]
    lower_sqrt_nm1=sqrt_interval(n-1,scale)[0]
    interior=6/lower_sqrt_nm1
    radial_shift=upper_sqrt2-lower_bulk_root
    tail=F(1,8)+F(1,64*n*n)
    if tail>=F(1,4):raise ArithmeticError('uniform interpolation tail margin')
    return {'n':n,'epsilon_upper':eps,'certified_bulk_squared':bulk_squared,
            'poisson_combined_error_upper':tail,
            'whole_background_covering_strict_upper':interior+radial_shift,
            'simple_covering_upper':interior+upper_sqrt2*eps,
            'scale':'Original z, sum |z_i|^2=n',
            'scope':'Every global optimizer, all centers |u|<=sqrt(2); O(sqrt(log n/n)). No boundary O(n^-1/2) constant claimed.'}
