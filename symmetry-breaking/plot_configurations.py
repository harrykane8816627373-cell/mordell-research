"""Optional Matplotlib illustration; the exact certificate uses integer data."""
from pathlib import Path
import argparse
import json
import math
import sys


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--module-dir',type=Path)
    args=parser.parse_args()
    if args.module_dir:
        sys.path.insert(0,str(args.module_dir))
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    here=Path(__file__).resolve().parent
    old=json.loads((here.parent/'radial-angular/population_search_verification.json').read_text(encoding='utf-8'))
    new=json.loads((here/'integer_witnesses.json').read_text(encoding='utf-8'))
    certificates=json.loads((here/'witness_verification.json').read_text(encoding='utf-8'))
    olds={r['n']:r for r in old['rows']}
    news={r['n']:r for r in new['rows']}
    certs={r['n']:r for r in certificates['rows']}
    fig,axes=plt.subplots(2,2,figsize=(10.8,10.5),layout='constrained')
    fig.suptitle('Mordell configurations: releasing the regular-ring constraint',fontsize=16,weight='bold')
    for row,n in enumerate((27,100)):
        prior=olds[n]
        p=[[0.,0.]] if prior['center'] else []
        for (m,b),sign in zip(prior['rings'],prior['signs']):
            r=math.sqrt(b/(n-1))
            angles=(2*np.pi*np.arange(m)+(np.pi if sign<0 else 0))/m
            p.extend(np.column_stack((r*np.cos(angles),r*np.sin(angles))).tolist())
            axes[row,0].add_patch(Circle((0,0),r,fill=False,lw=0.7,ls='--',color='#abb7c8'))
        p=np.array(p)
        z=np.array(news[n]['integer_points'],dtype=float)
        z-=z.mean(axis=0)
        z*=math.sqrt(n/np.sum(z*z))
        for col,points,color in [(0,p,'#62758e'),(1,z,'#126e6e')]:
            ax=axes[row,col]
            ax.scatter(points[:,0],points[:,1],s=28 if n==100 else 45,c=color,zorder=3)
            ax.set(xlim=(-1.6,1.6),ylim=(-1.6,1.6),aspect='equal')
            ax.axhline(0,lw=.5,color='#d5dce5')
            ax.axvline(0,lw=.5,color='#d5dce5')
            ax.set_xticks([-1,0,1]);ax.set_yticks([-1,0,1])
            for spine in ax.spines.values():spine.set_color('#ccd5df')
        axes[row,0].set_title(f'n = {n} | previous regular-ring candidate',loc='left',fontsize=12)
        gain=certs[n]['relative_gain_basis_points_floor']/100
        axes[row,1].set_title(f'n = {n} | certified lower bound +{gain:.2f}%',loc='left',fontsize=12)
    fig.supxlabel('Centered and normalized: sum |z_i|^2 = n. New points are rounded integer witnesses, not certified optimizers.',fontsize=10)
    fig.savefig(here/'configurations.png',dpi=170,facecolor='white')
    fig.savefig(here/'configurations.svg',facecolor='white',metadata={'Date':'2026-10-28'})
    svg=here/'configurations.svg'
    svg.write_bytes(b'\n'.join(line.rstrip() for line in svg.read_bytes().splitlines())+b'\n')


if __name__=='__main__':main()
