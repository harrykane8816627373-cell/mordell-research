"""Optional NumPy diagnostics of pinned certified centers; numerical only."""
from pathlib import Path
import json
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent


def normalized(points):
    z=np.asarray(points,dtype=float).reshape(-1,2).copy()
    z-=z.mean(axis=0)
    return z*np.sqrt(len(z)/np.sum(z*z))


def force_jacobian(points):
    z=normalized(points);n=len(z);lam=(n-1)/2
    jac=-lam*np.eye(2*n)
    for i in range(n):
        for j in range(i):
            delta=z[i]-z[j];d2=np.dot(delta,delta)
            block=np.eye(2)/d2-2*np.outer(delta,delta)/d2**2
            for a,b,sign in ((i,i,1),(j,j,1),(i,j,-1),(j,i,-1)):
                jac[2*a:2*a+2,2*b:2*b+2]+=sign*block
    return jac


def cluster_labels(values,tolerance=1e-8):
    centers=[];labels=[]
    for value in values:
        value=np.array(value,dtype=float,ndmin=1)
        matches=[j for j,c in enumerate(centers) if np.max(np.abs(value-c))<tolerance]
        if matches:labels.append(matches[0])
        else:centers.append(value);labels.append(len(centers)-1)
    return centers,labels


def diagnostics(points):
    z=normalized(points);n=len(z);q=z[:,0]+1j*z[:,1];r2=np.sum(z*z,axis=1)
    radial,labels=cluster_labels(r2)
    contacts=[]
    for i in range(n):
        for j in range(i):
            a,b=sorted((r2[i],r2[j]));t=float(np.dot(z[i],z[j]));contacts.append((a,b,t))
    contacttypes,clabels=cluster_labels(contacts)
    distances,dlab=cluster_labels([sum((z[i]-z[j])**2) for i in range(n) for j in range(i)])
    rotations=[]
    for k in range(2,13):
        transformed=q*np.exp(2j*np.pi/k)
        d2=np.abs(transformed[:,None]-q[None,:])
        indices=np.argmin(d2,axis=1)
        residual=float(np.max(np.min(d2,axis=1)))
        if len(set(indices.tolist()))==n and residual<1e-8:
            rotations.append({'order':k,'max_matching_distance':residual,'permutation':indices.tolist()})
    radius=max(abs(q));moments=[]
    for k in range(1,19):
        value=np.sum((q/radius)**k)/n
        moments.append({'k':k,'normalized_abs':float(abs(value)),
                        'radial_weight1_abs':float(abs(np.sum((q/radius)**k*r2)/n)),
                        'radial_weight2_abs':float(abs(np.sum((q/radius)**k*r2**2)/n))})
    radialvector=z.ravel();rotation=np.column_stack((-z[:,1],z[:,0])).ravel()
    _,_,vh=np.linalg.svd(np.vstack((radialvector,rotation)),full_matrices=True)
    tangent=vh[2:].T
    stiffness=-tangent.T@force_jacobian(z)@tangent
    eigenvalues=np.linalg.eigvalsh(stiffness)
    shell=[{'squared_radius':float(c[0]),'population':labels.count(i)} for i,c in enumerate(radial)]
    shell.sort(key=lambda r:r['squared_radius'])
    return {'n':n,'diagnostic_scope':'floating geometry and eigenvalues only; no new exact symmetry or Hessian certificate',
            'radial_shells_at_tolerance_1e_8':shell,
            'radial_shell_count':len(radial),'pair_count':n*(n-1)//2,
            'SO2_pair_contact_types_at_tolerance_1e_8':len(contacttypes),
            'distance_squared_types_at_tolerance_1e_8':len(distances),
            'rotation_matches':rotations,'normalized_complex_moments':moments,
            'sphere_tangent_minus_force_Jacobian_eigenvalues':eigenvalues.tolist(),
            'sphere_tangent_min_curvature':float(eigenvalues[0]),
            'sphere_tangent_max_curvature':float(eigenvalues[-1]),
            'sphere_tangent_condition_number':float(eigenvalues[-1]/eigenvalues[0]) if eigenvalues[0]>0 else None}


def main():
    data=json.loads((ROOT/'validated-maxima/certificates.json').read_text())
    rows=[]
    for row in data['rows']:
        exponent=data['center_exponent']
        v=[float(q/2**exponent) for q in row['center_numerators']]
        points=np.array(v[:1]+[0.0]+v[1:]).reshape(row['n'],2)
        result=diagnostics(points)
        result['existing_exact_classification']=row['classification']
        rows.append(result)
    out={'date':'2026-11-08','numpy_version':np.__version__,'rows':sorted(rows,key=lambda r:r['n'])}
    (HERE/'geometry_report.json').write_text(json.dumps(out,indent=2)+'\n')
    for r in out['rows']:
        print(r['n'],'rotation orders',[a['order'] for a in r['rotation_matches']],
              'radial types',r['radial_shell_count'],'pairtypes',r['SO2_pair_contact_types_at_tolerance_1e_8'],
              'mincurvature',r['sphere_tangent_min_curvature'],flush=True)


if __name__=='__main__':main()
