import numpy as np
import search as s

a,b,c,r,x=.11183669980128009,1.0074790085402798,1.188144595492442,1.096269848082956,-.07294054523381732
A,B,C=np.sqrt([a,b,c])
p,q=np.sqrt([(r+x)/2,(r-x)/2])
z=np.array([A,-A,1j*B,-1j*B,C,-C,p+1j*q,-p-1j*q,p-1j*q,-p+1j*q])
v=np.array([-630j,630j,63,-63,224j,-224j,-118-105j,118+105j,118-105j,-118+105j])/1000
d=z[:,None]-z[None,:]+np.eye(10)
k=v[:,None]-v[None,:]
J=s.jacobian(z/np.sqrt(.9))/.9
gp=(4*np.sum((k*k/d**3),axis=1).conjugate()).view(float)
w=(-np.linalg.pinv(J,rcond=1e-9)@gp).view(complex)
print('w',w, np.linalg.norm(J@w.view(float)+gp))
roundw=np.round(w*1e3)/1e3
print('rounded w', roundw,np.linalg.norm(J@roundw.view(float)+gp))
l=roundw[:,None]-roundw[None,:]
print('Q0',-np.sum((k/d)**2).real-10*np.vdot(v,v).real)
print('bounds: radius, row sum remainder, spectral remainder')
for t in (.01,.012,.015,.018,.02,.025,.03,.035,.04,.05):
    W=(abs(l*d-6*k*k)+abs(l)*np.sqrt(2)*t)/(abs(d)-np.sqrt(2)*t)**4
    U=np.diag(W.sum(axis=1))-W
    print(t,2*max(W.sum(axis=1))*t*t,np.linalg.eigvalsh(U)[-1]*t*t)
