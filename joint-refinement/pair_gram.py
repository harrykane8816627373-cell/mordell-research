"""All complex two-column principal minors of 16(I-11*/8)-B*B."""
import force_enclosure as F

E=F.E


def exclusion(qbox,inverse=None):
    if inverse is None:
        try:inverse,distances=E.invert_edges(qbox)
        except ZeroDivisionError:return None
    else:distances=[E.cnorm(E.cs(qbox[i],qbox[j])) for i,j in E.PAIRS]
    b=[[E.ci(0) for _ in range(8)] for _ in range(8)]
    off_norm_lower=[0]*8
    for (i,j),zinv,distance in zip(E.PAIRS,inverse,distances):
        w=E.cm(zinv,zinv)
        b[i][i]=E.ca(b[i][i],w);b[j][j]=E.ca(b[j][j],w)
        b[i][j]=b[j][i]=E.cn(w)
        if distance[1]<=0:return None
        lower=E.S**3//(distance[1]*distance[1])
        off_norm_lower[i]+=lower;off_norm_lower[j]+=lower
    upper=[14*E.S-E.cnorm(b[i][i])[0]-off_norm_lower[i] for i in range(8)]
    for i,value in enumerate(upper):
        if value<0:return {'reason':'Hessian_Gram_diagonal','point':i,'diagonal_upper':value}
    for i,j in E.PAIRS:
        gram=E.ci(0)
        for k in range(8):gram=E.ca(gram,E.cm(E.cc(b[k][i]),b[k][j]))
        entry=E.cs(E.ci(-2),gram)
        lower=E.cnorm(entry)[0]
        if lower*E.S>upper[i]*upper[j]:
            return {'reason':'Hessian_Gram_pair','pair':[i,j],
                'off_diagonal_norm_lower':lower,'diagonal_uppers':[upper[i],upper[j]]}
    return None
