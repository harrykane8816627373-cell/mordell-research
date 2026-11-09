"""Interval construction of the candidate-saturating polynomial and its local covers."""
from pathlib import Path
import importlib.util
import json
from fractions import Fraction as F
from math import comb, factorial
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
p=ROOT/'twelve-point-sharp-contacts/verify_all.py'
sp=importlib.util.spec_from_file_location('sharp_contacts',p)
sharp=importlib.util.module_from_spec(sp);sp.loader.exec_module(sharp)
ZERO,ONE,SCALE,Ball=sharp.ZERO,sharp.ONE,sharp.SCALE,sharp.Ball
exact=sharp.exact


def polyproduct(a,b):
    result=[ZERO]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            result[i+j]+=x*y
    return result


def radial_features(k,nodes,means,ratios):
    raw=[[exact(comb(j,i)*(-1)**(j-i)) for i in range(j+1)] for j in range(5)]
    if k==0:
        q=[ONE]
        for node in nodes:
            q=polyproduct(q,[-node,ONE])
        return [[-means[j-1]]+[ZERO]*(j-1)+[ONE] for j in range(1,5)]+[[10*x for x in polyproduct(q,q)]]
    if k in (3,6,9):
        for i in range(1,5):
            raw[i][0]-=ratios[k][i]
        return raw[1:]
    return raw


def pair_polynomial(B,nodes,means,ratios,weights):
    H=[{(0,0,0):1},{(0,0,1):1}]
    for k in range(2,12):
        current={}
        for (i,j,l),a in H[-1].items(): current[i,j,l+1]=current.get((i,j,l+1),0)+2*a
        for (i,j,l),a in H[-2].items(): current[i+1,j+1,l]=current.get((i+1,j+1,l),0)-a
        H.append(current)
    coefficients={(0,0,0):B/66-1,(1,0,0):exact(F(6,11)),(0,1,0):exact(F(6,11))}
    def add(key,value):
        coefficients[key]=coefficients.get(key,ZERO)+value
    def block(k,features,Q):
        for row,f in enumerate(features):
            for col,g in enumerate(features):
                a=Q[row][col]
                if a.lo==a.hi==0:continue
                for i,fi in enumerate(f):
                    for j,gj in enumerate(g):
                        value=a*fi*gj
                        for (p,q,l),b in H[k].items():add((i+p,j+q,l),-b*value)
                        add((i+j+k,0,0),-value/22)
                        add((0,i+j+k,0),-value/22)
    for k in range(12):
        features=radial_features(k,nodes,means,ratios)
        dim=len(features)
        Q=[[weights[f'Q{k}_{min(i,j)}{max(i,j)}'] for j in range(dim)] for i in range(dim)]
        block(k,features,Q)
    features=radial_features(1,nodes,means,ratios)
    A=[[weights.get(f'A1_{min(i,j)}{max(i,j)}',ZERO) for j in range(5)] for i in range(5)]
    block(1,features,A)
    for k,c in ((0,F(1,12)),(2,F(-3,4)),(3,F(-2,3))):
        a=[weights.get(f'A{k}_a{i}',ZERO) for i in range(5)]
        if k!=2:
            a[0]=-sharp.dot(a[1:],ratios[k][1:])
        v=[exact(c),ONE,ZERO,ZERO,ZERO]
        A=[[v[i]*a[j]+a[i]*v[j] for j in range(5)] for i in range(5)]
        block(k,features,A)
    return coefficients


def reconstruct():
    names,rows,matrix,contacts,ratios,nodes,means=sharp.build()
    minor=json.loads((p.parent/'minor.json').read_bytes())
    seed=json.loads((ROOT/'twelve-point-banded-upper/kernel.json').read_bytes())
    parameters=sharp.seed_parameters(seed,names)
    target,B=sharp.right_hand_side(contacts)
    residual=[b-sharp.dot(row,parameters) for row,b in zip(matrix,target)]
    chosen_rows, chosen_columns = minor['selected_rows'], minor['selected_columns']
    sharp.need(len(set(chosen_rows)) == len(set(chosen_columns)) == 44, '44 distinct rows and columns')
    sharp.need([rows[i] for i in chosen_rows] == minor['selected_row_names'] and
               [names[i] for i in chosen_columns] == minor['selected_column_names'],
               'fixed contact minor and coordinate names')
    pivots,determinant,correction=sharp.certify_minor(matrix,chosen_rows,chosen_columns,residual)
    sharp.need(sharp.abs_lower(determinant)>0,'sharp contact determinant excludes zero')
    for index,change in zip(minor['selected_columns'],correction):
        parameters[index]+=change
    for row,b in zip(matrix,target):
        enclosure=sharp.dot(row,parameters)-b
        sharp.need(enclosure.lo<=0<=enclosure.hi,'all contact equations have zero in their enclosure')
    weights=dict(zip(names,parameters))
    grams=sharp.check_grams(weights)
    curvatures,normal=sharp.check_contacts(contacts,nodes,means,ratios,names,parameters)
    C=pair_polynomial(B,nodes,means,ratios,weights)
    for c in contacts:
        args=[c[x] for x in 'uvt']
        d=args[0]+args[1]-2*args[2]
        v=sum((a*sharp.ball_power(args[0],i)*sharp.ball_power(args[1],j)*sharp.ball_power(args[2],k)
                for (i,j,k),a in C.items()),ZERO)
        log=sharp.log_ball(d)
        sharp.need(v.lo<=log.hi and log.lo<=v.hi,'expanded polynomial and exact contact overlap')
    for u,v,t in ((F(1,5),F(1),F(-1,4)),(F(6,5),F(13,10),F(1,3)),
                  (F(9,50),F(4,25),F(-1,9)),(F(7,5),F(11,10),F(-3,5))):
        point=dict(u=exact(u),v=exact(v),t=exact(t))
        column_names,jets=sharp.contact_columns(point,nodes,means,ratios)
        sharp.need(column_names==names,'independent evaluator uses the same ordered parameters')
        direct=B/66-1+exact(F(6,11)*(u+v))+sharp.dot(parameters,[jet.v for jet in jets])
        powers=[{j:sharp.ball_power(exact(a),j) for j in range(20)} for a in (u,v,t)]
        expanded=sum((a*powers[0][i]*powers[1][j]*powers[2][k] for (i,j,k),a in C.items()),ZERO)
        sharp.need(expanded.lo<=direct.hi and direct.lo<=expanded.hi,'independent expanded coefficient evaluation')
    proof={'contact_minor_order':len(pivots),'exact_positive_definite_Gram_blocks':len(grams),
           'contacts_with_strict_point_curvature':len(curvatures),
           'min_gram_pivot_strict_lower_power_of_two':sharp.lower_power(min(grams))-1,
           'inward_contact_derivative_strict_lower_power_of_two':sharp.lower_power(normal)-1,
           'independent_polynomial_interval_crosschecks':18}
    return C,contacts,proof


def translate(poly,center):
    result=dict(poly)
    for axis,c in enumerate(center):
        powers=[ONE]
        for i in range(1,20):powers.append(powers[-1]*c)
        output={}
        for key,value in result.items():
            m=key[axis]
            for j in range(m+1):
                index=key[:axis]+(j,)+key[axis+1:]
                output[index]=output.get(index,ZERO)+value*comb(m,j)*powers[m-j]
        result=output
    return result


def diff(poly,axis):
    return {p[:axis]+(p[axis]-1,)+p[axis+1:]:a*p[axis] for p,a in poly.items() if p[axis]}


def bound(poly,windows):
    powers=[]
    for value in windows:
        line=[ONE]
        for i in range(1,20):line.append(line[-1]*value)
        powers.append(line)
    return sum((a*powers[0][i]*powers[1][j]*powers[2][k] for (i,j,k),a in poly.items()),ZERO)


def symmetric(r):
    z=exact(r)
    return Ball(-z.hi,z.hi)


def add(poly,key,value):
    poly[key]=poly.get(key,ZERO)+value
def polynomial_cubic_gap(C,contact):
    center=[contact[x] for x in 'uvt']
    d=center[0]+center[1]-2*center[2]
    poly=translate(C,center)
    # The constant is exactly zero after subtracting the contact logarithm.
    poly[(0,0,0)]=ZERO
    e=[1,1,-2]
    for degree,sign in ((1,-1),(2,1),(3,-1)):
        scale=sign/(degree*sharp.ball_power(d,degree))
        for i in range(degree+1):
            for j in range(degree-i+1):
                k=degree-i-j
                coef=factorial(degree)//(factorial(i)*factorial(j)*factorial(k))*e[2]**k
                add(poly,(i,j,k),coef*scale)
    return poly,center


def certify_local_boxes(coefficients, contacts):
    records=json.loads((HERE/'local_boxes.json').read_bytes())
    sharp.need(len(records)==len(contacts)==14,'one local certificate for every contact')
    boxes=[]
    pivots=[]
    for index,(record,contact) in enumerate(zip(records,contacts)):
        sharp.need(record['index']==index,'ordered contact certificate')
        radius=F(record['radius'])
        sharp.need(radius>0 and radius.denominator & (radius.denominator-1)==0,
                   'positive dyadic local radius')
        polynomial,center=polynomial_cubic_gap(coefficients,contact)
        gradient=[diff(polynomial,i) for i in range(3)]
        hessian=[[diff(gradient[i],j) for j in range(3)] for i in range(3)]
        widths=[radius,radius,2*radius if contact['boundary'] else radius]
        windows=list(map(symmetric,widths))
        u,v,t=[a+b for a,b in zip(center,windows)]
        sharp.need(min(u.lo,v.lo)>0 and (u+v-2*t).lo>0,'local continuation remains in positive distance')
        if contact['boundary']:
            sharp.need(record['preconditioner'] is None,'boundary uses the tangent coordinates directly')
            gt=bound(gradient[2],windows)
            sharp.need(gt.lo>0,'cubic gap increases inward throughout the boundary box')
            face=-sharp.sym.sqrt_ball(u*v,Ball,SCALE)
            t_margin=2*exact(radius).lo
            sharp.need(center[2].hi-t_margin<=face.lo<=face.hi<=center[2].lo+t_margin,
                       'Gram projection belongs to the larger box for the exact contact')
            matrix=[[bound(hessian[i][j],windows) for j in range(3)] for i in range(3)]
            tangent=[[ONE,ZERO,face/(2*u)],[ZERO,ONE,face/(2*v)]]
            mixed=face/(4*u*v)
            second=[[-face/(4*u.square()),mixed],[mixed,-face/(4*v.square())]]
            lower=[[sum((tangent[i][a]*matrix[a][b]*tangent[j][b]
                         for a in range(3) for b in range(3)),ZERO)+gt*second[i][j]
                    for j in range(2)] for i in range(2)]
        else:
            W=[list(map(sharp.dyadic,row)) for row in record['preconditioner']]
            sharp.need(len(W)==3 and all(len(row)==3 for row in W),'3 by 3 congruence witness')
            J=[[{} for j in range(3)] for i in range(3)]
            for i in range(3):
                for j in range(3):
                    for a in range(3):
                        for b in range(3):
                            for key,value in hessian[a][b].items():
                                add(J[i][j],key,value*W[a][i]*W[b][j])
            lower=[[bound(J[i][j],windows) for j in range(3)] for i in range(3)]
        pivots.append(sharp.ldl_lower(lower,f'cubic gap on full local box {index}'))
        # Intersection of cubes for all possible centers. Each such box is
        # contained in the cube about the exact root certified above.
        box=tuple((F(x.hi,SCALE)-radius,F(x.lo,SCALE)+radius) for x in center)
        sharp.need(all(lo<hi for lo,hi in box),'nonempty rational local cover')
        boxes.append(box)
    return boxes,{'local_boxes':len(boxes),'minimum_nominal_local_radius':'1/2048',
                  'local_cubic_hessian_min_pivot_strict_lower_power_of_two':sharp.lower_power(min(pivots))-1}
