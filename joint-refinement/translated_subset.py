"""Fixed-centre subset variance cuts after the accepted ellipse-Schur kernel."""
import ellipse_schur_radius as K
import subset_budget as D

E=K.E


def bound(box,qbox,members,outside,coefficients):
    """An outward upper bound on the subset distance sum, in squared mesh units."""
    m=len(members)
    center=tuple(sum(qbox[i][axis][0]+qbox[i][axis][1] for i in members)//(2*m)
        for axis in range(2))
    fixed=tuple(E.point(v) for v in center)
    shifted_mass=sum(E.cnorm(E.cs(qbox[i],fixed))[1] for i in members)
    total=D.centroid(box,qbox,members,outside,coefficients)
    remainder=E.cs(total,tuple(E.point(m*v) for v in center))
    upper=m*shifted_mass-E.cnorm(remainder)[0]
    return {'members':list(members),'center':list(center),'upper':upper,
        'strict_threshold_lower':m*m*(m-1)*E.S//8}


def exclusion(box,qbox,settings):
    if settings['scope']!='global_maximum':
        raise ValueError('strict translated-subset bounds require global maximization')
    for members,outside,coefficients in D.RECORDS:
        try:evidence=bound(box,qbox,members,outside,coefficients)
        except ValueError:
            return {'reason':'translated_subset_geometry','members':list(members)}
        if evidence['upper']<=evidence['strict_threshold_lower']:
            return dict(evidence,reason='translated_subset_budget')
    return None


def contract(original,settings):
    if settings['scope']!='global_maximum':
        raise ValueError('translated-subset composition requires global maximization')
    box,reason,aux=K.contract(original,settings)
    if box is None:return box,reason,aux
    # The inherited final interval-engine result is (point_boxes, distances).
    evidence=exclusion(box,aux[0],settings)
    if evidence is not None:return None,evidence['reason'],evidence
    return box,None,aux
