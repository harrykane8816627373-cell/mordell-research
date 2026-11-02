"""Compose shared subset mass/centroid budgets with the accepted force kernel."""
import gram_force as K
import subset_budget as D

E=K.E


def contract(original,settings,passes=2):
    box=[tuple(p) for p in original]
    for step in range(passes):
        before=list(box)
        box,reason,aux=E.contract(box,settings,passes=5)
        if box is None:return box,reason,aux
        box,reason,evidence=D.propagate(box,aux[0],settings)
        if box is None:return box,reason,evidence
        box,reason,aux=K.contract(box,settings)
        if box is None:return box,reason,aux
        if box==before:break
    return box,None,aux
