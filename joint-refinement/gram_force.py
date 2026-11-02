"""Retain phase-aware force contraction and add full-stability pair Gram cuts."""
import phase_force as P
import pair_gram as G

E=P.E


def contract(original,settings):
    box,reason,aux=P.contract(original,settings)
    if box is None:return box,reason,aux
    evidence=G.exclusion(aux[0],aux[1])
    if evidence is not None:return None,evidence['reason'],evidence
    return box,None,aux
