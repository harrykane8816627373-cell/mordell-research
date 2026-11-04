"""Radius 23/16 for the global-maximizer branch with the old moment ellipse."""
import subset_force as K

E=K.E


def model(settings,*,old_ellipse=False):
    if settings['scope']!='global_maximum' or old_ellipse is not True:
        raise ValueError('ellipse-Schur radius requires global maximization and the old moment ellipse')
    result=dict(settings)
    result['radius2']=min(settings['radius2'],E.rat(23,16)[1])
    result['radius']=min(settings['radius'],E.sqrt_hi(result['radius2']))
    return result


def contract(original,settings):
    # This fixed remaining-branch contractor includes the old-ellipse filter in
    # the inherited force kernel. Its feasible domain requires that ellipse;
    # it is not a radius theorem for every global maximizer.
    return K.contract(original,model(settings,old_ellipse=True),passes=1)
