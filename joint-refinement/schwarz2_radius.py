"""Second-order Schwarz radius bound for an eight-point global maximizer."""
import subset_force as K

E=K.E


def model(settings):
    if settings['scope']!='global_maximum':
        raise ValueError('second-order Schwarz radius is for the global-maximizer reduction')
    result=dict(settings)
    result['radius2']=min(settings['radius2'],E.rat(231,160)[1])
    result['radius']=min(settings['radius'],E.sqrt_hi(result['radius2']))
    return result


def contract(original,settings):
    return K.contract(original,model(settings),passes=1)
