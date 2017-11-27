import numpy as np

def gavrikov_model(delta, Ea_R, Tps, Tvn, T0):

    # constants
    a = -0.007843787493
    b = 0.1777662961
    c = 0.02371845901
    d = 1.477047968
    e = 0.1545112957
    f = 0.01547021569
    g = -1.446582357
    h = 8.730494354
    i = 4.599907939
    j = 7.443410379
    k = 0.4058325462
    m = 1.453392165
    
    X = Ea_R/Tps
    Y = Tvn/T0
    
    lnLambda_Delta = Y * (a*Y-b) + X*(c*X-d+(e-f*Y)*Y) + g*np.log(Y) + h*np.log(X) + Y*(i/X-k/(X**m)) - j
    lambda_delta = np.exp(lnLambda_Delta)
    
    return lambda_delta*delta
    
