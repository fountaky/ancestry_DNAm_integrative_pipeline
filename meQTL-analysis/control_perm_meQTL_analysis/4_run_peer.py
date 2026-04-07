import peer
import scipy as SP
import numpy as np
import os

# Now run analysis in a loop
n = int(os.getenv("n"))
num = n + 1

for i in range(1, num): 
    
    file = 'methylation_to_peer_{}.csv'.format(i) # Load methylation matrix to obtain residuals for
    
    expr = SP.loadtxt(file, delimiter=' ')
    cov = SP.loadtxt('tensor_covariates_eqtls.csv', delimiter=' ') # Load file with covariates of interest to control
    
    expr.shape
    cov.shape

    model = peer.PEER()

    model.setPhenoMean(expr)
    model.getPhenoMean().shape

    # Set covariates
    model.setCovariates(cov)

    # Set number of PEER factors
    model.setNk(3) 

    # Build model
    model.update()

    factors = model.getX()

    weights = model.getW()

    precision = model.getAlpha()

    ## Isolate the residuals array
    residuals = model.getResiduals()

    result = 'tensor_meth_residuals_3_{}.txt'.format(i)
    
    np.savetxt(result, residuals)
