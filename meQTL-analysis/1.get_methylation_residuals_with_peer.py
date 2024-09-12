# Print requirements for this analysis
print("To run this analysis you need to have packages scipy and PEER installed")

# Set up paths and parameters
path_input_expr = input("Enter path to methylation file (in .csv format):")
path_iput_cov = input("Enter path to covariates file (in .csv format):")
no_factors = input("Enter number of PEER factors to include:")
path_out = input("Enter path to save methylation redisuals array (in .txt format):")

# Run analysis
import peer
import scipy as SP
import numpy as np

expr = SP.loadtxt(path_input_expr, delimiter=' ')
cov = SP.loadtxt(path_iput_cov, delimiter=' ')

# Set first verisuib if the model
model = peer.PEER()

model.setPhenoMean(expr)

## Set covariates
model.setCovariates(cov)

## Choose no of PEER factors you wanna include
model.setNk(no_factors)

## Build the model
model.update()

factors = model.getX()

weights = model.getW()

precision = model.getAlpha()

## Isolate the residuals array
residuals = model.getResiduals()

# Save the redisuals array
np.savetxt(path_out, residuals)