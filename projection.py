# projection.py

import numpy as np

def projectBinary(x):
    for row in range(len(x)):
        if x[row] <= 0.5:
            x[row] = 0
        else:
            x[row] = 1
    return x

def projectHyperplane(x,A,b):
    for row in range(A.shape[0]):
        if np.dot(A[row,:],x) - b[row] > 0:
            x = x - ((( (np.dot(np.array(A[row,:]),x))[0] - b[row] )/ ((np.array(np.square(A).sum(axis=1))[row])[0]) ) * np.array(A[row,:].flatten())[0])
    return x
