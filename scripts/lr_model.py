import pandas as pd
import numpy as np

def log_factorial(x):
    return np.sum(np.log(np.arange(1,x+1)))

def logprob(x,lmbda):
    logfac = np.zeros((len(x)))
    for i in np.arange(len(x)):
        logfac[i] = log_factorial(x[i])
    return np.sum(x*np.log(lmbda)-lmbda-logfac)

def log_likelihood(X,lmbdas,pies):
    K = lmbdas.shape[0]
    N = X.shape[0]
    logprobs = np.zeros((K,N))
    loglik = 0
    labels = np.zeros(N)
    for n in np.arange(N):
        for k in np.arange(K):
            x = X[n,:]
            l = lmbdas[k,:]
            logprobs[k,n] = logprob(x,l) + np.log(pies[k])
        loglik += np.sum(logprobs[:,n])
        labels[n] += np.argmax(logprobs[:,n])
    return (logprobs,loglik,labels)

#test function
x = np.array( [[1, 1, 8,7], [8, 9, 1,1], [1, 1, 9,8]], dtype='float32' )
l = np.array( [[1, 1, 7,7], [8, 9, 1,1]], dtype='float32' )
p = np.array([0.33, 0.33])        
a,b,c = log_likelihood(x,l,p)
print(a)
print(b)
print(c)



