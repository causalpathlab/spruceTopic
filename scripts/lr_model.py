import pandas as pd
import numpy as np


class poisson_mixture():
    def __init__(self):
        super(poisson_mixture, self).__init__()
    
    def log_factorial(self,x):
        return np.sum(np.log(np.arange(1,x+1)))
    
    def logsum(self,lp):
        m = np.max(lp)    
        return np.log(np.sum(np.exp(lp-m))) + m

    def logprob(self,x,lmbda):
        logfac = np.zeros((len(x,)))
        for i in np.arange(len(x)):
            logfac[i] = self.log_factorial(x[i])
        return np.sum(x*np.log(lmbda)-lmbda-logfac)

    def log_likelihood(self,X,lmbdas,pies):
        K = lmbdas.shape[0]
        N = X.shape[0]
        logprobs = np.zeros((K,N))
        loglik = 0
        labels = np.zeros(N)
        for n in np.arange(N):
            for k in np.arange(K):
                x = X[n,:]
                l = lmbdas[k,:]
                logprobs[k,n] = self.logprob(x,l) + np.log(pies[k])
            docloglik = self.logsum(logprobs[:,n])    
            loglik = loglik + docloglik
            logprobs[:,n] -= docloglik
            labels[n] = np.argmax(logprobs[:,n])
        
        return logprobs,loglik,labels

    def update_params(self,xprob,x,cons = 0.001):
        lmbdas = (np.dot(xprob,x)+cons)/(cons+np.sum(xprob,1)[:,np.newaxis])
        pies = np.sum(xprob,1) / x.shape[0]
        return lmbdas,pies

    def fit(self,xs,k,niter=1):

        np.random.seed(0)

        nlmbdas = xs.shape[1]
        ndata = xs.shape[0]
        
        lmbdas = np.mean(xs,0)*(1.0 + 0.5*(np.random.rand(k,nlmbdas)-0.5)) 
        pies = [1./k]*k
        
        logliks = []
        
        for s in range(niter):

            logprobs,loglik,labels = self.log_likelihood(xs,lmbdas,pies)
            qs = np.exp(logprobs)
            logliks.append(loglik)
            lmbdas,pies = self.update_params(qs,xs)
            print(s,loglik)
        
        return labels



    
    
    
# data = "wine_data.csv"
# df = pd.read_csv("../output/etm_other_data/wine_data/"+data,sep="\t")

# drop_columns_sample = [ col for col,val  in df.iloc[:,1:-1].sum(axis=0).iteritems() if val < 1 ]
# df = df.drop(drop_columns_sample,axis=1)

# X = df.iloc[:,1:-1].values
X = np.array( [[1,2,8,9], [8, 9, 1,2], [1, 1,9, 9]], dtype='float32' )

model = poisson_mixture()
labels = model.fit(X,3)
print(labels)

