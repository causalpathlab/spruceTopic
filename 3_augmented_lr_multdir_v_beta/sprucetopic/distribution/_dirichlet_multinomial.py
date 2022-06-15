import torch as t

def reparameterize(mean,lnvar):
	sig = t.exp(lnvar/2.)
	eps = t.randn_like(sig)
	return eps.mul_(sig).add_(mean)

def log_likelihood(x,alpha):
	a = t.lgamma(alpha.sum(1)) - t.lgamma(alpha).sum(1)
	b = t.lgamma(x + alpha).sum(1) - t.lgamma( (x + alpha).sum(1))
	return a + b 

def kl_loss(mean,lnvar):
	return  -0.5 * t.sum(1. + lnvar - t.pow(mean,2) - t.exp(lnvar), dim=-1)


