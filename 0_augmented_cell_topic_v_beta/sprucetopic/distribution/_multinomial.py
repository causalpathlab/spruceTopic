import torch

def reparameterize(mean,lnvar):
	sig = torch.exp(lnvar/2.)
	eps = torch.randn_like(sig)
	return eps.mul_(sig).add_(mean)

def multinm_log_likelihood(xx,pr, eps=1e-8):
	return torch.sum(xx * torch.log(pr+eps),dim=-1)

def multi_dir_log_likelihood(x,alpha):
	a = torch.lgamma(alpha.sum(1)) - torch.lgamma(alpha).sum(1)
	b = torch.lgamma(x + alpha).sum(1) - torch.lgamma( (x + alpha).sum(1))
	return a + b 


def kl_loss(mean,lnvar):
	return  -0.5 * torch.sum(1. + lnvar - torch.pow(mean,2) - torch.exp(lnvar), dim=-1)
