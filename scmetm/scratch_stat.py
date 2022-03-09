
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

def bivariate_gibbs():
    #Conditionals of Multivariate Gaussians
    def p_x_given_y(y, mus, sigmas):
        mu = mus[0] + sigmas[1, 0] / sigmas[0, 0] * (y - mus[1])
        sigma = sigmas[0, 0] - sigmas[1, 0] / sigmas[1, 1] * sigmas[1, 0]
        return np.random.normal(mu, sigma)

    def p_y_given_x(x, mus, sigmas):
        mu = mus[1] + sigmas[0, 1] / sigmas[1, 1] * (x - mus[0])
        sigma = sigmas[1, 1] - sigmas[0, 1] / sigmas[0, 0] * sigmas[0, 1]
        return np.random.normal(mu, sigma)
        
    def gibbs_sampler(mus, sigmas, n_iter=5000):    
        samples = []
        y = mus[1]
        for _ in range(n_iter):
            x = p_x_given_y(y, mus, sigmas)
            y = p_y_given_x(x, mus, sigmas)
            samples.append([x, y])
        return samples
        
    # mus = np.asarray([0, 0])
    # sigmas = np.asarray([[1, .001], [.001, 1]])
    # samples = gibbs_sampler(mus, sigmas)
    # print(samples[:5])
    # burn = 100
    # x, y = zip(*samples[burn:])
    # sns.jointplot(x, y, kind='hex')
    # plt.savefig("../output/gibbsbivarn.png");plt.close()

class EM:
    def __init__(self,x,k):
        self.x = x
        self.k = k
        self.phi = np.ones(k)/k
        self.mu = np.random.uniform(x.min(), x.max(), size=self.k)
        self.std = np.ones(k)
        self.w = np.zeros((self.k, x.shape[0]))

    def e_step(self):
        for z_i in range(self.k):
            self.w[z_i] = stats.norm(self.mu[z_i],self.std[z_i]).pdf(self.x) * self.phi[z_i]
        self.w /= self.w.sum(0)

    def m_step(self):
        self.phi = self.w.mean(1)
        self.mu = (self.w * self.x).sum(1)/ self.w.sum(1)
        self.std = (( self.w * (self.x - self.mu[:,None])**2).sum(1) / self.w.sum(1))**0.5

    def fit(self,iters=300):
        for i in range(iters):
            self.e_step()
            self.m_step()
    
    def plot(self):
        fitted_m = [stats.norm(mu, std) for mu, std in zip(self.mu, self.std)]

        plt.figure(figsize=(16, 6))
        x = np.linspace(-5, 12, 150)
        plt.plot(x, fitted_m[0].pdf(x))
        plt.plot(x, fitted_m[1].pdf(x))
        plt.savefig("../output/gmm_em.png");plt.close()

    # np.random.seed(654)
    # # Draw samples from two Gaussian w.p. z_i ~ Bernoulli(phi)
    # generative_m = np.array([stats.norm(2, 1), stats.norm(5, 1.8)])
    # z_i = stats.bernoulli(0.75).rvs(100)
    # x_i = np.array([g.rvs() for g in generative_m[z_i]])
    # m =  EM(x_i,2)
    # m.fit()
    # m.plot()