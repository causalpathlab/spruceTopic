import numpy as np

class ldasample():
    def __init__(self,count_matrix, n_topics,alpha,beta):
        self.count_matrix = count_matrix
        self.n_docs,self.vocab_size = count_matrix.shape
        self.n_topics  = n_topics
        self.alpha = alpha
        self.beta = beta

        self.doc_topic = np.zeros((self.n_docs,self.n_topics))
        self.topic_word = np.zeros((self.n_topics,self.vocab_size))
        self.doc_word = np.zeros((self.n_docs,self.vocab_size))
        
        self.topic_per_doc = np.zeros(self.n_docs)
        self.topic_count = np.zeros(self.n_topics)

    def initialize(self):
        for d in range(self.n_docs):
            for w in range(self.vocab_size):
                z = np.random.randint(self.n_topics)
                self.doc_word[d,w] = z
                self.doc_topic[d,z] += 1
                self.topic_word[z,w] += 1
                self.topic_per_doc[d] += 1
                self.topic_count[z]  += 1

    def conditional_distribution(self, d, w):
        left = (self.topic_word[:,w] + self.beta) / (self.topic_count + self.beta * self.vocab_size)
        right = (self.doc_topic[d,:] + self.alpha) / (self.topic_per_doc[d] + self.alpha * self.n_topics)
        p_z = left * right
        p_z /= np.sum(p_z)
        return p_z

    def run(self,maxiter=300):

        for i in range(maxiter):
            for d in range(self.n_docs):
                for w in range(self.vocab_size):
                    z = int(self.doc_word[d,w])
                    self.doc_topic[d,z] -= 1
                    self.topic_per_doc[d] -= 1
                    self.topic_word[z,w] -= 1
                    self.topic_count[z] -= 1

                    p_z = self.conditional_distribution(d, w)
                    z = np.random.multinomial(1,p_z).argmax()

                    self.doc_word[d,w] = z
                    self.doc_topic[d,z] += 1
                    self.topic_per_doc[d] += 1
                    self.topic_word[z,w] += 1
                    self.topic_count[z] += 1
        
              
    def get_doc_topic(self,title):
        import pandas as pd
        import matplotlib.pylab as plt
        import seaborn as sns
        from matplotlib.cm import ScalarMappable

        df_z = pd.DataFrame(self.doc_topic)
        df_z.columns = ["hh"+str(i)for i in df_z.columns]
        
        data_color = range(len(df_z.columns))
        data_color = [x / max(data_color) for x in data_color] 
        custom_map = plt.cm.get_cmap('coolwarm') 
        custom = custom_map(data_color)  
        df_z.plot(kind='bar', stacked=True, color=custom,figsize=(25,10))
        plt.ylabel("hidden state proportion", fontsize=18)
        plt.xlabel("samples", fontsize=22)
        plt.xticks([])
        plt.title(title,fontsize=25)
        plt.savefig("../output/gibbs_sample_"+title+".png");plt.close()



