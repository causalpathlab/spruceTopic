import pandas as pd
from sampler import gibbs_lda
import importlib

data = "nydata.csv"
df = pd.read_csv("../output/news_data/"+data,sep="\t")
df = df.sample(n=200)

count_matrix = df.iloc[:,1:-1].values

n_topics = 16
alpha = 0.001
beta = 0.001
m = gibbs_lda.ldasample(count_matrix, n_topics,alpha,beta)
m.initialize()
m.run()
m.get_doc_topic("lda")