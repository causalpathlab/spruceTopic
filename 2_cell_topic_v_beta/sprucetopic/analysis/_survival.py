import pandas as pd
import numpy as np

def generate_data_it(spr,df,dfm):
    dfr = spr.interaction_topic.beta_rm.copy()
    dfl = spr.interaction_topic.beta_lm.copy()

    receptors = [ x for x in dfr.columns if x in df.columns]
    ligands = [ x for x in dfl.columns if x in df.columns]

    df = df[receptors+ligands]
    df_rl = pd.concat([dfr,dfl],axis=1)
    df_rl = df_rl[receptors+ligands]

    df_score = pd.DataFrame(np.dot(df,df_rl.T),index=df.index)
    df_score = df_score.reset_index()

    df_score = pd.merge(df_score,dfm[['icgc_donor_id','donor_vital_status','overall_time']],on='icgc_donor_id',how='left')

    df_score = df_score[['icgc_donor_id','donor_vital_status','overall_time',2,4,7,10,18,22,24]]

    for x in [2,4,7,10,18,22,24]: 
        l = list(df_score[x].sort_values())
        # n = len(l)  
        # low = np.max(l[0 : int(n/3)])
        # mid = np.max(l[int(n/3) : 2*int(n/3)])
        # df_score['topic'+str(x)] = ['low' if y<low else 'mid' if (y>low and y<mid) else 'high' for y in df_score[x]]
        lmedian = np.median(l)
        df_score['topic'+str(x)] = ['low' if y<lmedian else 'high' for y in df_score[x]]

    df_score = df_score.drop(columns=[2,4,7,10,18,22,24])

    df_score['donor_vital_status'] = ['FALSE' if x=='alive' else 'TRUE' for x in df_score['donor_vital_status']]

    return df_score

def generate_data_ct(spr,df,dfm):

    dfg = spr.cell_topic.beta_mean.copy()

    genes = [ x for x in dfg.columns if x in df.columns]

    df = df[genes]
    dfg = dfg[genes]

    df_score = pd.DataFrame(np.dot(df,dfg.T),index=df.index)
    df_score = df_score.reset_index()

    df_score = pd.merge(df_score,dfm[['icgc_donor_id','donor_vital_status','overall_time']],on='icgc_donor_id',how='left')

    for x in range(50): 
        l = list(df_score[x].sort_values())
        n = len(l)  
        low = np.max(l[0 : int(n/3)])
        mid = np.max(l[int(n/3) : 2*int(n/3)])
        df_score['topic'+str(x)] = ['low' if y<low else 'mid' if (y>low and y<mid) else 'high' for y in df_score[x]]

    df_score = df_score.drop(columns=[x for x in range(50)])

    df_score['donor_vital_status'] = ['FALSE' if x=='alive' else 'TRUE' for x in df_score['donor_vital_status']]

    return df_score