# import libraries used throughout this notebook
import numpy as np
from sklearn.metrics import precision_recall_curve, average_precision_score, log_loss, roc_auc_score, make_scorer
import keras.backend as K
from keras.layers import Input, Dense, Conv2D, Flatten, GlobalMaxPooling2D, AveragePooling2D, MaxPooling2D, Dropout, Activation
from keras.models import Model, model_from_json
from keras.activations import softmax, softplus, softsign, relu
from keras.callbacks import EarlyStopping
from keras import regularizers
import tensorflow as tf

import plotly
import plotly.graph_objs as go
import json, requests, re
from uuid import uuid4 
import pickle



def make_ohe(seq, struct):
    '''
        function returns the data in ohe shape. The columns correspond to the lexicon.
        INPUT: sequence. Sequence of amino acids or secondary structure (ss) elements.
               lexicon. Ordered list of all 20 amino acids or ss elements.
        OUTPUT: ohe_data (shape = (1, len(lexicon))
        e.g. of lexicon for ss: ["E","H","-"] --> beta, alpha, coil

        NOTE: This function can be vectorized since it will constitute a ufunc 
              and the result matrix should have a shape = (len(sequences), len(lexicon))
    '''
    # one hot encode data
    aa = ['R','H','K','D','E','S','T','N','Q','A','V','L','I','M','F' ,'Y', 'W', 'C','G','P']
    ss = ['E','H','-'] # list of secondary structure elements

    # initialize tensors
    ohe_seq = np.zeros(shape=(len(seq), 20))
    ohe_ss = np.zeros(shape=(len(struct),3))

    # encode sequence and secondary structure
    for n in range(len(seq)):
        ohe_seq[n,aa.index(seq[n])] = 1
        ohe_ss[n, ss.index(struct[n])] = 1

    # join botho tensor 
    ohe = np.vstack([ohe_seq.T, ohe_ss.T]).T #.reshape(1,len(seq),23,1)

    return ohe


def predict_full(seq):

    # predict secondary structure
    struct, random_id = predict_ss(seq)

    # extend adapters for the extremes
    seq = ''.join(['G']*15) + seq + ''.join(['G']*15)
    struct = ''.join(['-']*15) + struct + ''.join(['-']*15)

    # encode for keras and initialize results
    ohe = make_ohe(seq,struct)
    results = np.zeros(len(seq)-30)

    # roll window of predictions
    for n in range(results.shape[0]):
        results[n] = ADPred.predict(ohe[n:n+30].reshape(1,30,23,1))[0][0]
        print(results[n])
    
    # save the csv data
    csv_file = 'predictions/' + random_id + '.csv'
    fasta = 'fastas/' + random_id + '.fasta'
    with open(csv_file,'w') as f:
        f.write(','.join([str(i) for i in results]))

    # save into file the smoothed data
    y = np.array([i if i>0.8 else 0 for i in results])
    results_smooth = np.convolve(y, np.ones(20)/20, "same")    
    csv_file = 'predictions/' + random_id + '_smooth.csv'
    with open(csv_file,'w') as f:
        f.write(','.join([str(i) for i in results_smooth]))

    return results, fasta, random_id+'.csv', random_id+'_smooth.csv'


def identifier2fasta(sequence):
    page1 = 'https://www.uniprot.org/uniprot/'+ sequence.replace(' ','').replace('\n','') +'.fasta'
    page2 = 'https://www.uniprot.org/uniprot/?query='+ sequence.replace(' ','').replace('\n','') +'&sort=score'

    # case is a uniprot systematic name 
    try:
        page = requests.get(page1).text 
    except Exception as e:
        return -1
        #print('fasta page could not be downloaded in the first exception',str(e))
    
    # case is a common name (e.g. gcn4)
    if page[0] == ">":
        return clean_input(page)
    else:
        try:
            page = requests.get(page2).text
            identifier = re.search("<tr id=\".{1,10}\"", page).group()[7:].replace('"','')
            return clean_input(requests.get('https://www.uniprot.org/uniprot/'+ identifier +'.fasta').text)

        except Exception as e:
            return -1
            #print('protein name could not be extracted from uniprot site',str(e))
    
    return -1


def clean_input(fasta):
    if fasta[0]==">":
        fasta = ''.join(fasta.split('\n')[1:])
    return fasta.replace('\n','').replace('\r','').upper()
