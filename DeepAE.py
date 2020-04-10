#import tensorflow as tf
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()
import numpy as np
import random
from scipy.stats import pearsonr
from analyze_predictions import *
from scipy.spatial import distance
#import seaborn as sns
import matplotlib.pyplot as plt
import time
import pandas as pd
import sys, getopt
from sklearn.metrics import mean_absolute_error
from scipy import stats
from scipy.spatial.distance import cdist
from scipy.spatial import distance
from scipy.stats import spearmanr, pearsonr



def Euclidean_dist_new(A, B):
    dis = cdist(A, B, metric='euclidean')
    return np.mean(dis)
def Euclidean_dist(A, B):
    C = A - B
    return sum(map(sum, C * C)) ** 0.5


def MAE(A, B):  ## Mean Absolute Error
    C = A - B
    return sum(map(sum, C * C)) / (C.shape[0] * C.shape[1])


def random_split_train_test(X0, training_dictionary_fraction, seed, dictionary_size=0.5, biased_training=0.):
    training_dictionary_size = max(int(training_dictionary_fraction * X0.shape[1]), 5)
    if dictionary_size < 1:
        dictionary_size = dictionary_size * training_dictionary_size
    dictionary_size = int(dictionary_size)
    xi = np.zeros(X0.shape[1], dtype=np.bool)
    if biased_training > 0:
        np.random.seed(seed)
        i = np.random.randint(len(xi))
        dist = distance.cdist([X0[:, i]], X0.T, 'correlation')[0]
        didx = np.argsort(dist)[1:int(biased_training * training_dictionary_size) + 1]
    else:
        didx = []
    xi[didx] = True
    if biased_training < 1:
        remaining_idx = np.setdiff1d(range(len(xi)), didx)
        np.random.seed(seed)
        xi[np.random.choice(remaining_idx, training_dictionary_size - xi.sum(), replace=False)] = True
    xa = X0[:, xi]
    xb = X0[:, np.invert(xi)]
    return xa, xb


if __name__ == "__main__":
    seed_all = {"GSE71858": [272, 781, 692, 219, 292],  #
                "GSE60361": [283, 446, 562, 114, 739],  #
                "GSE62270": [629, 685, 953, 595, 378],  #
                "GSE48968": [623, 19, 621, 802, 557],  #
                "GSE52529": [550, 939, 76, 260, 328],  #
                "GSE77564": [475, 649, 316, 639, 741],
                "GSE78779": [152, 866, 808, 796, 184],  #
                "GSE10247": [702, 217, 944, 338, 701],  #
                "GSE69405": [317, 470, 798, 283, 695],
                "GSE45235": [282, 713, 521, 717, 517],  #
                "GSE25038": [480, 402, 413, 64, 574],
                "mass_cytomatry": [943, 800, 175, 486, 749],
                "GSE66525": [943, 800, 175, 486, 749],
                "GSE84133": [943, 800, 175, 486, 749]
                }
    filename = {"data1": "GSE60361.npy",  #
                #"data2": "GSE62270.npy",  #
                #"data3": "GSE48968.npy",  #
                #"data4": "GSE52529.npy",  #
                #"data5": "GSE78779.npy",  #
                #"data6": "GSE10247.npy",  #
                #"data7": "GSE69405.npy",
                #"data8": "GSE25038.txt",
                #"data9": "GSE45235.txt",
                #"data10": "mass_cytomatry.txt",
                "data11": "GSE66525.txt",
                "data12": "GSE84133.txt"
                }


    for _, f_name in filename.items():
        #lead GSE data
        data_path = "./Data/" + f_name
        if f_name[-3:] == "txt":
            X = np.loadtxt(data_path)
        else:
            X = np.load(data_path)
        f_name = f_name[:-4]
        print(f_name)

        start = time.time()
        tf.set_random_seed(1)
        # Hyper Parameters
        LR = 0.0001  # learning rate
        Dropout_rate = 0.5
        # GSE Data

        if(f_name == "GSE25038" or "GSE45235"):
            X = np.delete(X, -1, axis=1)  ##### just for metabolic profiling data
        elif(f_name == "mass_cytomatry"):
            X = np.delete(X, (0, 1, 2), axis=1)  ## for mass_cytomatry

        training_dictionary_fraction = 0.05
        #training_dictionary_fraction = 0.20 #just for small dataset
        genes, samples = X.shape


        seeds = seed_all[f_name]
        print(seeds)
        measurement_all = (1, 2, 5, 10, 25, 50, 100)

        for i in range(len(measurement_all)):
            measurements = measurement_all[i]
            ############################# Define architectures ##################################
            # tf placeholder
            tf_x = tf.placeholder(tf.float32, [None, genes])  # value in the range of (0, 1)

            # encoder
            # Dn0 = tf.layers.dropout(tf_x, rate=Dropout_rate, training=True)
            en0 = tf.layers.dense(tf_x, 1280, tf.nn.leaky_relu)
            en1 = tf.layers.dense(en0, 640, tf.nn.leaky_relu)
            en2 = tf.layers.dense(en1, 256, tf.nn.leaky_relu)

            encoded = tf.layers.dense(en2, measurements)

            # decoder
            de0 = tf.layers.dense(encoded, 256, tf.nn.leaky_relu)
            de1 = tf.layers.dense(de0, 640, tf.nn.leaky_relu)
            de2 = tf.layers.dense(de1, 1280, tf.nn.leaky_relu)

            decoded = tf.layers.dense(de2, genes, tf.nn.leaky_relu)

            loss = tf.losses.mean_squared_error(labels=tf_x, predictions=decoded)
            train = tf.train.AdamOptimizer(LR).minimize(loss)
            ############################# Running ##################################
            sys.stdout = open('./Data/DeepAE/' + f_name +'_zscore_new.log', 'a')
            print('measurement:' + str(measurements))
            itr = 0
            while itr < 5:
                Results = {}
                sess = tf.Session()
                sess.run(tf.global_variables_initializer())
                #if tf.test.gpu_device_name():
                #    print('Default GPU Device: {}'.format(tf.test.gpu_device_name()))
                #else:
                #    print("Please install GPU version of TF")

                X_train, X_test = random_split_train_test(X, training_dictionary_fraction, seed=seeds[itr])

                print(X.shape)  #
                print(X_train.shape)  #
                print(X_test.shape)  #
                #print(X_train[0, 0:10])

                X_train = np.transpose(X_train)
                X_test = np.transpose(X_test)
                X_train = stats.zscore(X_train, axis=1)
                X_test = stats.zscore(X_test, axis=1)

                for step in range(400):
                    b_x = X_train
                    _, encoded_, decoded_, loss_ = sess.run([train, encoded, decoded, loss], {tf_x: b_x})

                    if step % 80 == 0:
                        # print('------------------Step: %d' % step + '---------------')
                        decoded_data_train = sess.run(decoded, {tf_x: b_x})
                        train_pp = pearsonr(X_train.flatten(), decoded_data_train.flatten())[0]
                        train_ED = Euclidean_dist(X_train, decoded_data_train)
                        train_MAE = MAE(X_train, decoded_data_train)
                        train_spearman = spearmanr(X_train.flatten(), decoded_data_train.flatten())[0]
                        encode = sess.run(encoded, {tf_x: b_x})


                        decoded_data_testing = sess.run(decoded, {tf_x: X_test})
                        test_pp = pearsonr(X_test.flatten(), decoded_data_testing.flatten())[0]
                        test_ED = Euclidean_dist_new(X_test, decoded_data_testing)
                        test_MAE = mean_absolute_error(X_test, decoded_data_testing)
                        test_spearman = spearmanr(X_test.flatten(), decoded_data_testing.flatten())[0]


                print(decoded_data_testing.shape)

                result_train = 'DeepAE4 (training)_' + str(itr+1)
                result_test = 'DeepAE4 (testing )_' + str(itr+1)
                Results[result_train] = [train_pp, train_ED, train_MAE, train_spearman]
                Results[result_test] = [test_pp, test_ED, test_MAE, test_spearman]


                print(data_path)
                for k, v in sorted(Results.items()):
                    print('\t'.join([k] + [str(x) for x in v]))

                end = time.time()
                print(end-start)
                itr += 1
                print('----------------End Iteration: %d' % itr + '------------------------')
