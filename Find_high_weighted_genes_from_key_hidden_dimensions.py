import tensorflow as tf
import numpy as np
import random
from scipy.stats import pearsonr
from scipy.spatial import distance
# import seaborn as sns
import matplotlib.pyplot as plt
import tensorflow.compat.v1 as tf

tf.disable_v2_behavior()
import heapq
import csv
import math


# deduplicate
def func1(one_list):
    return list(set(one_list))


def frequency_sort(one_list):
    a = {}
    one_list = one_list.tolist()
    for i in one_list:
        # if top_3.count(i)>1:
        a[i] = one_list.count(i)
    a = sorted(a.items(), key=lambda item: item[1], reverse=True)
    return a


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


def compare_results(A, B):
    results = list((1 - distance.correlation(A.flatten(), B.flatten())))
    results += list(Euclidean_dist(A, B))
    results += list(MAE(A, B))
    return results


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
            "GSE84133": [943, 800, 175, 486, 749]}
tf.set_random_seed(1)

# Hyper Parameters
LR = 0.0001  # learning rate
Dropout_rate = 0.5
# GSE Data
data_path = "./Data/GSE60361.npy"
X = np.load(data_path)
training_dictionary_fraction = 0.05
genes, samples = X.shape
seeds = seed_all['GSE66525']
############################# Define architectures ##################################
# tf placeholder
tf_x = tf.placeholder(tf.float32, [None, genes])  # value in the range of (0, 1)

# encoder
# Dn0 = tf.layers.dropout(tf_x, rate=Dropout_rate, training=True)
en0 = tf.layers.dense(tf_x, 1280, tf.nn.leaky_relu)
en1 = tf.layers.dense(en0, 640, tf.nn.leaky_relu)
en2 = tf.layers.dense(en1, 256, tf.nn.leaky_relu)

encoded = tf.layers.dense(en2, 10)

# decoder
de0 = tf.layers.dense(encoded, 256, tf.nn.leaky_relu)
de1 = tf.layers.dense(de0, 640, tf.nn.leaky_relu)
de2 = tf.layers.dense(de1, 1280, tf.nn.leaky_relu)

decoded = tf.layers.dense(de2, genes, tf.nn.leaky_relu)

loss = tf.losses.mean_squared_error(labels=tf_x, predictions=decoded)
train = tf.train.AdamOptimizer(LR).minimize(loss)

############################# Running ##################################
Results = {}
# seeds = random.sample(range(0, 1000), 5)
# seeds = [283, 446, 562, 114, 739]

print(seeds)
for i in range(1):

    sess = tf.Session()
    sess.run(tf.global_variables_initializer())

    X_train, X_test = random_split_train_test(X, training_dictionary_fraction, seed=seeds[i])
    # print(xi)
    # np.savetxt("GSE60361_Xi.csv", xi, delimiter=',')
    print(X.shape)  #
    print(X_train.shape)  #
    print(X_test.shape)  #
    print(X_train[0, 0:10])
    X_train = np.transpose(X_train)
    X_test = np.transpose(X_test)

    for step in range(500):
        b_x = X_train
        _, encoded_, decoded_, loss_ = sess.run([train, encoded, decoded, loss], {tf_x: b_x})

        if step % 100 == 0:
            # print('------------------Step: %d' % step + '---------------')
            # print('train loss: %.4f' % loss_)
            # plotting decoded image (second row)
            decoded_data_train = sess.run(decoded, {tf_x: b_x})
            # train_p = (1 - distance.correlation(X_train.flatten(), decoded_data_train.flatten()))
            train_pp = pearsonr(X_train.flatten(), decoded_data_train.flatten())[0]
            train_ED = Euclidean_dist(X_train, decoded_data_train)
            train_MAE = MAE(X_train, decoded_data_train)
            # print('train Pearson: %.4f' % train_p)
            # print('train Pearson_: %.4f' % train_pp)
            # print('train Euclidean_dist: %e' % train_ED)
            # print('train MAE: %.4f' % train_MAE)

            encode = sess.run(encoded, {tf_x: b_x})
            # print(encode.shape)
            # print('------------------Test---------------')
            decoded_data_testing = sess.run(decoded, {tf_x: X_test})
            encoded_data = sess.run(encoded, {tf_x: X_test})
            # test_p = (1 - distance.correlation(X_test.flatten(), decoded_data.flatten()))
            test_pp = pearsonr(X_test.flatten(), decoded_data_testing.flatten())[0]
            test_ED = Euclidean_dist(X_test, decoded_data_testing)
            test_MAE = MAE(X_test, decoded_data_testing)
            # print('test Pearson: %.4f' % test_p)
            # print('test Pearson_: %.4f' % test_pp)
            # print('test Euclidean_dist: %e' % test_ED)
            # print('test MAE: %.4f' % test_MAE)
            # print('----------------------------------------')
    #       Result = compare_results(X_test, decoded_data)
    #       print(Result)
    decoded_data_testing = sess.run(decoded, {tf_x: X_test})

    print(decoded_data_testing.shape)

    result_train = 'DeepAE4 (training)_' + str(i)
    result_test = 'DeepAE4 (testing )_' + str(i)
    Results[result_train] = [train_pp, train_ED, train_MAE]
    Results[result_test] = [test_pp, test_ED, test_MAE]
    print('----------------End Iteration: %d' % i + '------------------------')

print(data_path)
for k, v in sorted(Results.items()):
    print('\t'.join([k] + [str(x) for x in v]))


top = []
for k in range(10):
    Weights_4 = tf.get_default_graph().get_tensor_by_name('dense_3/kernel:0')
    Weights_4 = sess.run(Weights_4)
    print(Weights_4.shape)
    print(math.floor(Weights_4.shape[0] * 0.1))
    select_4 = math.floor(Weights_4.shape[0] * 0.1)
    print(select_4)
    w4 = Weights_4[:, k]
    top_4 = heapq.nlargest(select_4, range(len(w4)), w4.take)

    Weights_3 = tf.get_default_graph().get_tensor_by_name('dense_2/kernel:0')
    Weights_3 = sess.run(Weights_3)
    print(Weights_3.shape)
    print(math.floor(Weights_4.shape[0] * 0.1))
    select_3 = math.ceil(Weights_3.shape[0] * 0.1)
    print(select_3)
    top_3 = []
    for i in range(select_4):
        index = top_4[i]
        top = heapq.nlargest(select_3, range(len(Weights_3[:, index])), Weights_3[:, index].take)
        top_3 = np.hstack((top_3, top))
    print(top_3)
    print(top_3.shape)
    top_3 = frequency_sort(top_3)
    len(top_3)

    Weights_2 = tf.get_default_graph().get_tensor_by_name('dense_1/kernel:0')
    Weights_2 = sess.run(Weights_2)
    print(Weights_2.shape)
    print(math.floor(Weights_3.shape[0] * 0.1))
    select_2 = math.floor(Weights_2.shape[0] * 0.1)
    print(select_2)
    top_2 = []
    for i in range(math.floor(len(top_3) * 0.1)):
        index = int(top_3[i][0])
        top = heapq.nlargest(select_2, range(len(Weights_2[:, index])), Weights_2[:, index].take)
        top_2 = np.hstack((top_2, top))
    print(top_2)
    print(top_2.shape)
    top_2 = frequency_sort(top_2)
    len(top_2)

    Weights_1 = tf.get_default_graph().get_tensor_by_name('dense/kernel:0')
    Weights_1 = sess.run(Weights_1)
    print(Weights_1.shape)
    print(math.floor(Weights_2.shape[0] * 0.1))
    select_1 = math.floor(Weights_1.shape[0] * 0.1)
    print(select_1)
    top_1 = []
    for i in range(math.floor(len(top_2) * 0.1)):
        index = int(top_2[i][0])
        top = heapq.nlargest(select_1, range(len(Weights_1[:, index])), Weights_1[:, index].take)
        top_1 = np.hstack((top_1, top))
    print(top_1)
    print(top_1.shape)
    top_1 = frequency_sort(top_1)
    len(top_1)

    top_10 = []
    for i in range(1997):
        index = int(top_1[i][0])
        top_10.append(index)
    file = "GSE60361_top_new_" + str(k) + ".csv"
    np.savetxt(file, top_10, delimiter=',')
    top = np.hstack((top, top_10))

np.savetxt("GSE60361_top_new.csv", top, delimiter=',')
print(top.shape)
