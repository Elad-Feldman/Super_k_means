import numpy as np
import argparse
import sys
import time
import math
import spkmeans
#python3 kmeans_pp.py 3 100 input_1_db_1.txt input_1_db_2.txt

## NOTE - WORK with python 3 only


def load_data_to_dots(filename): # load data into 2d matrix
    dots_list = []
    file = open(filename, 'r')
    Lines = file.readlines()

    for line in Lines:
        dot = [float(word) for word in line.split(sep=",")]
        dots_list.append(dot)
    return dots_list




def find_index_nearest_cluster(w):
    # returns a random index by the chances provided by w
    return int(np.random.choice(range(len(w)), 1, replace=False, p=w))



def get_cluster_distance(dot, cluster):
    #calcute distance between dot and cluster
    a = np.array(cluster)
    b = np.array(dot)
    dis = np.linalg.norm(a-b)
    dis = math.pow(dis, 2)
    return dis


def find_initial_clusters(n_dots, k):
    # this is kmeans++ find initail clusters
    n = len(n_dots)
    random_index = np.random.choice(n)
    first_cluster = n_dots[random_index]
    clusters_list = [first_cluster]
    indices = [random_index]

    distances_list = [float("inf")] * n
    for z in range(1, k):
        for j, dot in enumerate(n_dots):
            d = get_cluster_distance(dot, clusters_list[z - 1])
            distances_list[j] = min(distances_list[j], d)

        sum_d = sum(distances_list)
        weights = [d/sum_d for d in distances_list ]
        index = find_index_nearest_cluster(weights)
        indices.append(index)
        clusters_list.append(n_dots[index])
    assert (k == len(indices))
    return indices

def get_cluster_list (observations, indices):
    # get first culsters from index list
    clusters = []
    for ind in indices:
        clusters.append(observations[ind])
    return clusters


def check_param(k, goal,filename):
    # chack input
    assert type(k) is int, "Invalid Input!"
    assert type(goal) is str, "Invalid Input!"
    assert type(filename) is str, "Invalid Input!"

def parse2():
    # pasre arguments
    if len(sys.argv)== 4:
        k = int( sys.argv[1])
        goal = str(sys.argv[2])
        filename = str(sys.argv[3])
    else:
        assert False,"Invalid Input!"

    check_param(k, goal, filename)
    return k, goal, filename



if __name__ == "__main__":
    T0 = time.process_time()
    np.random.seed(0)
    k, goal, filename = parse2()
    observations = load_data_to_dots(filename)
    assert len(observations) > k, "An Error Has Occured"
    T,k = spkmeans.get_flag(goal, k, observations)



    if goal=="spk":
        indices = find_initial_clusters( T, k)
        clusters = get_cluster_list( T, indices)
        spkmeans.fit(T, clusters, indices)



