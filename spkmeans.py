import numpy as np
import pandas as pd
import argparse
import sys
import time
import math
import spkmeans
#python3 kmeans_pp.py 3 100 input_1_db_1.txt input_1_db_2.txt

## NOTE - WORK with python 3 only

def smart_print(msg):
    print_time = True
    if print_time:
        print(msg)

def load_data_to_dots(filename):
    dots_list = []
    file = open(filename, 'r')
    Lines = file.readlines()

    for line in Lines:
        dot = [float(word) for word in line.split(sep=",")]
        dots_list.append(dot)
    return dots_list

def print_dot_values(dot):
    s=""
    for num in dot:
        s+=str(num)+","
    print(s[:-1])

def print_dot_list_values(dot_lst):
    for dot in dot_lst:
        print_dot_values(dot)



def find_index_nearest_cluster(w):# returns a random index by the chances provided by w
    return int(np.random.choice(range(len(w)), 1, replace=False, p=w))



def get_cluster_distance(dot, cluster):
    a = np.array(cluster)
    b = np.array(dot)
    d = np.linalg.norm(a-b)
    d = math.pow(d, 2)
    return d


def find_initial_clusters(n_dots, k):
    n = len(n_dots)
    random_index = np.random.randint(n-1)
    first_cluster = n_dots[random_index]
    clusters_list = [first_cluster]
    indices = [random_index]

    distances_list = [float("inf")] * n
    for z in range(1, k):
        for j, dot in enumerate(n_dots):
            d = get_cluster_distance(dot, clusters_list[z - 1])
            distances_list[j] = min(distances_list[j], d)
        sum_d = sum(distances_list)
        weights = [d/sum_d for d in distances_list]
        index = find_index_nearest_cluster(weights)
        indices.append(index)
        clusters_list.append(n_dots[index])
    return indices

def get_cluster_list (observations, indices):
    clusters = []
    for ind in indices:
        clusters.append(observations[ind])
    return clusters


def check_param(k, goal,filename):
    assert type(k) is int, "k must be an integer"
    assert k >= 0, "k must be non negative"
    assert type(goal) is str, "goal must be a string"
    assert type(filename) is str, "file path must be a string"

def parse2():
    smart_print(f'Number of arguments:{len(sys.argv)} arguments.')
    smart_print(f'Argument List: {str(sys.argv)}')
    if len(sys.argv)== 4:
        k = int( sys.argv[1])
        goal = str(sys.argv[2])
        filename = str(sys.argv[3])
    else:
        assert False,"Number of arguments is wrong"

    check_param(k, goal, filename)
    return k, goal, filename



def main():
    T0 = time.process_time()
    np.random.seed(0)
    k, goal, filename = parse2()
    observations = load_data_to_dots(filename)
    assert len(observations) > k, "k must be smaller than number of input vectors"

    T1 = time.process_time()
    smart_print(f"time load data:{T1 - T0}")
    T_and_k = spkmeans.get_flag(goal,k,observations)
    T2 = time.process_time()
    print("done !")
    smart_print(f"time for {goal}:{T2 - T1}")
    if goal!="spk":
        return
    indices = find_initial_clusters(T_and_k[0], k)
    clusters = get_cluster_list(observations, indices)
    t1 = time.process_time()
    smart_print(f"time to read files:{time.process_time() - t1}")
    spkmeans.fit(T_and_k[0], clusters, indices, observations)
    return
    clusters = np.array()
    clusters = np.round(clusters, 4)
    smart_print(f"time to kmean files:{time.process_time() - t1}")
    smart_print("-----RESULTS:------")
    print_dot_values(indices)
    print_dot_list_values(clusters)

main()