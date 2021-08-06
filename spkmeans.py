import numpy as np
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
        s+= "{:10.4f}".format(num) + ","
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
    assert (k== len(indices))
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
    if len(sys.argv)== 4:
        k = int( sys.argv[1])
        goal = str(sys.argv[2])
        filename = str(sys.argv[3])
    else:
        assert False,"Number of arguments is wrong"

    check_param(k, goal, filename)
    return k, goal, filename


def save_to_out_out(T,flag,filename):
    i = filename.split("_")
    print("this is i",i)
    i = i[2].split(".")
    print(i)
    i = int(i[0])
    result_name = f'Test_files/output_{i}_{flag}.txt'
    A = np.array(T)
    A = np.round(A,4)
    A = A.tolist()

    with open(result_name, 'w') as filehandle:
        for listitem in A:
            row = str(listitem)[1:-1]
            row = row.replace(" ", "")
            print(row)
            filehandle.write('%s\n' % row)



def main():
    T0 = time.process_time()
    np.random.seed(0)
    k, goal, filename = parse2()
    observations = load_data_to_dots(filename)
    assert len(observations) > k, "k must be smaller than number of input vectors"
    T, k = spkmeans.get_flag(goal, k, observations)





    indices = find_initial_clusters(T, k)
    clusters = get_cluster_list(T, indices)
    t1 = time.process_time()
    print(f"K:{k}")
    spkmeans.fit(T, clusters, indices, observations)
    print("done !")
    return

main()