
import numpy as np
import math
import  json
import time

N = 5
d=5
np.random.seed(0)
np.set_printoptions(precision=4)
ind = 0
failed_tests = dict()

def mat_to_dic(mat1,mat2):
    N = mat1.shape[0]
    mat1 = np.round(mat1,decimals=4)
    mat1.tolist()
    mat2 = np.round(mat1, decimals=4)
    mat2.tolist()
    res = []
    for i in range(N):
        listToStr1 ="[" + ','.join([str(elem) for elem in mat1[i]]) +"] "
        listToStr2 = "[" + ','.join([str(elem) for elem in mat2[i]]) + "]"
        res.append(listToStr1 +listToStr2)

    return res

def sign(x):
    if x==0:
        return 1
    return np.sign(x)

def no_dignal(x):
    N=x.shape[0]
    inf_diag = np.matrix(np.eye(N) * 50000)
    return x - inf_diag

def is_dignal(x):
    c = np.count_nonzero(x - np.diag(np.diagonal(x)))
    return c==0


def create_weigh_adj_mat_2():
    N = 100
    d = 60
    mat = np.random.rand(N, d)
    w = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            norm = np.linalg.norm(mat[i]-mat[j])
            w[i][j] =math.exp( - norm / 2 )

    return w


def create_weigh_adj_mat_3(mat):
    N = 5
    d = 3
    w = np.zeros((N,N))
    for i in range(N):
        for j in range(i):
            norm = np.linalg.norm(mat[i]-mat[j])
            w[i][j] =math.exp( - norm / 2 )
            w[j][i] =w[i][j]
    return w

def create_weigh_adj_mat():
    N = 100
    d = 60
    mat = np.random.rand(N, d)
    w = np.zeros((N,N))
    for i in range(N):
        for j in range(i):
            norm = np.linalg.norm(mat[i]-mat[j])
            w[i][j] =math.exp( - norm / 2 )
            w[j][i] =w[i][j]
    return w

def create_dig_dgr_mat(w):
    N = w.shape[0]
    D = np.zeros((N, N))
    for i in range(N):
        d = sum(w[i])
        D[i][i] = 1 / math.sqrt(d)
    return D


def test_function(func,mat_start,mat_result):
    global ind, failed_tests
    ind += 1
    func_res = func(mat_start)
    if   not np.array_equal(func_res,mat_result):
        print("test fail:"+str(ind))
        failed_tests[ind] =(mat_to_dic(func_res,mat_result))

    else:
        print("test pass:" + str(ind))


def run_tests():
    for i in range(8):
        N = 5
        d = 3
        a = np.random.rand(N, d)
        w = create_weigh_adj_mat_3(a)
        D = create_dig_dgr_mat(w)
        test_function(create_weigh_adj_mat_3,a,w)
        test_function(create_weigh_adj_mat_3, a, D)
    with open('data.json', 'w') as fp:
        json.dump(failed_tests, fp,  indent=1)


def meas_avg_time(string, reapet=5):
    s = 0
    for i in range(reapet):
        t0 = time.time()
        eval(string)
        t1 = time.time()
        s+= t1-t0
    print(f"avg time: {s/reapet}")


def test_time():
    meas_avg_time("create_weigh_adj_mat_2()")
    meas_avg_time("create_weigh_adj_mat()")

run_tests()





