import random
import os
import numpy as np
import time
N = 900
D = 5
MID = 400

TEST_NUM = 3

np.random.seed(42)
random.seed(42)
os.system("python  setup.py build_ext --inplace")

def print_if_files_are_diff(fn1, fn2):
    msg = f"NOT STABLE! {fn1} and {fn2} differ at "
    l1 = load_to_list(fn1)
    l2 = load_to_list(fn2)

    if len(l1) != len(l2):
        print( msg + " length")
        return False

    for i, (a, b) in enumerate(zip( l1, l2)):
        if a != b:
            print(msg + f" line {i}")
            return False
    return True

def check_test_stabilty( filename, flag, k,ITER = 1):
    Total_time = 0
    for i in range(TEST_NUM):
        curr_name = f"_{i}_.txt"
        cmd = f"python spkmeans.py {k} {flag} {filename}  > {curr_name}"
        start = time.time()
        os.system(cmd)
        end = time.time()
        Total_time += end - start
        if i==0:
            continue
        prev_name = f"_{i-1}_.txt"
        print_if_files_are_diff(curr_name,prev_name)
        os.remove(prev_name)
    print ("\t{:10.2f} [sec]".format(Total_time / ITER))




def print_vector(dot):
    s=""
    for num in dot:
        s+= "{:10.4f}".format(num) + ","
    print(s[:-1])

def print_matrix(mat):
    for row in mat:
        print_vector(row)

def create_matrix_4digits(n,d):
    B = (0.5 - np.random.rand(n, d))
    A = B
    if n == d:
         A = np.tril(B) + np.tril(B, -1).T
    return A.round(4).tolist()

def create_test_file(filename,is_symatric):
    n = random.randint(4, N)
    d = random.randint(4, D)
    d =  n if is_symatric  else d
    A = create_matrix_4digits(n,d)
    print(f"\tn= {n}, d={d}")


    with open(filename, 'w') as filehandle:
        for listitem in A:
            row = str(listitem)[1:-1]
            row = row.replace(" ", "")
            filehandle.write('%s\n' % row)

def create_output():
    print(f"N={N}, D={D}")
    flags = ["wam","ddg","lnorm","spk"]
    for i in range(TEST_NUM):
        filename = f'Test_files/input_{i}.txt'
        print(f"input_{i}.txt: ")
        create_test_file(filename, False)
        k = 3
        print(f"\tFlag\tk \t\tavg time")
        for flag in flags:
            print(f"\t {flag} \t{k}",end=" ")
            check_test_stabilty(filename,flag,k)


def create_output_symatric():
    flag ="jacobi"
    for i in range(10):
        filename = f'Test_files/input_J_{i}.txt'
        create_test_file(filename, True)
        result_name = f'Test_files/output_J_{i}.txt'
        cmd = f"python spkmeans.py 3 {flag} {filename} > {result_name} "
        print(flag, end=" ")
        start = time.time()
        os.system(cmd)
        end = time.time()
        dif = round(end - start, 2)
        print(f"time: {dif}")


def compre_Mat(A,B):
    if len(A) != len(B):
        return False
    for row_a,row_b in zip(A,B):
        if len(row_a) != len(row_b):
            return False

        for elm_a, elm_b in zip(row_a, row_b):
            if elm_a != elm_b:
                return False
    return True

def load_to_list(filename):
    with open(filename, 'r')as file:
        lines = file.readlines()
    return lines


def test_3_basic(is_py):
    for i in range(ITER):
        filename = f'Test_files/input_{i}.txt'
        for flag in ["spk"]:
            if is_py:
                cmd = f"python spkmeans.py 1 {flag}  {filename} "
            else:
                os.system("gcc spkmeans.c && gcc  -o spkmeans spkmeans.c")
                cmd = f"spkmeans 1 {flag}  {filename} "
            os.system(cmd + "> tmp.txt")
            dots = load_data_to_dots("tmp.txt")
            for j in range(3):
                os.system(cmd + "> tmp.txt")
                dots_i = load_data_to_dots("tmp.txt")
                STATUS = "FAIL"
                print_matrix(dots)
                print("=====================")
                print_matrix(dots_i)
                if  compre_Mat(dots,dots_i):
                    STATUS = "PASS"
                print(f"{STATUS}:  flag: {flag}, {i}: {j}")

#test_3_basic(True)
#create_output_symatric()

create_output()



