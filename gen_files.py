import random
import os
import numpy as np
import time
N = 20
D = 9
MID = 400

ITER = 5

np.random.seed(42)
random.seed(42)
os.system("python  setup.py build_ext --inplace")

def is_zero_EPS(num,eps=0.0000005):
    return num < eps and -num > -eps


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
    print(f"n= {n}, d={d}",end=" ")


    with open(filename, 'w') as filehandle:
        for listitem in A:
            row = str(listitem)[1:-1]
            row = row.replace(" ", "")
            filehandle.write('%s\n' % row)

def create_output():
    flags = ["wam","ddg","lnorm","spk"]
    for i in range(ITER):
        filename = f'Test_files/input_{i}.txt'
        create_test_file(filename, False)

        for flag in flags:

            result_name = f'Test_files/output_{i}_{flag}.txt'

            cmd = f"python spkmeans.py 0 {flag} {filename} > {result_name} "
            print(flag,end =" ")
            start = time.time()
            os.system(cmd)
            end = time.time()
            dif = round(end- start,2)
            print(f"time: {dif }")


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

def load_data_to_dots(filename):
    dots_list = []
    file = open(filename, 'r')
    Lines = file.readlines()

    for line in Lines:
        if "," not in line:
            continue
        dot = [float(word) for word in line.split(sep=",")]
        dots_list.append(dot)
    return dots_list

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

#create_output()

create_output_symatric()
