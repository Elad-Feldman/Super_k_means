import random
import os
import numpy as np
import time
N = 12
D = 9
MID = 400

np.random.seed(42)
random.seed(42)
os.system("python  setup.py build_ext --inplace")
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


    with open(filename, 'w') as filehandle:
        for listitem in A:
            row = str(listitem)[1:-1]
            row = row.replace(" ", "")
            filehandle.write('%s\n' % row)

def create_output():
    flags = ["wam","ddg","lnorm","spk"]
    for i in range(10):
        filename = f'Test_files/input_{i}.txt'
        create_test_file(filename, False)

        for flag in flags:

            result_name = f'Test_files/output_{i}_{flag}.txt'
            cmd = f"python spkmeans.py 3 {flag} {filename} > {result_name} "
            print(cmd)


def create_output_symatric():
    flag ="jacobi"
    for i in range(10):
        filename = f'Test_files/input_J_{i}.txt'
        create_test_file(filename, True)
        result_name = f'Test_files/output_J_{i}.txt'
        cmd = f"python spkmeans.py 3 {flag} {filename} > {result_name} "
        print(cmd)
        os.system(cmd)


create_output()
create_output_symatric()
