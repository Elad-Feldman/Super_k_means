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
def create_test_file(filename):


    n = random.randint(4, N)
    d = random.randint(4, D)
    m =  1 # random.randint(1, MID)


    print(f"n:{n}, d:{d}")
    B = (0.5 - np.random.rand(n, d))
    A= B * m
    # A = np.tril(B) + np.tril(B, -1).T
    A = A.round(4).tolist()

    # create input file
    with open(filename, 'w') as filehandle:
        for listitem in A:
            row = str(listitem)[1:-1]
            row = row.replace(" ", "")
            filehandle.write('%s\n' % row)

def create_output():
    flags = ["wam","ddg","lnorm"]
    for i in range(10):
        filename = f'Test_files/input_{i}.txt'

        create_test_file(filename)
        for flag in flags:
            result_name = f'Test_files/output_{i}_{flag}.txt'
            print
            T0 = time.process_time()
            os.system(f"python spkmeans.py 3 {flag} {filename} > {result_name } ")
            T1 = time.process_time()
            print(f"for {i}  in {flag} , time:{T1-T0}")


create_output()
