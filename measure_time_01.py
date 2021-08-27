import random
import os
import numpy as np
import time

os.system("python  setup.py build_ext --inplace")
for flag in ["spk"]:
    cmd = f"python spkmeans.py 0 {flag} Test_files/input_0.txt "
    print(cmd)
    start = time.time()
    os.system(cmd)
    end = time.time()
    print(flag +". time :")
    print(end - start)