import os
os.system("python setup.py build_ext --inplace")
os.system("python3 spkmeans.py 3 wam dots_10.txt")