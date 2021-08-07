import os
import  numpy as np
os.system( "python  setup.py build_ext --inplace")
for i in range(1):
    print(f"=================={i}========================")
    os.system(f"python spkmeans.py 0 jacobi Test_files/input_J_{i}.txt")

print("done")


