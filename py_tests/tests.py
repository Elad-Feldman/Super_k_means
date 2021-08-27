import os
import  numpy as np
os.system( "python  setup.py build_ext --inplace")
cmd = f"python spkmeans.py 0 jacobi   Test_files/input_J_0.txt"
print(cmd)
os.system(cmd)
def test_loop():
    for i in range(10):
        print(f"=================={i}========================")
        cmd = f"python spkmeans.py 0 jacobi   Test_files/input_J_{i}.txt > Test_files/OutPut/output_J_{i}.txt "
        print(cmd)
        os.system(cmd)

    print("Done!!!")


def test_Yair():
    os.system(f"python spkmeans.py 3 lnorm YAIR_TEST/test2.csv > YAIR_TEST/MY_res_lnorm_2.txt")


#test_loop()