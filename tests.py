import os
import  numpy as np
import sys
PRINT_DIFF_ELEM = False
np.set_printoptions(threshold=sys.maxsize,precision= 4)

def compre_mats(my_result, test_results ):


    n1 = len(my_result)
    if n1==0:
        print ("EMPTY FILE !")
        return 0
    d1 = len(my_result[0])
    n2 = len(test_results)
    d2 = len(test_results[0])
    if (n1 != n2):
        print("number of rows differ")
        print(f"my: {n1}")
        print(f"test:{n2}")
        return 0
    if (d1 != d2):
        print("number of cols differ")
        print(f"my: {d1}")
        print(f"test:{d2}")
        return 0

    count = 0
    for i in range(n1):
        for j in range(d1):
            if my_result[i][j] != test_results[i][j]:
                    count +=1
                    if PRINT_DIFF_ELEM:
                        print(f"i={i}, j={j}, diffre:")
                        print(f"my: {my_result[i][j]}")
                        print(f"test:{test_results[i][j]}")
   ##  print(my_result)
    ## print("----------------------------------------------------------")
    ## print(test_results)
    if count >0:
        print(f"differ in {count} / {n1*d1} elements")
    else:
        print( "PASS !")

def run_and_compre(flag,test_ind,k):
    i = test_ind


    if flag == "jacobi":
        test_results_file = f"tests/reference/jacobi/test{i}_{flag}_output_P.txt"
        my_result_file = f"tests/reference/my_output/output_{i}_{flag}_output_P.txt"
        folder = "jacobi"
    else:
        test_results_file = f"tests/reference/general/test{i}_{flag}_{k}_output_P.txt"
        my_result_file = f"tests/reference/my_output/output_{i}_{flag}_{k}_output_P.txt"
        folder = "spk"

    cmd = f"python spkmeans.py {k} {flag}   tests/test_data/{folder}_tests/test{i}.csv >{my_result_file} "
    print(cmd)
    os.system(cmd)


    my_result = np.loadtxt(my_result_file, delimiter=',')
    test_results = np.loadtxt(test_results_file, delimiter=',')
    compre_mats(my_result, test_results)


#os.system( "python  setup.py build_ext --inplace")
def test_loop():
    flags =["wam","ddg","lnorm","jacobi"] # TODO check spk
    for i in range(10):
        print(f"=================={i}========================")
        for flag in flags:
                run_and_compre(flag, i, 1)
    print("Done!!!")



test_loop()


def get_spk_tests():
    spk_arg_list = []
    files = os.listdir('tests/reference/general/')
    for file in files:
        if "spk"  not in file:
            continue
        if "_P"  not in file:
            continue
        if file[5:6] =="_":
            i = int( file[4:5])
        else:
            i = int(file[4:6])
        file_sp = file.split("_")
        k=file_sp[2]
        args = (k,i)
        spk_arg_list.append(args)
    return spk_arg_list

def run_C_tests():
    spk_arg_list = get_spk_tests()
    for k,i in spk_arg_list:
        args = f" {k} spk  tests/test_data/spk_tests/test{i}.csv"
        my_result_file = f"tests/reference/my_spk/output_{i}_spk_{k}_C_ELAD.txt"
        cmd = f"spkmeans {args} >{my_result_file} "
        print(cmd)
        os.system(cmd)


#run_C_tests()