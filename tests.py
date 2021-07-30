import os
os.system("python setup.py build_ext --inplace")


#os.system("python3 spkmeans.py 3 wam dots_10.txt")
#os.system("python3 spkmeans.py 3 ddg dots_10.txt")
#os.system("python3 spkmeans.py 3 lnorm dots_10.txt")
#os.system("python3 spkmeans.py 3 jacobi dots_10.txt")
os.system("python3 spkmeans.py 3 spk dots_10.txt")