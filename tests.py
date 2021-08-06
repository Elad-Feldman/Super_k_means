import os
import  numpy as np
os.system( "python  setup.py build_ext --inplace")

os.system(f"python spkmeans.py 3 spk dots_10.txt")
#for flag in ["wam","ddg","lnorm","jacobi","spk"]:
 #   print(flag)
  #  os.system(f"python spkmeans.py 3 {flag} dots_10.txt")


