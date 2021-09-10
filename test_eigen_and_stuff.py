
def is_vectors_equals(a,b):
    if len(a) != len(b):
        return False
    for e,v in zip(a,b):
        if abs(e-v) > 0.0001:
            return False
    return True

def compare_eigen_vectors(A,vectors,values):
    n = len(A)
    np_A = np.matrix(A)
    for i in range(n):  # A * v = lambda * v
        A_v = np.asarray(np_A @ (vectors[i])).reshape(-1)
        lambda_v = ( values[i] *vectors[i])
        if not is_vectors_equals(A_v, lambda_v):
            print("Diff eigen Vectors !")
            print_vector(A_v)
            print_vector(lambda_v)

def test_eigen(A,T):
    n = len(A)
    np_val,np_vec = np.linalg.eigh(A)
    # print(np_vec.round(4).transpose())

    # Comapre eigen values
    np_val_s = np.sort(np_val).round(4)
    my_val = np.sort(T[0]).round(4)
    if not is_vectors_equals(np_val_s,my_val):
        print("Diff eigen values !")
        print_vector(np_val_s)
        print_vector(my_val)

    print("start numpy values:")
    compare_eigen_vectors(A,np_vec.transpose(),np_val) # compre numpy values, sainty check
    print("Done numpy values ! ")
    print("start our values:")
    compare_eigen_vectors(A, np.array(T[1:]), T[0])  # compre numpy values, sainty check // CHECK SORT !
    print("Done our values ! ")