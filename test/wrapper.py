import ctypes
import numpy as np
from scipy.io import FortranFile

libname = "libfort.so"

def call_fort(n, A: np.ndarray):
    # コメント
    f = np.ctypeslib.load_library(libname, ".")

    f.test.argtypes = [
        ctypes.POINTER(ctypes.c_int32),
        np.ctypeslib.ndpointer(dtype=np.float64, shape=A.shape)
    ]
    f.test.restype = ctypes.c_void_p
    fn = ctypes.byref(ctypes.c_int32(n))
    B = A.T.copy()
    f.test(fn, B)
    return B

if __name__ == "__main__":
    print("** calling fortran from python.")

    n = 3
    A = np.array([[1,2,3],[4,5,6],[7,8,9]], dtype=np.float64)

    f = FortranFile("input.dat", "w")
    f.write_record(A)

    print("** write from python")
    print(n)
    print(A)

    B = call_fort(n, (A))

    print(B)

    f = FortranFile("output.dat", "r")
    dat = ((f.read_reals(dtype=np.float64)))
    print(dat)