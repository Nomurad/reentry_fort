import time
import ctypes
import numpy as np
from scipy.io import FortranFile

libname = "libfort.so"
tsize = 100000

def call_fort(n, A: np.ndarray):
    # コメント
    f = np.ctypeslib.load_library(libname, ".")

    f.test.argtypes = [
        ctypes.POINTER(ctypes.c_int32),
        np.ctypeslib.ndpointer(dtype=np.float64, shape=A.shape)
    ]
    f.test.restype = ctypes.c_void_p
    fn = ctypes.byref(ctypes.c_int32(n))
    # B = A.T.copy()
    B = np.asanyarray(A, order="F")  # transpose for fortran code
    f.test(fn, B)
    return B

class trj_results(ctypes.Structure):
    _fields_ = [
        ("r_hist"        , ctypes.c_double*tsize),
        ("theta_hist"    , ctypes.c_double*tsize),
        ("phi_hist"      , ctypes.c_double*tsize),
        ("rho_hist"      , ctypes.c_double*tsize),
        ("V_hist"        , ctypes.c_double*tsize),
        ("gamma_hist"    , ctypes.c_double*tsize),
        ("psi_hist"      , ctypes.c_double*tsize),
        ("theta_dot_hist", ctypes.c_double*tsize),
        ("phi_dot_hist"  , ctypes.c_double*tsize)
    ]

def call_reentry(posi, weight, ref_area, initV, gamma, psi, len_t, ts, F_t):
    f = np.ctypeslib.load_library("libreentry_sub.so",".")
    f.reentry_from_py_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        
        ctypes.POINTER(ctypes.c_int32),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64)
        # ctypes.POINTER(trj_results)
        # ctypes.POINTER(ctypes.c_double)
    ]
    f.reentry_from_py_.restype = ctypes.c_void_p
    pweight   = ctypes.byref(ctypes.c_double(weight))
    pref_area = ctypes.byref(ctypes.c_double(ref_area)) 
    pinitV    = ctypes.byref(ctypes.c_double(initV)) 
    pgamma    = ctypes.byref(ctypes.c_double(gamma)) 
    ppsi      = ctypes.byref(ctypes.c_double(psi)) 
    plen_t    = ctypes.byref(ctypes.c_int(len_t))
    # res = trj_results()
    res = ctypes.byref(ctypes.c_double(1.0))
    # print(res._fields_)

    posi = np.asanyarray(posi, order="F")
    ts = np.asanyarray(ts, order="F")
    F_t = np.asanyarray(F_t, order="F")
    f.reentry_from_py_(posi, pweight, pref_area, pinitV, pgamma, ppsi, plen_t, ts, F_t)

    return res


if __name__ == "__main__":

    # n = 3
    # A = np.array([[1,2,3],[4,5,6],[7,8,9]], dtype=np.float64)
    # f = FortranFile("input.dat", "w")
    # f.write_record(A)

    # print("** write from python")
    # print(n)
    # print(A)

    # B = call_fort(n, (A))

    posi = np.array([6378.0+400.0, 0.0, 0.0])
    wei = 100.0
    ref = np.pi
    iniV = 7.74
    gamma = 0
    psi = 0
    dt = 1e-2
    ts = np.arange(dt*tsize, step=dt)
    F_t = np.zeros((3, tsize))
    F_t[0,:] = 1e-12
    F_t[1,:] = -1e0
    
    sttime = time.time()
    res = call_reentry(posi, wei, ref, iniV, gamma, psi, tsize, ts, F_t)
    print("calc trj time: ", time.time()-sttime)

    # print(res)
    # print(B)

    # f = FortranFile("output.dat", "r")
    # dat = ((f.read_reals(dtype=np.float64)))
    # print(dat)

    sttime = time.time()
    f = FortranFile("trj_calc_rslt.dat",  "r")
    datsize = 0
    dat = f.read_record(f"2f8")
    f.close()
    print("read time: ",time.time()-sttime)
    print("landing time[s]: ", dat[-1, 0])

    