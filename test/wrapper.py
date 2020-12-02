import time
import ctypes
import numpy as np
import astropy.time as att
from scipy.io import FortranFile

libname = "libfort.so"
tsize = 300000

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
    import pandas as pd
    from reentry_py import sattelite_trj_calc as sat
    from reentry_py.calc_all_range import vel_eci2ecef
    from reentry_py import calcutil as util

    tsize = 10**7
    # n = 3
    # A = np.array([[1,2,3],[4,5,6],[7,8,9]], dtype=np.float64)
    # f = FortranFile("input.dat", "w")
    # f.write_record(A)

    # print("** write from python")
    # print(n)
    # print(A)

    # B = call_fort(n, (A))

    deb    = sat.trj_calc_fromTLE(sat.l1, sat.l2)
    pos    =  []
    vs     = []
    pos_2D = []
    r0     = deb.get_position()
    v0     = deb.get_velocity_vec()
    _pos   = list(r0*1000)
    stdate = att.Time(deb.epoch, format="mjd")
    r0     = np.array(sat.pm.eci2ecef(*_pos, stdate))/1000
    v0     = vel_eci2ecef(v0, stdate)
    rotmat = util.euler_angle([deb.raan, deb.inc, deb.AoP],[3,1,3])
    rotmat_T = rotmat.T

    # posi  = np.array([-6637.45552308, 1418.74612805, 7.33199191])
    posi  = r0.copy()
    wei   = 100.0
    ref   = np.pi
    iniR  = np.linalg.norm(posi)
    iniV  = np.linalg.norm(v0)
    iniR  = np.linalg.norm(r0)
    iniV  = np.linalg.norm(v0)
    gamma = np.pi/2.0 - np.arccos((r0.dot(v0))/(iniR*iniV))
    # gamma = 0
    print("r0 = ", posi, iniR)
    print("v0 = ", v0, iniV)
    print("gamma=",util.rad2deg(gamma))
    # input()
    psi   = deb.inc
    dt    = 1e-2
    ts    = np.arange(dt*tsize, step=dt)
    F_t   = np.zeros((tsize, 3))
    F_t[:,0] = -0.001
    F_t[:,1] = -0.01
    # print(F_t)
    
    sttime = time.time()
    res = call_reentry(posi, wei, ref, iniV, gamma, psi, tsize, ts, F_t)
    print("R(py) = ", iniR)
    print("gamma=",util.rad2deg(gamma))
    print("psi=",util.rad2deg(psi))
    print("calc trj time: ", time.time()-sttime)

    # print(res)
    # print(B)

    # f = FortranFile("output.dat", "r")
    # dat = ((f.read_reals(dtype=np.float64)))
    # print(dat)

    sttime = time.time()
    f = FortranFile("trj_calc_rslt.dat",  "r")
    datsize = 0
    dat = f.read_record(f"3f8")
    f.close()
    pd.DataFrame(dat).to_csv("result.csv")
    print("initial r = ", iniR)
    print("last r = ",dat[-1,1] )
    print("read time: ",time.time()-sttime)
    print("landing time[s]: ", dat[-1, 0])

    