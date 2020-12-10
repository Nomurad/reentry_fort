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
    import matplotlib.pyplot as plt
    import numpy as np
    from numpy import sin, cos
    from reentry_py.calcutil import rad2deg, deg2rad
    from reentry_py import sattelite_trj_calc as sat
    from reentry_py.calc_all_range import vel_eci2ecef
    from reentry_py import calcutil as util

    # n = 3
    # A = np.array([[1,2,3],[4,5,6],[7,8,9]], dtype=np.float64)
    # f = FortranFile("input.dat", "w")
    # f.write_record(A)

    # print("** write from python")
    # print(n)
    # print(A)

    # B = call_fort(n, (A))

    # tsize = 10**6

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
    wei   = 1000.0
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
    F_t   = np.zeros((tsize))
    # F_t[:] = -10.0
    # F_t[:,1] = -0.01
    # print(F_t)
    # gamma = -0.0
    
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
    dat = f.read_record(f"11f8")
    f.close()
    print("landing time[s]: ", dat[-1, 0])
    # pd.DataFrame(dat).to_csv("result.csv")
    # print("initial r = ", iniR)
    print("last r = ",dat[-1,1] )
    print("read time: ",time.time()-sttime)


       ########## save ##########
    t_s = dat[:,0]
    r_s = dat[:,1]
    theta_s = dat[:,2]
    phi_s = dat[:,3]
    gamma_s = (dat[:,6])
    v_s = dat[:,5]
    r_dot_s = dat[:,8]
    theta_dot_s = r_s*(dat[:,9])
    phi_dot_s = r_s*(dat[:,10])
    ### set path
    # if not os.path.exists("result"):
    #     os.makedirs("result")
    pos_array2 = []
    for r, theta, phi in zip(r_s, theta_s, phi_s) :
        x = r*cos(theta)*cos(phi)
        y = r*sin(theta)*cos(phi)
        z = r*sin(phi)
        position = [x, y, z]
        pos_array2.append(position)
    pos_array2 = np.array(pos_array2)
    print(F_t)

    
    # savedata = {
    #     "pos_1st_stage": pos_1st_stage,
    #     "pos_2nd_stage": pos_2nd_stage,
    #     "r_s": r_s,
    #     "theta_s": theta_s,
    #     "phi_s": phi_s
    # }

    # for data_key, data in savedata.items():
    #     fname = f"result_reentry_{data_key}.csv"
    #     result_path = os.path.join("result", fname)
    #     np.savetxt(result_path, data, delimiter=",")

    ########## Visualize ##########
    # plot graph
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # earth & trajectory
    r_earth = 6378.0 # set semi-axis
    r = 6378.0 # set semi-axis
    theta_1_0 = np.linspace(0, 2*np.pi, 100) 
    theta_2_0 = np.linspace(0, 2*np.pi, 100) 
    theta_1, theta_2 = np.meshgrid(theta_1_0, theta_2_0) # convert to 2-dim array
    x = np.cos(theta_2)*np.sin(theta_1) * r 
    y = np.sin(theta_2)*np.sin(theta_1) * r 
    z = np.cos(theta_1) * r 
    ax.plot_surface(x,y,z, alpha=0.3)  #plot earth

    pos_2nd_stage = (pos_array2)
    ax.scatter(pos_2nd_stage[:,0], pos_2nd_stage[:,1], pos_2nd_stage[:,2],
               c="red", s=2)
    plt.gca().set_aspect(aspect="auto")
    fig.tight_layout()

    # 2D trajectory
    alpha = np.sqrt(theta_s**2 + phi_s**2)
    downrange = (alpha - alpha[0])*r_earth
    # print(downrange)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title("2D trajectory")
    ax.scatter(downrange, r_s - r_earth)
    ax.set_xlabel("downrange [km]")
    ax.set_ylabel("Height [km]")
    ax.set_ylim(0, ax.get_ylim()[-1])
    fig.tight_layout()

    # fig = plt.figure()
    # plt.scatter(t_s, r_dot_s)
    # fig = plt.figure()
    # plt.scatter(t_s, rad2deg(gamma_s))
    print(gamma_s)
    # plt.scatter(t_s, r_dot_s)
    # plt.scatter(t_s, np.sqrt(theta_dot_s**2+phi_dot_s**2))


    plt.show()
