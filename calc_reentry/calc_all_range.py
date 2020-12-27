import time
import sys
import logging
import numpy as np
from numpy import sin, cos 
from numpy import linalg as LA
from scipy.io import FortranFile
from datetime import datetime
import astropy.time as att
from astropy.utils import iers
from astropy.utils.data import download_file
from astropy.utils.iers import IERS_A_URL, IERS_A_URL_MIRROR
iers.conf.iers_auto_url = "ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all"
iers.conf.iers_auto_url = IERS_A_URL_MIRROR
# download_file(iers.conf.iers_auto_url)
# iers.iers.download_file(iers.conf.iers_auto_url, "update")
# iers.conf.auto_download = False # Cancel Auto_Downloading for working offline 
# from pymap3d.eci import eci2ecef

from reentry_py.calcutil import *
from reentry_py.spacecraft import SpaceCraft
from reentry_py.reentry_calc import ReEntry_calc
import reentry_py.sattelite_trj_calc as h_trj
from reentry_py.sattelite_trj_calc import DAY2SEC, SEC2DAY

import wrapper

# TLE DATA
#ISS
l1 = "1 25544U 98067A   13158.50723059  .00016717  00000-0  10270-3 0  9040"
l2 = "2 25544  51.6453 246.3482 0010707  41.5842 318.6122 15.50690189 33209"
# l2 = "2 25544  51.6453 246.3482 7130000  41.5842 318.6122 15.50690189 33209"
#MERIDIAN 6
# l2 = "2 38995 062.9125 153.0958 7134453 295.0386 009.2449 02.00613016  4036"
# l1 = "1 38995U 12063A   13154.54631441  .00000267  00000-0  00000+0 0  1242"

# H-2A R/B
# l1 = "1 41341U 16012E   19008.20966682  .00000428  00000-0  22403-4 0  9999"
# l2 = "2 41341  30.6275  17.9986 0018131 128.9892 231.2355 15.04728935159089"

def _rottrip(angle: np.ndarray) -> np.ndarray:
    angle = angle.squeeze()
    if angle.size > 1:
        raise ValueError("only one angle allowed at a time")
    """ported from:
    https://github.com/dinkelk/astrodynamics/blob/master/rot3.m
    """
    return np.array([[np.cos(angle), np.sin(angle), 0],
                     [-np.sin(angle), np.cos(angle), 0],
                     [0, 0, 1]])

def vel_eci2ecef(vel_eci: np.ndarray, 
                 time: datetime) -> np.ndarray:
    """convert velocity(eci) to velocity(ecef)
    
    Arguments:
        vel_eci {np.ndarray} -- velocity seen from eci location
        time {datetime.datetime} -- datetime
    
    Returns:
        vel_ecef {np.ndarray} -- velocity seen from ecef location
    """
    gst = att.Time(time).sidereal_time("apparent", "greenwich").radian
    gst = np.atleast_1d(gst)
    assert gst.ndim == 1 and isinstance(gst[0], float) # must be in radians
    
    vel_eci = np.atleast_2d(vel_eci)
    assert vel_eci.shape[0] == gst.size, 'length of time does not match number of ECI positions'

    N, trip = vel_eci.shape 
    if vel_eci.ndim > 2 or trip != 3:
        raise ValueError("eci triplets must be shape (N,3)")

    vel_ecef = np.empty_like(vel_eci)

    for i in range(N):
        vel_ecef[i, :] = _rottrip(gst[i]) @ vel_eci[i, :]

    return vel_ecef.squeeze()


def calc_all_range(l1: str, 
                   l2: str,
                   delta_day: float,
                   delta_V: float,
                   angle: float,
                   **kwargs
                   ) -> (np.ndarray, np.ndarray, dict, float):
    """[summary]
    
    Arguments:
        l1 {str} -- TLE string line1
        l2 {str} -- TLE string line2
        delta_day {float} -- Days elapsed from TLE day
        delta_V {float} -- dV
        angle {float} -- dV angle
    
    Returns:
        pos_array -- re-entry trajectory above re-entry height
        pos_array2 -- re-entry trajectory under re-entry height
        ret_dict -- results(r, theta, phi, etc...)
        total_flighttime -- total flight time[sec]
    """

    lg = logging.getLogger(name=__name__)
    verbose = kwargs.get("verbose", logging.WARNING)
    tmpfile = "tmp.log"
    tmpfstream = open(tmpfile, "w")
    fstream = kwargs.get("logfile", None)
    if fstream is None:
        fstream = tmpfstream

    ########## init trj_calc ##########
    deb = h_trj.trj_calc_fromTLE(l1, l2)
    r0 = deb.get_position()
    v0 = deb.get_velocity_vec()

    reentry_height = 122.0
    # delta_day = 0.5
    stay_time = (1/deb.MM)*delta_day

    num_of_trials = 1
    iter_max = 1000*3
    # iter_max = 1
    pos_array = []
    v_array = []
    # pos_2D = []

    ########## calc trajectory above re-entry height ##########
    deb.move(stay_time, "dt")
    reentry_st_date = deb.epoch
    vel = deb.get_velocity_vec()

    # deceleration = 0.9
    # dv = (1.0 - deceleration)*v0
    # v = deceleration*v0
    rotmat = euler_angle([deb.raan, deb.inc, deb.AoP],[3,1,3]) #ijk -> pqw
    rotmat_T = rotmat.T

    v_pqw = deb.get_velocity_vec(_2D=True)
    # v_pqw_direction = v_pqw/LA.norm(v_pqw)
    # dv_pqw = delta_V*v_pqw_direction
    dv_pqw = rot_matrix(angle, 3).dot((v_pqw*delta_V))
    dv = rotmat_T.dot(dv_pqw)
    v = vel + dv
    # print("v0", v0)
    # print("dv = ", dv)
    # print("dv norm = ", LA.norm(dv))

    r = deb.get_position(_2D=False)
    k_ele, _ = deb.get_orbital_element(r, v)
    if k_ele[1] > 1.0:
        return False, False, False, False
    deb.param_set(list(k_ele))

    flag_sccess_reentry = False

    # deceleration = 0.9
    # dv += (1.0 - deceleration)*v
    # v = deceleration*v
    # r = deb.get_position()
    # k_ele, _ = deb.get_orbital_element(r, v)
    # deb.param_set(list(k_ele))

    v_norm = LA.norm(v)
    v0_norm = LA.norm(v0)
    # print(k_ele)
    # print("V0,V :",v0_norm,v_norm , "dv :", dv, file=fstream)
    def _logprint(s: str, end="\n", file=sys.stdout):
        print(s, end=end, file=file)

    logprint = _logprint

    for i in range(iter_max):
        # deb.move((1/deb.MM)/iter_max, "dt")
        deb.move((deb.T)*SEC2DAY/iter_max, "dt_T")
        position = deb.get_position()
        velocity = deb.get_velocity_vec()
        pos_array.append(position)
        v_array.append(velocity)
        r = LA.norm(position)
        v = LA.norm(velocity)
        # print("pos = ", position)
        logprint(f"r(1st stage) = {r:.3f}, v={v:.3f} time={(deb.T)*SEC2DAY/iter_max*i}", 
                end="\r", file=fstream)
        # if velocity[0] > 0:
        #     continue
        if r <= (ReEntry_calc.r_earth + reentry_height) and r > ReEntry_calc.r_earth:
            date = att.Time(deb.epoch, format="mjd").to_datetime()
            print("re entry start height :",LA.norm(position), file=fstream)
            print("epoch: ", date, file=fstream)
            reentry_st_date  = att.Time(deb.epoch, format="mjd")
            _pos = list(position*1000)
            pos_ecef = np.array(h_trj.pm.eci2ecef(*_pos, reentry_st_date))/1000

            flag_sccess_reentry = True
            break
        # print("r = ", LA.norm(position))

    # if r > (ReEntry_calc.r_earth + reentry_height) or deb.ecc >= 1.0:
    # flag_sccess_reentry = True
    if not flag_sccess_reentry:
        print("\nfailed reentry, r(last), e = ", r, deb.ecc)
        return min([LA.norm(pos) for pos in pos_array]), False, False, False

    reentry_st_date = att.Time(deb.epoch, format="mjd")
    # reentry_st_date = reentry_st_date.to_datetime()
    position = deb.get_position()
    _pos = list(position*1000)
    pos_ecef = np.array(h_trj.pm.eci2ecef(*_pos, reentry_st_date))/1000
    # pos_array.append(position)
    # v_array.append(deb.get_velocity_vec())

    pos_array = np.array(pos_array)
    v_array = np.array(v_array)
    try:
        date = att.Time(deb.epoch, format="mjd").to_datetime()
    except:
        print("epoch", deb.epoch)
    flighttime_1st = date - deb.Date


    ########## init re-entry ##########
    initial_state = np.array([ReEntry_calc.r_earth+120, 0, 0])
    craft = SpaceCraft(initial_state, 100, np.pi)
    re_entry = ReEntry_calc(craft, verbose=verbose)
    # logging.basicConfig(level=kwargs.get("verbose", 100))
    # re_entry.craft.C_D = 1.0
    x, y, z = pos_ecef
    # r = LA.norm(pos_ecef)
    # theta = np.arccos(x/LA.norm([x,y]))
    # phi = np.arcsin(z/r)
    re_entry.set_r_vec(pos_ecef)
    r_E = re_entry.r 
    theta_E = re_entry.theta
    phi_E = re_entry.phi
    v_E_eci = v_array[-1] # velocity seen from eci location
    v_E = vel_eci2ecef(v_E_eci, reentry_st_date.to_datetime())
    v = LA.norm(v_E)
    vx, vy, vz = v_E
    print("v_xyz_eci = ", v_E_eci, v, file=fstream)
    print("v_xyz = ", v_E, v, file=fstream)

    rtp    = [r_E, theta_E, phi_E]
    # vr_xy  = vx*cos(rtp[1]) + vy*sin(rtp[1])
    # vr_zx  = vx*cos(rtp[2]) + vz*sin(rtp[2])
    # vtheta =  np.sqrt((vx**2 + vy**2 - vr_xy**2))
    # vphi   =  np.sqrt((vx**2 + vz**2 - vr_zx**2))
    # vr     =  np.sqrt(v**2 - vtheta**2 - vphi**2)
    v_E = rot_matrix(phi_E, "x").dot(rot_matrix(theta_E,"z").dot(v_E))
    # v_E = [vr, vtheta, vphi]
    print(r_E, np.rad2deg(theta_E), np.rad2deg(phi_E), file=fstream)
    print("v_rtp = ", v_E[:], LA.norm(v_E), file=fstream)
    # v_E = v_E_eci
    re_entry.set_v_vec(v_E)
    # re_entry.gamma = deg2rad(-15)
    re_entry.psi = deb.inc
    print("\nr =", pos_ecef, file=fstream)
    print("v =", v_E, file=fstream)

    ### init iteration
    x0 = np.array([r_E, theta_E, phi_E])
    time_max = 10000
    iter_max = time_max*10
    time_list = np.linspace(0, time_max, iter_max+1, endpoint=True)

    ### calc re-entry trajectory 
    print("dt = ", time_list[1]-time_list[0], file=fstream)
    print("gamma = ", rad2deg(re_entry.gamma), file=fstream)
    print("V = ", np.linalg.norm(v_E), file=fstream)
    F_t = np.zeros((3,time_list.shape[0]))
    # Ft = np.zeros((time_list.shape[0]))
    # Ft[:20000] = -( 0.8)
    sttime = time.time()
    res = wrapper.call_reentry(x0, craft.weight, craft.ref_area, 
                            LA.norm(v_E), re_entry.gamma, re_entry.psi, 
                            len(time_list), time_list, F_t)
    # ret = re_entry.trj_calc_3d(x0, time_list, F_t)
    print("calctime=", time.time()-sttime, file=fstream)

    f = FortranFile("trj_calc_rslt.dat",  "r")
    datsize = 0
    dat = f.read_record(f"11f8")
    f.close()



    ########## post processing ##########
    ret_dict = {}
    labels = ["r", "theta", "phi", 
              "rhos", "V", "gamma", "psi",
              "theta_dot", "phi_dot"]
    for item, label in zip(dat.T[1:], labels):
    # for item, label in zip(ret, labels):
        if not hasattr((item), "__iter__"):
            item = [item]
        ret_dict[label]= item

    # print("ret_dict", ret_dict)
    r_s = ret_dict["r"]
    theta_s = ret_dict["theta"]
    phi_s = ret_dict["phi"]
    # r_s = dat[:,1]
    # theta_s = dat[:,2]
    # phi_s = dat[:,3]

    end_idx = (len(r_s))
    if end_idx >= iter_max:
        end_idx -= 1
    flighttime_2nd = time_list[end_idx]
    # flighttime_2nd = dat[-1, 0]
    print(f"re-entry time(2nd) = {flighttime_2nd} [s]")


    pos_array2 = []
    for r, theta, phi in zip(r_s, theta_s, phi_s) :
        x = r*cos(theta)*cos(phi)
        y = r*sin(theta)*cos(phi)
        z = r*sin(phi)
        position = [x, y, z]
        pos_array2.append(position)

    pos_array2 = np.array(pos_array2)
    if r_s[-1] > re_entry.r_earth:
        return pos_array, pos_array2, ret_dict, False

    total_flighttime = flighttime_1st.total_seconds() + flighttime_2nd
    print(f"re-entry time(total) = {total_flighttime} [s]", file=fstream)
    # print(f"re-entry time = {total_flighttime/60} [min]")

    return pos_array, pos_array2, ret_dict, total_flighttime


if __name__ == "__main__":
    import os, shutil
    import random
    import matplotlib.pyplot as plt 
    from matplotlib.pyplot import rc
    from mpl_toolkits.mplot3d import Axes3D

    verbose = logging.DEBUG
    logs = None
    logs = sys.stdout

    ########## calculate Debris's trajectory ##########
    dv = 0.0
    delta_day = 0
    dv = random.random()*0.3
    dv = 0.2
    print("dv", dv)
    iter_max = 1
    for i in range(iter_max):
        angle = deg2rad(80)
        res = calc_all_range(l1, l2, delta_day, dv, angle, 
                             logfile=logs, verbose=verbose)
        pos_1st_stage, pos_2nd_stage, ret_dict, totaltime = res
        if res[-1] is not False:
            break
        elif i >= iter_max - 1:
            # print(res)
            print("not conversion")
            # exit()
        else:
            print("add dv")

    # print("v = ", ret_dict["V"][0])
    # print("r = ", ret_dict["r"][-1])
    # print("theta = ", ret_dict["theta"][-1])
    # print("phi = ", ret_dict["phi"][-1])

    ########## save ##########
    r_s = np.array(ret_dict["r"])
    theta_s = np.array(ret_dict["theta"])
    phi_s = np.array(ret_dict["phi"])
    ### set path
    if not os.path.exists("result"):
        os.makedirs("result")

    
    savedata = {
        "pos_1st_stage": pos_1st_stage,
        "pos_2nd_stage": pos_2nd_stage,
        "r_s": r_s,
        "theta_s": theta_s,
        "phi_s": phi_s
    }

    for data_key, data in savedata.items():
        fname = f"result_reentry_{data_key}.csv"
        result_path = os.path.join("result", fname)
        np.savetxt(result_path, data, delimiter=",")

    ########## Visualize ##########
    # matplotlib init
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams['axes.formatter.use_mathtext'] = True
    rc('mathtext', **{
        'rm': 'Times New Roman',
        'it': 'Times New Roman',
        'bf': 'Times New Roman',
        'fontset': 'stix'
        })
    plt.rcParams["lines.linewidth"] = 3
    plt.rcParams["font.size"] = 25

    # plot graph
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # earth & trajectory
    r = ReEntry_calc.r_earth # set semi-axis
    theta_1_0 = np.linspace(0, 2*np.pi, 100) 
    theta_2_0 = np.linspace(0, 2*np.pi, 100) 
    theta_1, theta_2 = np.meshgrid(theta_1_0, theta_2_0) # convert to 2-dim array
    x = np.cos(theta_2)*np.sin(theta_1) * r 
    y = np.sin(theta_2)*np.sin(theta_1) * r 
    z = np.cos(theta_1) * r 
    ax.plot_surface(x,y,z, alpha=0.3)  #plot earth

    ax.scatter(pos_1st_stage[0,0], pos_1st_stage[0,1], pos_1st_stage[0,2],c="green")
    ax.scatter(pos_1st_stage[:,0], pos_1st_stage[:,1], pos_1st_stage[:,2],
               c="red", s=2)
    # ax.scatter(pos_2nd_stage[:,0], pos_2nd_stage[:,1], pos_2nd_stage[:,2],
    #            c="red", s=2)
    plt.gca().set_aspect(aspect="auto")
    fig.tight_layout()

    # 2D trajectory
    alpha = np.sqrt(theta_s**2 + phi_s**2)
    downrange = (alpha - alpha[0])*ReEntry_calc.r_earth
    # print(downrange)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title("2D trajectory")
    ax.plot(downrange, r_s - ReEntry_calc.r_earth)
    ax.set_xlabel("downrange [km]")
    ax.set_ylabel("Height [km]")
    ax.set_ylim(0, ax.get_ylim()[-1])
    fig.tight_layout()


    plt.show()
