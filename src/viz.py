import matplotlib.pyplot as plt
import numpy as np

from math import sqrt
from random import seed

from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter

from classes.mpoly import MPoly
from classes.point import Point
from distance.mpoly_mpoly_dist import mpoly_mpoly_cfs
from distance.v_clip import compute_s

fig, ax = plt.subplots()
pln, = ax.plot([], [], '-r.')
rln, = ax.plot([], [], '-b.')
cfln, = ax.plot([], [], '-go')
txt = ax.text(-8, 42, '', fontsize=15)

# For reproducibility
seed(1)
n = 10
mp = MPoly.random(2 * n, 2 * n, 3 * n, 3 * n, n, n, n)
mr = MPoly.random(0, 0, n, n, n, n, n)
mp.start_pose.theta = 0
mp.end_pose.theta = 3.14
mr.start_pose.theta = 0
mr.end_pose.theta = 3.14
cfs, _, _ = mpoly_mpoly_cfs(mp, mr)

intervals = 200
ts = [i / intervals for i in range(intervals + 1)]


def init():
    ax.set_xlim(-n, 5 * n)
    ax.set_ylim(-n, 5 * n)
    return pln, rln, cfln, txt


def find_cf(cfs, t):
    for i, (t_i, cf_p, cf_r) in enumerate(cfs):
        if i == len(cfs) - 1:
            return cf_p, cf_r
        else:
            if t_i <= t < cfs[i + 1][0]:
                return cf_p, cf_r


def get_points(poly):
    xsp = [p.x for p in poly.pts]
    ysp = [p.y for p in poly.pts]
    return xsp, ysp


def update(t):
    p = mp.at(t)
    r = mr.at(t)
    xsp, ysp = get_points(p)
    pln.set_data(xsp, ysp)
    xsr, ysr = get_points(r)
    rln.set_data(xsr, ysr)
    cf_p, cf_r = find_cf(cfs, t)
    if cf_p % 2 == 0 and cf_r % 2 == 0:
        vp = p.get_v(cf_p // 2)
        vr = r.get_v(cf_r // 2)
    elif cf_p % 2 == 0 and cf_r % 2 == 1:
        vp = p.get_v(cf_p // 2)
        vr1 = r.get_v(cf_r // 2)
        vr2 = r.get_v(cf_r // 2 + 1)
        s = compute_s(vp, vr1, vr2)
        vr = Point.interpolate(vr1, vr2, s)
    elif cf_p % 2 == 1 and cf_r % 2 == 0:
        vp1 = p.get_v(cf_p // 2)
        vp2 = p.get_v(cf_p // 2 + 1)
        vr = r.get_v(cf_r // 2)
        s = compute_s(vr, vp1, vp2)
        vp = Point.interpolate(vp1, vp2, s)
    cfln.set_data([vp.x, vr.x], [vp.y, vr.y])
    dist = sqrt(pow(vp.x - vr.x, 2) + pow(vp.y - vr.y, 2))
    txt.set_text(f"T = {round(t, 2)}\nDistance = {round(dist, 1)}")
    return pln, rln, cfln, txt


ani = FuncAnimation(
    fig, update, frames=ts, interval=30, init_func=init, blit=True)
plt.show()

# writergif = PillowWriter(fps=60)
# ani.save("../img/poly_poly_dist.gif", writer=writergif)