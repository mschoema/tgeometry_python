import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
params = {
    'legend.fontsize': 25,
    'axes.labelsize': 25,
    'axes.titlesize': 25,
    'xtick.labelsize': 25,
    'ytick.labelsize': 25
}
pylab.rcParams.update(params)

from math import pi

from classes.mpoint import MPoint
from classes.mpoly import MPoly
from classes.point import Point
from classes.poly import Poly
from classes.pose import Pose

from distance.mpoint_mpoly_dist import mpoint_mpoly_cfs
from distance.mpoint_mpoly_dist import mpoint_mpoly_verify_cfs
from distance.mpoly_mpoly_dist import mpoly_mpoly_cfs
from distance.mpoly_mpoly_dist import mpoly_mpoly_verify_cfs

THETA_LIST = [i * pi / 10 for i in range(1, 11)]
N_LIST = [n for n in range(5, 500, 45)]


def non_optimized():
    point_poly_kt_list = []
    for theta in THETA_LIST:
        for n in N_LIST:
            while True:
                poly = Poly.random(0, 0, n, n, n, center=False, digits=12)
                mr = MPoly(poly, Pose(), Pose(Point(0, n), theta))
                p1 = Point.random(2 * n, n, 3 * n, 2 * n)
                mp = MPoint(p1, Point(p1.x, p1.y - n))
                t = 0
                cfs, _, _ = mpoint_mpoly_cfs(mp, mr)
                if not mpoint_mpoly_verify_cfs(mp, mr, cfs):
                    continue
                for i in range(10):
                    _, d_0, d = mpoint_mpoly_cfs(mp, mr)
                    t += d
                point_poly_kt_list.append((len(cfs), t / 10))
                break
    point_poly_kt_list.sort(key=lambda tup: tup[0])
    poly_poly_kt_list = []
    for theta in THETA_LIST:
        for n in N_LIST:
            while True:
                poly1 = Poly.random(0, 0, n, n, n, center=False, digits=12)
                mr1 = MPoly(poly1, Pose(), Pose(Point(0, n), theta))
                poly2 = Poly.random(
                    2 * n, n, 3 * n, 2 * n, n, center=False, digits=12)
                mr2 = MPoly(poly2, Pose(), Pose(Point(0, -n), theta))
                t = 0
                cfs, _, _ = mpoly_mpoly_cfs(mr1, mr2)
                if not mpoly_mpoly_verify_cfs(mr1, mr2, cfs):
                    continue
                for i in range(10):
                    _, d_0, d = mpoly_mpoly_cfs(mr1, mr2)
                    t += d
                poly_poly_kt_list.append((len(cfs), t / 10))
                break
    poly_poly_kt_list.sort(key=lambda tup: tup[0])
    return point_poly_kt_list, poly_poly_kt_list


def optimized():
    point_poly_kt_list = []
    for n in N_LIST:
        while True:
            poly = Poly.random(0, 0, n, n, n, center=False, digits=12)
            mr = MPoly(poly, Pose(), Pose(Point(0, n), 0))
            p1 = Point.random(2 * n, n, 3 * n, 2 * n)
            mp = MPoint(p1, Point(p1.x, p1.y - n))
            t = 0
            cfs, _, _ = mpoint_mpoly_cfs(mp, mr)
            if not mpoint_mpoly_verify_cfs(mp, mr, cfs):
                continue
            for i in range(10):
                _, d_0, d = mpoint_mpoly_cfs(mp, mr)
                t += d
            point_poly_kt_list.append((len(cfs), t / 10))
            break
    point_poly_kt_list.sort(key=lambda tup: tup[0])
    poly_poly_kt_list = []
    for n in N_LIST:
        while True:
            poly1 = Poly.random(0, 0, n, n, n, center=False, digits=12)
            mr1 = MPoly(poly1, Pose(), Pose(Point(0, n), 0))
            poly2 = Poly.random(
                2 * n, n, 3 * n, 2 * n, n, center=False, digits=12)
            mr2 = MPoly(poly2, Pose(), Pose(Point(0, -n), 0))
            t = 0
            cfs, _, _ = mpoly_mpoly_cfs(mr1, mr2)
            if not mpoly_mpoly_verify_cfs(mr1, mr2, cfs):
                continue
            for i in range(10):
                _, d_0, d = mpoly_mpoly_cfs(mr1, mr2)
                t += d
            poly_poly_kt_list.append((len(cfs), t / 10))
            break
    poly_poly_kt_list.sort(key=lambda tup: tup[0])
    return point_poly_kt_list, poly_poly_kt_list


def run(save_path=None, show_plot=True):
    pt_pl_kts, pl_pl_kts = non_optimized()
    avg_pt_pl = sum([t / k for k, t in pt_pl_kts]) / len(pt_pl_kts)
    avg_pl_pl = sum([t / k for k, t in pl_pl_kts]) / len(pl_pl_kts)
    print(f"Avg pt-pl non_optimized (ms): {avg_pt_pl}")
    print(f"Avg pl-pl non_optimized (ms): {avg_pl_pl}")

    labels = [r"point-poly (ms)", r"poly-poly (ms)"]
    plt.plot(
        [kt[0] for kt in pt_pl_kts], [kt[1] for kt in pt_pl_kts],
        label=labels[0])
    plt.plot(
        [kt[0] for kt in pl_pl_kts], [kt[1] for kt in pl_pl_kts],
        label=labels[1])
    plt.xlabel('k')
    plt.ylabel('t', rotation=0, labelpad=15)
    plt.legend()
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path + 'distance_t_k_non_opt.png', bbox_inches='tight')
    if show_plot:
        plt.show()
    else:
        plt.clf()

    pt_pl_kts, pl_pl_kts = optimized()
    avg_pt_pl = sum([t / k for k, t in pt_pl_kts]) / len(pt_pl_kts)
    avg_pl_pl = sum([t / k for k, t in pl_pl_kts]) / len(pl_pl_kts)
    print(f"Avg pt-pl optimized (ms): {avg_pt_pl}")
    print(f"Avg pl-pl optimized (ms): {avg_pl_pl}")

    labels = [r"point-poly (ms)", r"poly-poly (ms)"]
    plt.plot(
        [kt[0] for kt in pt_pl_kts], [kt[1] for kt in pt_pl_kts],
        label=labels[0])
    plt.plot(
        [kt[0] for kt in pl_pl_kts], [kt[1] for kt in pl_pl_kts],
        label=labels[1])
    plt.xlabel('k')
    plt.ylabel('t', rotation=0, labelpad=15)
    plt.legend()
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path + 'distance_t_k_opt.png', bbox_inches='tight')
    if show_plot:
        plt.show()
    else:
        plt.clf()