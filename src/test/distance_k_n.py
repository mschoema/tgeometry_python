import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
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
from distance.mpoint_mpoly_dist import mpoint_mpoly_verify_cfs_fast
from distance.mpoly_mpoly_dist import mpoly_mpoly_cfs
from distance.mpoly_mpoly_dist import mpoly_mpoly_verify_cfs_fast

THETA_LIST = [pi, pi / 2, 0]
N_LIST = [n for n in range(5, 500, 5)]


def point_to_poly():
  k_lists = []
  for i, theta in enumerate(THETA_LIST):
    k_lists.append([])
    for n in N_LIST:
      c = 0
      k = 0
      while c < 5:
        poly = Poly.random(0, 0, n, n, n, center=False, digits=12)
        mr = MPoly(poly, Pose(), Pose(Point(0, n), theta))
        p1 = Point.random(2 * n, n, 3 * n, 2 * n)
        mp = MPoint(p1, Point(p1.x, p1.y - n))
        cfs, _, _ = mpoint_mpoly_cfs(mp, mr)
        if mpoint_mpoly_verify_cfs_fast(mp, mr, cfs):
          k += len(cfs)
          c += 1
      k_lists[i].append(k / c)
  return k_lists


def poly_to_poly():
  k_lists = []
  for i, theta in enumerate(THETA_LIST):
    k_lists.append([])
    for n in N_LIST:
      c = 0
      k = 0
      while c < 5:
        poly1 = Poly.random(0, 0, n, n, n, center=False, digits=12)
        mr1 = MPoly(poly1, Pose(), Pose(Point(0, n), theta))
        poly2 = Poly.random(2 * n, n, 3 * n, 2 * n, n, center=False, digits=12)
        mr2 = MPoly(poly2, Pose(), Pose(Point(0, -n), theta))
        cfs, _, _ = mpoly_mpoly_cfs(mr1, mr2)
        if mpoly_mpoly_verify_cfs_fast(mr1, mr2, cfs):
          k += len(cfs)
          c += 1
      k_lists[i].append(k / c)
  return k_lists


def run(save_path=None, show_plot=True):
  k_lists_point_poly = point_to_poly()

  labels = [
      r"k(n,$\theta$=$\pi$)", r"k(n,$\theta$=$\pi/2$)", r"k(n,$\theta$=0)"
  ]
  for i, k_list in enumerate(k_lists_point_poly):
    plt.plot(N_LIST, k_list, label=labels[i])
  plt.xlabel('n')
  plt.ylabel('k', rotation=0, labelpad=15)
  plt.legend()
  plt.tight_layout()
  if save_path is not None:
    plt.savefig(save_path + 'distance_k_n_point_poly.png', bbox_inches='tight')
  if show_plot:
    plt.show()
  else:
    plt.clf()

  k_lists_poly_poly = poly_to_poly()

  labels = [
      r"k(n+m,$\theta$=$\pi$)", r"k(n+m,$\theta$=$\pi/2$)",
      r"k(n+m,$\theta$=0)"
  ]
  for i, k_list in enumerate(k_lists_poly_poly):
    plt.plot([2 * n for n in N_LIST], k_list, label=labels[i])
  plt.xlabel('n + m')
  plt.ylabel('k', rotation=0, labelpad=15)
  plt.legend()
  plt.tight_layout()
  if save_path is not None:
    plt.savefig(save_path + 'distance_k_n_poly_poly.png', bbox_inches='tight')
  if show_plot:
    plt.show()
  else:
    plt.clf()