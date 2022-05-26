from math import sqrt

from distance.v_clip import compute_dist2
from distance.v_clip import point_poly_v_clip


def point_poly_dist(p, r):
    cf = point_poly_v_clip(p, r)
    v = cf // 2
    if cf % 2 == 0:
        q = r.get_v(v)
        dist = q.dist(p)
    else:  # cf % 2 == 1
        q = r.get_v(v)
        r = r.get_v(v + 1)
        dist = sqrt(compute_dist2(p, q, r))
    return dist