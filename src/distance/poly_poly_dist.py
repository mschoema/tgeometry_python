from math import sqrt

from distance.v_clip import compute_dist2
from distance.v_clip import poly_poly_v_clip


def poly_poly_dist(r1, r2):
  cf1, cf2 = poly_poly_v_clip(r1, r2)
  v1 = cf1 // 2
  v2 = cf2 // 2
  if cf1 % 2 == 0 and cf2 % 2 == 0:
    p = r1.get_v(v1)
    q = r2.get_v(v1)
    dist = q.dist(p)
  elif cf1 % 2 == 0 and cf2 % 2 == 1:
    p = r1.get_v(v1)
    q = r2.get_v(v2)
    r = r2.get_v(v2 + 1)
    dist = sqrt(compute_dist2(p, q, r))
  elif cf1 % 2 == 1 and cf2 % 2 == 0:
    p = r2.get_v(v2)
    q = r1.get_v(v1)
    r = r1.get_v(v1 + 1)
    dist = sqrt(compute_dist2(p, q, r))
    pass
  else:  # cf1 % 2 == 1 and cf2 % 2 == 1
    print("This should never happen")
    dist = -1
  return dist