from time import process_time

from classes.constants import EPSILON
from classes.constants import LOOP_MAX_ITERS
from distance.v_clip import compute_s
from distance.v_clip import mpoly_mpoly_v_clip

EQ_0 = True
EQ_1 = False


def f_mpoly_mpoly(mr1, mr2, v1, v2, t, eq_0):
  p = mr1.at_v(v1, t)
  q = mr2.at_v(v2, t)
  r = mr2.at_v(v2 + 1, t)
  if eq_0:
    return (p.x - q.x) * (r.x - q.x) + (p.y - q.y) * (r.y - q.y)
  else:
    return (p.x - r.x) * (r.x - q.x) + (p.y - r.y) * (r.y - q.y)


def solve_s_mpoly_mpoly(mr1, mr2, v1, v2, t_prev, eq_0):

  # No rotation during the movement
  if mr1.linear and mr2.linear:
    p = mr1.at_v(v1, 0)
    q = mr2.at_v(v2, 0)
    r = mr2.at_v(v2 + 1, 0)
    dx1 = mr1.end_pose.pos.x - mr1.start_pose.pos.x
    dy1 = mr1.end_pose.pos.y - mr1.start_pose.pos.y
    dx2 = mr2.end_pose.pos.x - mr2.start_pose.pos.x
    dy2 = mr2.end_pose.pos.y - mr2.start_pose.pos.y
    discr = (dx1 - dx2) * (r.x - q.x) + (dy1 - dy2) * (r.y - q.y)
    if abs(discr) < EPSILON:
      return 2
    if eq_0:
      t = ((q.x - p.x) * (r.x - q.x) + (q.y - p.y) * (r.y - q.y)) / discr
    else:
      t = ((r.x - p.x) * (r.x - q.x) + (r.y - p.y) * (r.y - q.y)) / discr
    if (t > t_prev + EPSILON and t < 1 - EPSILON):
      return t
    return 2

  ts = t_prev
  te = 1
  vl = f_mpoly_mpoly(mr1, mr2, v1, v2, ts, eq_0)
  v0 = f_mpoly_mpoly(mr1, mr2, v1, v2, (ts + te) / 2, eq_0)
  vr = f_mpoly_mpoly(mr1, mr2, v1, v2, te, eq_0)
  if abs(vl) > EPSILON and vl * v0 < 0:
    tl = ts
    tr = (ts + te) / 2
    vr = v0
  elif v0 * vr < 0:
    tl = (ts + te) / 2
    tr = te
    vl = v0
  else:
    return 2

  i = 0
  while abs(tr - tl) >= 2 * EPSILON and i < 100:
    i += 1
    t0 = (tl * vr - tr * vl) / (vr - vl)
    v0 = f_mpoly_mpoly(mr1, mr2, v1, v2, t0, eq_0)
    if abs(v0) < EPSILON:
      break
    if vl * v0 <= 0:
      tr = t0
      vr = v0
    else:
      tl = t0
      vl = v0
  return t0


def f_parallel_mpoly_mpoly(mr1, mr2, v1, v2, t):
  ps = mr1.at_v(v1, t)
  pe = mr1.at_v(v1 + 1, t)
  qs = mr2.at_v(v2, t)
  qe = mr2.at_v(v2 + 1, t)
  return (pe.x - ps.x) * (qe.y - qs.y) - (pe.y - ps.y) * (qe.x - qs.x)


def solve_parallel_edges_mpoly_mpoly(mr1, mr2, v1, v2, t_prev):

  # No rotation during the movement
  # Edges do not rotate, so no need to solve this
  if mr1.linear and mr2.linear:
    return 2

  ts = t_prev
  te = 1
  vl = f_parallel_mpoly_mpoly(mr1, mr2, v1, v2, ts)
  v0 = f_parallel_mpoly_mpoly(mr1, mr2, v1, v2, (ts + te) / 2)
  vr = f_parallel_mpoly_mpoly(mr1, mr2, v1, v2, te)
  if abs(vl) > EPSILON and vl * v0 < 0:
    tl = ts
    tr = (ts + te) / 2
    vr = v0
  elif v0 * vr < 0:
    tl = (ts + te) / 2
    tr = te
    vl = v0
  else:
    return 2

  i = 0
  while abs(tr - tl) >= EPSILON and i < 100:
    i += 1
    t0 = (tl * vr - tr * vl) / (vr - vl)
    v0 = f_parallel_mpoly_mpoly(mr1, mr2, v1, v2, t0)
    if abs(v0) < EPSILON:
      break
    if vl * v0 <= 0:
      tr = t0
      vr = v0
    else:
      tl = t0
      vl = v0
  return t0


def mpoly_mpoly_cfs_bf(mr1, mr2, n=1000):
  cf1, cf2 = mpoly_mpoly_v_clip(mr1, mr2, 0)
  cfs = [(0, cf1, cf2)]
  for i in range(1, n + 1):
    t = i / n
    _, cf1_prev, cf2_prev = cfs[-1]
    cf1, cf2 = mpoly_mpoly_v_clip(mr1, mr2, t, cf1_prev, cf2_prev)
    if cf1 != cf1_prev or cf2 != cf2_prev:
      cfs.append((t, cf1, cf2))
  return cfs


def mpoly_mpoly_verify_cfs(mr1, mr2, cfs):
  correct = True
  for i in range(len(cfs)):
    ts, cf1, cf2 = cfs[i]
    te = 1 if i == len(cfs) - 1 else cfs[i + 1][0]
    cf1_bf, cf2_bf = mpoly_mpoly_v_clip(mr1, mr2, (ts + te) / 2, cf1, cf2)
    if cf1 != cf1_bf or cf2 != cf2_bf:
      correct = False
      break
  return correct


def mpoly_mpoly_verify_cfs_fast(mr1, mr2, cfs):
  t, cf1, cf2 = cfs[-1]
  cf1_bf, cf2_bf = mpoly_mpoly_v_clip(mr1, mr2, (1 + t) / 2, cf1, cf2)
  return cf1 == cf1_bf and cf2 == cf2_bf


# Assumptions:
# - No special cases at t=0
# - No intersection during the movement
def mpoly_mpoly_cfs(mr1, mr2):
  start_0 = process_time()
  cf1, cf2 = mpoly_mpoly_v_clip(mr1, mr2, 0)
  end_0 = start = process_time()
  cfs = [(0, cf1, cf2)]
  cf1_prev, cf2_prev = -1, -1
  i = 0
  max_i = 4 * (mr1.poly.n + mr2.poly.n)
  while True:
    i += 1
    v1 = cf1 // 2
    v2 = cf2 // 2
    t_1 = t_2 = t_3 = t_4 = 2
    t_prev, _, _ = cfs[-1]

    if cf1 % 2 == 1 and cf2 % 2 == 1:
      print("This should never happen")
      break

    if cf1 % 2 == 0 and cf2 % 2 == 0:
      # Compute solutions of equations
      if (cf2_prev == -1 or (cf2_prev + 1) % (2 * mr2.poly.n) == cf2
          or (cf2_prev + 2) % (2 * mr2.poly.n) == cf2):
        t_1 = solve_s_mpoly_mpoly(mr1, mr2, v1, v2, t_prev, EQ_0)
      if (cf2_prev == -1 or (cf2_prev - 1) % (2 * mr2.poly.n) == cf2
          or (cf2_prev - 2) % (2 * mr2.poly.n) == cf2):
        v2_prev = (v2 - 1) % mr2.poly.n
        t_2 = solve_s_mpoly_mpoly(mr1, mr2, v1, v2_prev, t_prev, EQ_1)
      if (cf1_prev == -1 or (cf1_prev + 1) % (2 * mr1.poly.n) == cf1
          or (cf1_prev + 2) % (2 * mr1.poly.n) == cf1):
        t_3 = solve_s_mpoly_mpoly(mr2, mr1, v2, v1, t_prev, EQ_0)
      if (cf1_prev == -1 or (cf1_prev - 1) % (2 * mr1.poly.n) == cf1
          or (cf1_prev - 2) % (2 * mr1.poly.n) == cf1):
        v1_prev = (v1 - 1) % mr1.poly.n
        t_4 = solve_s_mpoly_mpoly(mr2, mr1, v2, v1_prev, t_prev, EQ_1)
      # Change cfs as needed
      if t_1 <= 1 and t_1 < t_2 and t_1 < t_3 and t_1 < t_4:
        cf2_prev = cf2
        cf2 = (cf2 + 1) % (2 * mr2.poly.n)
        cfs.append((t_1, cf1, cf2))
      elif t_2 <= 1 and t_2 < t_1 and t_2 < t_3 and t_2 < t_4:
        cf2_prev = cf2
        cf2 = (cf2 - 1) % (2 * mr2.poly.n)
        cfs.append((t_2, cf1, cf2))
      elif t_3 <= 1 and t_3 < t_1 and t_3 < t_2 and t_3 < t_4:
        cf1_prev = cf1
        cf1 = (cf1 + 1) % (2 * mr1.poly.n)
        cfs.append((t_3, cf1, cf2))
      elif t_4 <= 1 and t_4 < t_1 and t_4 < t_2 and t_4 < t_3:
        cf1_prev = cf1
        cf1 = (cf1 - 1) % (2 * mr1.poly.n)
        cfs.append((t_4, cf1, cf2))
      else:
        break
    else:
      if cf1 % 2 == 0 and cf2 % 2 == 1:
        tmp_mr1, tmp_mr2 = mr1, mr2
        tmp_cf1, tmp_cf2 = cf1, cf2
        tmp_cf1_prev, tmp_cf2_prev = cf1_prev, cf2_prev
      else:  # cf1 % 2 == 1 and cf2 % 2 == 0
        tmp_mr1, tmp_mr2 = mr2, mr1
        tmp_cf1, tmp_cf2 = cf2, cf1
        tmp_cf1_prev, tmp_cf2_prev = cf2_prev, cf1_prev
      tmp_v1 = tmp_cf1 // 2
      tmp_v2 = tmp_cf2 // 2

      # Compute solutions of equations
      if (tmp_cf2_prev == -1
          or (tmp_cf2_prev + 1) % (2 * tmp_mr2.poly.n) == tmp_cf2
          or (tmp_cf2_prev + 2) % (2 * tmp_mr2.poly.n) == tmp_cf2):
        t_1 = solve_s_mpoly_mpoly(tmp_mr1, tmp_mr2, tmp_v1, tmp_v2, t_prev,
                                  EQ_1)
      if (tmp_cf2_prev == -1
          or (tmp_cf2_prev - 1) % (2 * tmp_mr2.poly.n) == tmp_cf2
          or (tmp_cf2_prev - 2) % (2 * tmp_mr2.poly.n) == tmp_cf2):
        t_2 = solve_s_mpoly_mpoly(tmp_mr1, tmp_mr2, tmp_v1, tmp_v2, t_prev,
                                  EQ_0)
      t_3 = solve_parallel_edges_mpoly_mpoly(tmp_mr1, tmp_mr2, tmp_v1, tmp_v2,
                                             t_prev)
      t_4 = solve_parallel_edges_mpoly_mpoly(
          tmp_mr1, tmp_mr2, (tmp_v1 - 1) % tmp_mr1.poly.n, tmp_v2, t_prev)

      # Change cfs as needed
      if t_1 <= 1 and t_1 < t_2 and t_1 < t_3 and t_1 < t_4:
        tmp_cf2_prev = tmp_cf2
        tmp_cf2 = (tmp_cf2 + 1) % (2 * tmp_mr2.poly.n)
        t_next = t_1
      elif t_2 <= 1 and t_2 < t_1 and t_2 < t_3 and t_2 < t_4:
        tmp_cf2_prev = tmp_cf2
        tmp_cf2 = (tmp_cf2 - 1) % (2 * tmp_mr2.poly.n)
        t_next = t_2
      elif t_3 <= 1 and t_3 < t_1 and t_3 < t_2 and t_3 < t_4:
        s12 = compute_s(
            tmp_mr1.at_v(tmp_v1 + 1, t_3), tmp_mr2.at_v(tmp_v2, t_3),
            tmp_mr2.at_v(tmp_v2 + 1, t_3))
        s21 = compute_s(
            tmp_mr2.at_v(tmp_v2, t_3), tmp_mr1.at_v(tmp_v1, t_3),
            tmp_mr1.at_v(tmp_v1 + 1, t_3))
        if 0 < s12 < 1:
          tmp_cf1_prev = tmp_cf1
          tmp_cf1 = (tmp_cf1 + 2) % (2 * tmp_mr1.poly.n)
        elif 0 < s21 < 1:
          tmp_cf1_prev, tmp_cf2_prev = tmp_cf1, tmp_cf2
          tmp_cf1 = (tmp_cf1 + 1) % (2 * tmp_mr1.poly.n)
          tmp_cf2 = (tmp_cf2 - 1) % (2 * tmp_mr2.poly.n)
        else:
          tmp_cf1_prev, tmp_cf2_prev = tmp_cf1, tmp_cf2
          tmp_cf1 = (tmp_cf1 + 2) % (2 * tmp_mr1.poly.n)
          tmp_cf2 = (tmp_cf2 - 1) % (2 * tmp_mr2.poly.n)
        t_next = t_3
      elif t_4 <= 1 and t_4 < t_1 and t_4 < t_2 and t_4 < t_3:
        s12 = compute_s(
            tmp_mr1.at_v(tmp_v1 - 1, t_4), tmp_mr2.at_v(tmp_v2, t_4),
            tmp_mr2.at_v(tmp_v2 + 1, t_4))
        s21 = compute_s(
            tmp_mr2.at_v(tmp_v2 + 1, t_4), tmp_mr1.at_v(tmp_v1 - 1, t_4),
            tmp_mr1.at_v(tmp_v1, t_4))
        if 0 < s12 < 1:
          tmp_cf1_prev = tmp_cf1
          tmp_cf1 = (tmp_cf1 - 2) % (2 * tmp_mr1.poly.n)
        elif 0 < s21 < 1:
          tmp_cf1_prev, tmp_cf2_prev = tmp_cf1, tmp_cf2
          tmp_cf1 = (tmp_cf1 - 1) % (2 * tmp_mr1.poly.n)
          tmp_cf2 = (tmp_cf2 + 1) % (2 * tmp_mr2.poly.n)
        else:
          tmp_cf1_prev, tmp_cf2_prev = tmp_cf1, tmp_cf2
          tmp_cf1 = (tmp_cf1 - 2) % (2 * tmp_mr1.poly.n)
          tmp_cf2 = (tmp_cf2 + 1) % (2 * tmp_mr2.poly.n)
        t_next = t_4
      else:
        break

      if cf1 % 2 == 0 and cf2 % 2 == 1:
        mr1, mr2 = tmp_mr1, tmp_mr2
        cf1, cf2 = tmp_cf1, tmp_cf2
        cf1_prev, cf2_prev = tmp_cf1_prev, tmp_cf2_prev
      else:  # cf1 % 2 == 1 and cf2 % 2 == 0
        mr2, mr1 = tmp_mr1, tmp_mr2
        cf2, cf1 = tmp_cf1, tmp_cf2
        cf2_prev, cf1_prev = tmp_cf1_prev, tmp_cf2_prev
      cfs.append((t_next, cf1, cf2))

    if i > max_i:
      print("Cfs: Max iterations reached")
      break
  end = process_time()
  return cfs, (end_0 - start_0) * 1000, (end - start) * 1000


# Find a change in closest features between the two given features
# Might return a redundant instant to be removed later
def mpoly_mpoly_cfs_change(mr1, mr2, cft_start, cft_end):
  t_start, cf1_start, cf2_start = cft_start
  t_end, cf1_end, cf2_end = cft_end
  t1, t2 = t_start, t_end
  cf1_middle, cf2_middle = cf1_start, cf2_start
  while (t2 - t1 > EPSILON):
    t_middle = (t2 + t1) / 2
    cf1_middle, cf2_middle = mpoly_mpoly_v_clip(mr1, mr2, t_middle, cf1_middle,
                                                cf2_middle)
    # If start == middle and end == middle, assume nothing changes
    if (cf1_start == cf1_middle and cf2_start == cf2_middle
        and cf1_end == cf1_middle and cf2_end == cf2_middle):
      return False, None, None, None
    # If start == middle, update t1
    elif (cf1_start == cf1_middle and cf2_start == cf2_middle):
      t1 = t_middle
    # If end == middle, update t2
    elif (cf1_end == cf1_middle and cf2_end == cf2_middle):
      t2 = t_middle
    # If start != middle and end != middle, return middle
    else:
      # The value t_middle might not correspond to the exact
      # transition time, but that will be determined in a future call
      # This tuple will thus most probably be redundant at the end
      return True, t_middle, cf1_middle, cf2_middle
  # If t1 == t_start or t2 == t_end, no change happens
  if (t1 == t_start or t2 == t_end):
    return False, None, None, None
  elif (cf1_end == cf1_middle and cf2_end == cf2_middle):
    return True, t_middle, cf1_end, cf2_end
  else:
    return True, t2, cf1_end, cf2_end


# Iteratively finds the changes in closest features
# In a binary search / bracketing approach
def mpoly_mpoly_cfs_iter(mr1, mr2):
  start_0 = process_time()
  cf1_start, cf2_start = mpoly_mpoly_v_clip(mr1, mr2, 0)
  cf1_end, cf2_end = mpoly_mpoly_v_clip(mr1, mr2, 1, cf1_start, cf2_start)
  end_0 = start = process_time()
  # Initial and final closest features
  cfs = [(0, cf1_start, cf2_start), (1, cf1_end, cf2_end)]
  i = 1
  loop = 0
  while i < len(cfs):
    # Detect a change of closest features between the given instants
    found, t, cf1, cf2 = mpoly_mpoly_cfs_change(mr1, mr2, cfs[i - 1], cfs[i])
    if (found):
      cfs.insert(i, (t, cf1, cf2))
    else:
      i += 1

    loop += 1
    if (loop == LOOP_MAX_ITERS):
      print(f"Mpoly-Mpoly Cfs: Cycle detected")
      break
  end = process_time()
  # Remove redundant instants
  for i in range(len(cfs) - 1, 0, -1):
    t, cf1, cf2 = cfs[i]
    t_prev, cf1_prev, cf2_prev = cfs[i - 1]
    if (cf1 == cf1_prev and cf2 == cf2_prev):
      cfs.pop(i)
  return cfs, (end_0 - start_0) * 1000, (end - start) * 1000