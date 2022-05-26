from time import process_time

from classes.constants import EPSILON
from distance.v_clip import point_poly_v_clip

EQ_0 = True
EQ_1 = False


def f_mpoint_mpoly(mp, mr, v, t, eq_0):
    p = mp.at(t)
    q = mr.at_v(v, t)
    r = mr.at_v(v + 1, t)
    if eq_0:
        return (p.x - q.x) * (r.x - q.x) + (p.y - q.y) * (r.y - q.y)
    else:
        return (p.x - r.x) * (r.x - q.x) + (p.y - r.y) * (r.y - q.y)


def solve_s_mpoint_mpoly(mp, mr, v, t_prev, eq_0):

    # No rotation during the movement
    if mr.linear:
        p = mp.at(0)
        q = mr.at_v(v, 0)
        r = mr.at_v(v + 1, 0)
        dxp = mp.end_point.x - mp.start_point.x
        dyp = mp.end_point.y - mp.start_point.y
        dxr = mr.end_pose.pos.x - mr.start_pose.pos.x
        dyr = mr.end_pose.pos.y - mr.start_pose.pos.y
        discr = (dxp - dxr) * (r.x - q.x) + (dyp - dyr) * (r.y - q.y)
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
    vl = f_mpoint_mpoly(mp, mr, v, ts, eq_0)
    v0 = f_mpoint_mpoly(mp, mr, v, (ts + te) / 2, eq_0)
    vr = f_mpoint_mpoly(mp, mr, v, te, eq_0)
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
        v0 = f_mpoint_mpoly(mp, mr, v, t0, eq_0)
        if abs(v0) < EPSILON:
            break
        if vl * v0 <= 0:
            tr = t0
            vr = v0
        else:
            tl = t0
            vl = v0
    return t0


def mpoint_mpoly_cfs_bf(mp, mr, n=1000):
    cf = point_poly_v_clip(mp.at(0), mr.at(0))
    cfs = [(0, cf)]
    for i in range(1, n + 1):
        t = i / n
        _, cf_prev = cfs[-1]
        cf = point_poly_v_clip(mp.at(t), mr.at(t), cf_prev)
        if cf != cf_prev:
            cfs.append((t, cf))
    return cfs


def mpoint_mpoly_verify_cfs(mp, mr, cfs):
    correct = True
    for i in range(len(cfs)):
        ts, cf = cfs[i]
        te = 1 if i == len(cfs) - 1 else cfs[i + 1][0]
        cf_bf = point_poly_v_clip(
            mp.at((ts + te) / 2), mr.at((ts + te) / 2), cf)
        if cf != cf_bf:
            correct = False
            break
    return correct


def mpoint_mpoly_verify_cfs_fast(mp, mr, cfs):
    t, cf = cfs[-1]
    cf_bf = point_poly_v_clip(mp.at(t), mr.at(t), cf)
    return cf != cf_bf


# Assumptions:
# - No special cases at t=0
# - No intersection during the movement
#   (intersection at edge is not tested)
def mpoint_mpoly_cfs(mp, mr):
    start_0 = process_time()
    cf = point_poly_v_clip(mp.at(0), mr.at(0))
    end_0 = start = process_time()
    cfs = [(0, cf)]
    t_prev = 0
    cf_prev = -1
    i = 0
    max_i = 4 * mr.poly.n
    while True:
        i += 1
        v = cf // 2
        t_1 = t_2 = 2
        if cf % 2 == 0:
            if cf_prev == -1 or (cf_prev + 1) % (2 * mr.poly.n) == cf:
                t_1 = solve_s_mpoint_mpoly(mp, mr, v, t_prev, EQ_0)
            if cf_prev == -1 or (cf_prev - 1) % (2 * mr.poly.n) == cf:
                v_prev = (v - 1) % mr.poly.n
                t_2 = solve_s_mpoint_mpoly(mp, mr, v_prev, t_prev, EQ_1)
        else:  # cf % 2 == 1
            if cf_prev == -1 or (cf_prev + 1) % (2 * mr.poly.n) == cf:
                t_1 = solve_s_mpoint_mpoly(mp, mr, v, t_prev, EQ_1)
            if cf_prev == -1 or (cf_prev - 1) % (2 * mr.poly.n) == cf:
                t_2 = solve_s_mpoint_mpoly(mp, mr, v, t_prev, EQ_0)

        cf_prev = cf

        # No change in cf
        if t_1 == 2 and t_2 == 2:
            break
        # Intersection through vertex
        elif cf % 2 == 0 and abs(t_1 - t_2) < EPSILON:
            print("Cfs: Intersection detected. #2")
            break
        # Precision error, skip the edge
        elif cf % 2 == 1 and abs(t_1 - t_2) < EPSILON:
            # print("Cfs: Skipping edge")
            if len(cfs) < 2:
                print("Cfs: Problem, can't skip edge")
                break
            if cfs[-2][1] == (cf - 1) % (2 * mr.poly.n):
                cf = (cf + 1) % (2 * mr.poly.n)
            else:  # cfs[-2][1] == (cf + 1) % (2 * mr.poly.n)
                cf = (cf - 1) % (2 * mr.poly.n)
            cfs.append((t_1, cf))
            t_prev = t_1
        # Go to next cf
        elif t_1 < t_2:
            cf = (cf + 1) % (2 * mr.poly.n)
            cfs.append((t_1, cf))
            t_prev = t_1
        # Go to previous cf
        else:  # t_2 < t_1
            cf = (cf - 1) % (2 * mr.poly.n)
            cfs.append((t_2, cf))
            t_prev = t_2

        if i > max_i:
            print("Cfs: Max iterations reached")
            break
    end = process_time()
    return cfs, (end_0 - start_0) * 1000, (end - start) * 1000