from math import inf

from classes.constants import EPSILON
from classes.point import Point


def compute_s(p, q, r):
    return ((p.x - q.x) * (r.x - q.x) + (p.y - q.y) *
            (r.y - q.y)) / (pow(r.x - q.x, 2) + pow(r.y - q.y, 2))


def compute_angle(p, q, r, is_ccw):
    angle = (p.x - q.x) * (r.y - q.y) - (p.y - q.y) * (r.x - q.x)
    return angle if is_ccw else -angle


def compute_dist2(p, q, r):
    s = compute_s(p, q, r)
    x = q.x + (r.x - q.x) * s
    y = q.y + (r.y - q.y) * s
    return pow(p.x - x, 2) + pow(p.y - y, 2)


def compute_point_edge_dist2(p, q, r):
    s = compute_s(p, q, r)
    if s <= 0:
        return pow(p.x - q.x, 2) + pow(p.y - q.y, 2)
    elif s >= 1:
        return pow(p.x - r.x, 2) + pow(p.y - r.y, 2)
    else:
        x = q.x + (r.x - q.x) * s
        y = q.y + (r.y - q.y) * s
        return pow(p.x - x, 2) + pow(p.y - y, 2)


def point_poly_cf_bf(p, r):
    d_min = inf
    cf = 0
    for v_r in range(r.n):
        r_curr = r.get_v(v_r)
        r_next = r.get_v(v_r + 1)
        if p.dist2(r_curr) <= d_min:
            d_min = p.dist2(r_curr)
            cf = v_r * 2
        if compute_point_edge_dist2(p, r_curr, r_next) <= d_min:
            d_min = compute_point_edge_dist2(p, r_curr, r_next)
            cf = v_r * 2 + 1
    return cf


def point_poly_v_clip(p, r, cf=0):
    i = 0
    cfs = set()
    while True:
        if cf in cfs:
            print("V-clip: looping")
            return point_poly_cf_bf(p, r)
        cfs.add(cf)
        i += 1
        v_r = cf // 2
        r_curr = r.get_v(v_r)
        r_next = r.get_v(v_r + 1)
        s_rcurr_rnext = compute_s(p, r_curr, r_next)
        if cf % 2 == 0:
            r_prev = r.get_v(v_r - 1)
            s_rprev_rcurr = compute_s(p, r_prev, r_curr)
            if (abs(s_rprev_rcurr - 1) < EPSILON
                    and abs(s_rcurr_rnext) < EPSILON):
                print("V-clip: Intersection detected. #1")
                break
            elif s_rprev_rcurr < 1:
                cf = (cf - 1) % (2 * r.n)
            elif s_rcurr_rnext > 0:
                cf = (cf + 1) % (2 * r.n)
            else:
                break
        else:  # cf % 2 == 1
            if s_rcurr_rnext <= 0:
                cf = (cf - 1) % (2 * r.n)
            elif s_rcurr_rnext >= 1:
                cf = (cf + 1) % (2 * r.n)
            elif compute_angle(p, r_curr, r_next, r.is_ccw) < EPSILON:
                # print("Escaping local min")
                dmax = -1
                for j in range(r.n):
                    v_r = j
                    r_curr = r.get_v(v_r)
                    r_next = r.get_v(v_r + 1)
                    if (compute_angle(p, r_curr, r_next, r.is_ccw) > EPSILON):
                        d = compute_dist2(p, r_curr, r_next)
                        if d > dmax:
                            dmax = d
                            cf = 2 * v_r + 1
                if dmax == -1:
                    print("V-clip: Intersection detected. #2")
                    break
            else:
                break
    return cf


def poly_poly_cf_bf(r1, r2):
    cf = None
    min_dist = inf
    for v_r1 in range(r1.n):
        for v_r2 in range(r2.n):
            r1_curr = r1.get_v(v_r1)
            r2_curr = r2.get_v(v_r2)
            r1_next = r1.get_v(v_r1 + 1)
            r2_next = r2.get_v(v_r2 + 1)
            d1 = r1_curr.dist2(r2_curr)
            d2 = compute_point_edge_dist2(r1_curr, r2_curr, r2_next)
            d3 = compute_point_edge_dist2(r2_curr, r1_curr, r1_next)
            if d1 == d2 and d1 == d3 and d1 <= min_dist:
                min_dist = d1
                cf = (2 * v_r1, 2 * v_r2)
            elif d2 < d1 and d2 < d3 and d2 < min_dist:
                min_dist = d2
                cf = (2 * v_r1, 2 * v_r2 + 1)
            elif d3 < d1 and d3 < d2 and d3 < min_dist:
                min_dist = d3
                cf = (2 * v_r1 + 1, 2 * v_r2)
    return cf


# TODO: check for intersections.
def poly_poly_v_clip(r1, r2, cf1=0, cf2=0):
    i = 0
    cfs = set()
    while True:
        if (cf1, cf2) in cfs:
            print("V-clip: Looping")
            return poly_poly_cf_bf(r1, r2)
        cfs.add((cf1, cf2))
        i += 1
        v1 = cf1 // 2
        v2 = cf2 // 2
        if cf1 % 2 == 0 and cf2 % 2 == 0:
            if compute_s(r1.get_v(v1), r2.get_v(v2 - 1), r2.get_v(v2)) < 1:
                cf2 = (cf2 - 1) % (2 * r2.n)
            elif compute_s(r1.get_v(v1), r2.get_v(v2), r2.get_v(v2 + 1)) > 0:
                cf2 = (cf2 + 1) % (2 * r2.n)
            elif compute_s(r2.get_v(v2), r1.get_v(v1 - 1), r1.get_v(v1)) < 1:
                cf1 = (cf1 - 1) % (2 * r1.n)
            elif compute_s(r2.get_v(v2), r1.get_v(v1), r1.get_v(v1 + 1)) > 0:
                cf1 = (cf1 + 1) % (2 * r1.n)
            else:
                break
        elif cf1 % 2 == 0 and cf2 % 2 == 1:
            r1_curr = r1.get_v(v1)
            s12_curr_next = compute_s(r1_curr, r2.get_v(v2), r2.get_v(v2 + 1))
            if s12_curr_next <= 0:
                cf2 = (cf2 - 1) % (2 * r2.n)
            elif s12_curr_next >= 1:
                cf2 = (cf2 + 1) % (2 * r2.n)
            else:
                r2_interp = Point.interpolate(
                    r2.get_v(v2), r2.get_v(v2 + 1), s12_curr_next)
                if compute_s(r2_interp, r1.get_v(v1 - 1), r1_curr) < 1:
                    cf1 = (cf1 - 1) % (2 * r1.n)
                elif compute_s(r2_interp, r1_curr, r1.get_v(v1 + 1)) > 0:
                    cf1 = (cf1 + 1) % (2 * r1.n)
                elif compute_angle(r1_curr, r2.get_v(v2), r2.get_v(v2 + 1),
                                   r2.is_ccw) < EPSILON:
                    # print("Escaping local min")
                    dmax = -1
                    for j in range(r2.n):
                        v2 = j
                        r2_curr = r2.get_v(v2)
                        r2_next = r2.get_v(v2 + 1)
                        if (compute_angle(r1_curr, r2_curr, r2_next, r2.is_ccw)
                                > EPSILON):
                            d = compute_dist2(r1_curr, r2_curr, r2_next)
                            if d > dmax:
                                dmax = d
                                cf2 = 2 * v2 + 1
                    if dmax == -1:
                        print("V-clip: Intersection detected. #1")
                        break
                else:
                    break
        elif cf1 % 2 == 1 and cf2 % 2 == 0:
            r2_curr = r2.get_v(v2)
            s21_curr_next = compute_s(r2_curr, r1.get_v(v1), r1.get_v(v1 + 1))
            if s21_curr_next <= 0:
                cf1 = (cf1 - 1) % (2 * r1.n)
            elif s21_curr_next >= 1:
                cf1 = (cf1 + 1) % (2 * r1.n)
            else:
                r1_interp = Point.interpolate(
                    r1.get_v(v1), r1.get_v(v1 + 1), s21_curr_next)
                if compute_s(r1_interp, r2.get_v(v2 - 1), r2_curr) < 1:
                    cf2 = (cf2 - 1) % (2 * r2.n)
                elif compute_s(r1_interp, r2_curr, r2.get_v(v2 + 1)) > 0:
                    cf2 = (cf2 + 1) % (2 * r2.n)
                elif compute_angle(r2_curr, r1.get_v(v1), r1.get_v(v1 + 1),
                                   r1.is_ccw) < EPSILON:
                    # print("Escaping local min")
                    dmax = -1
                    for j in range(r1.n):
                        v1 = j
                        r1_curr = r1.get_v(v1)
                        r1_next = r1.get_v(v1 + 1)
                        if (compute_angle(r2_curr, r1_curr, r1_next, r1.is_ccw)
                                > EPSILON):
                            d = compute_dist2(r2_curr, r1_curr, r1_next)
                            if d > dmax:
                                dmax = d
                                cf1 = 2 * v1 + 1
                    if dmax == -1:
                        print("V-clip: Intersection detected. #2")
                        break
                else:
                    break
        else:  # cf1 % 2 == 1 and cf2 % 2 == 1
            d1 = compute_point_edge_dist2(
                r1.get_v(v1), r2.get_v(v2), r2.get_v(v2 + 1))
            d2 = compute_point_edge_dist2(
                r1.get_v(v1 + 1), r2.get_v(v2), r2.get_v(v2 + 1))
            d3 = compute_point_edge_dist2(
                r2.get_v(v2), r1.get_v(v1), r1.get_v(v1 + 1))
            d4 = compute_point_edge_dist2(
                r2.get_v(v2 + 1), r1.get_v(v1), r1.get_v(v1 + 1))
            if d1 <= d2 and d1 <= d3 and d1 <= d4:
                cf1 = (cf1 - 1) % (2 * r1.n)
            elif d2 <= d3 and d2 <= d4:
                cf1 = (cf1 + 1) % (2 * r1.n)
            elif d3 <= d4:
                cf2 = (cf2 - 1) % (2 * r2.n)
            else:
                cf2 = (cf2 + 1) % (2 * r2.n)
        # if i > max_i:
        #     print("V-clip: Max iteration reached")
        #     break
    return cf1, cf2