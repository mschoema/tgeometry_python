from math import inf

from classes.constants import EPSILON
from classes.point import Point

VCLIP_MAX_ITERS = 10000

VCLIP_CONTINUE = 0
VCLIP_DISJOINT = 1
VCLIP_INTERSECT = 2


def compute_s(p, q, r):
  return ((p.x - q.x) * (r.x - q.x) + (p.y - q.y) *
          (r.y - q.y)) / (pow(r.x - q.x, 2) + pow(r.y - q.y, 2))


def compute_angle(p, q, r, is_ccw=True):
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


def point_poly_vertex_vertex(p, r, cf):
  v_r = cf // 2
  # Get enpoints of previous end next edge
  r_prev = r.get_v(v_r - 1)
  r_curr = r.get_v(v_r)
  # Check if the point is in r's Voronoi region
  s_prev = compute_s(p, r_prev, r_curr)
  if (s_prev < 1):
    # Go to the previous edge
    return VCLIP_CONTINUE, (cf - 1) % (2 * r.n)
  r_next = r.get_v(v_r + 1)
  s_next = compute_s(p, r_curr, r_next)
  if (s_next > 0):
    # Go the next edge
    return VCLIP_CONTINUE, (cf + 1) % (2 * r.n)

  # We found the closest feature
  if (s_prev > 1 or s_next < 0):
    return VCLIP_DISJOINT, cf
  # Point is on the vertex
  print("V-clip: Intersection detected")
  return VCLIP_INTERSECT, None


def point_poly_edge_vertex(p, r, cf):
  v_r = cf // 2
  # Get edge enpoints
  r_start = r.get_v(v_r)
  r_end = r.get_v(v_r + 1)
  # Check if the point is in the edge's Voronoin region
  s = compute_s(p, r_start, r_end)
  if (s < 0):
    # Go to the start vertex
    return VCLIP_CONTINUE, (cf - 1) % (2 * r.n)
  elif (s > 1):
    # Go to the end vertex
    return VCLIP_CONTINUE, (cf + 1) % (2 * r.n)

  # Check for local minimum
  if (compute_angle(p, r_start, r_end, r.is_ccw) < EPSILON):
    # Found local minimum
    dmax = -1
    r_start = r.get_v(0)
    for j in range(r.n):
      # Find the edge with the largest
      # positive distance to the given point
      r_end = r.get_v(j + 1)
      if (compute_angle(p, r_start, r_end, r.is_ccw) > EPSILON):
        dist = compute_dist2(p, r_start, r_end)
        if dist > dmax:
          dmax = dist
          cf = 2 * j + 1
      r_start = r_end

    # If no positive distance, point is inside polygon
    if dmax == -1:
      print("V-clip: Intersection detected")
      return VCLIP_INTERSECT, None
    return VCLIP_CONTINUE, cf
  # We found the closest feature
  return VCLIP_DISJOINT, cf


def point_poly_v_clip(p, r, cf=0):
  loop = 0
  state = VCLIP_CONTINUE
  while state == VCLIP_CONTINUE:
    # Compute next closest feature
    if (cf % 2 == 0):
      state, cf = point_poly_vertex_vertex(p, r, cf)
    else:  # (cf % 2 == 1)
      state, cf = point_poly_edge_vertex(p, r, cf)

    loop += 1
    if (loop == VCLIP_MAX_ITERS):
      print(f"V-clip: Cycle detected, current feature: {cf}")
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


def poly_poly_vertex_vertex(r1, r2, cf1, cf2):
  v_r1 = cf1 // 2
  v_r2 = cf2 // 2
  # Get enpoints of previous and next edges
  r1_prev = r1.get_v(v_r1 - 1)
  r1_curr = r1.get_v(v_r1)
  r2_curr = r2.get_v(v_r2)
  # Check if v2 is in v1's Voronoi region
  s1_prev = compute_s(r2_curr, r1_prev, r1_curr)
  if (s1_prev < 1):
    # Go to the previous edge of r1
    return VCLIP_CONTINUE, (cf1 - 1) % (2 * r1.n), cf2
  r1_next = r1.get_v(v_r1 + 1)
  s1_next = compute_s(r2_curr, r1_curr, r1_next)
  if (s1_next > 0):
    # Go to the previous edge of r1
    return VCLIP_CONTINUE, (cf1 + 1) % (2 * r1.n), cf2
  r2_prev = r2.get_v(v_r2 - 1)
  r2_next = r2.get_v(v_r2 + 1)
  # Check if v1 is in v2's Voronoi region
  s2_prev = compute_s(r1_curr, r2_prev, r2_curr)
  if (s2_prev < 1):
    # Go to the previous edge of r2
    return VCLIP_CONTINUE, cf1, (cf2 - 1) % (2 * r2.n)
  s2_next = compute_s(r1_curr, r2_curr, r2_next)
  if (s2_next > 0):
    # Go to the next edge of r2
    return VCLIP_CONTINUE, cf1, (cf2 + 1) % (2 * r2.n)

  # We found the closest feature
  if ((s1_prev > 1 or s1_next < 0) and (s2_prev > 1 or s2_next < 0)):
    return VCLIP_DISJOINT, cf1, cf2
  # Point is on the vertex
  print("V-clip: Intersection detected")
  return VCLIP_INTERSECT, None, None


def poly_poly_edge_vertex(r1, r2, cf1, cf2):
  v_r1 = cf1 // 2
  v_r2 = cf2 // 2
  # Get enpoints of current, previous and next edges
  r1_start = r1.get_v(v_r1)
  r1_end = r1.get_v(v_r1 + 1)
  r2_curr = r2.get_v(v_r2)

  # Check if v2 is in the Voronoi region of the edge of poly1
  s1 = compute_s(r2_curr, r1_start, r1_end)
  if (s1 < 0):
    # Go to the start vertex
    return VCLIP_CONTINUE, (cf1 - 1) % (2 * r1.n), cf2
  elif (s1 > 1):
    # Go to the end vertex
    return VCLIP_CONTINUE, (cf1 + 1) % (2 * r1.n), cf2

  # Check for local minimum
  if (compute_angle(r2_curr, r1_start, r1_end, r1.is_ccw) < 0):
    # Found local minimum
    dmax = -1
    r1_start = r1.get_v(0)
    for j in range(r1.n):
      # Find the edge with the largest
      # positive distance to the given point
      r1_end = r1.get_v(j + 1)
      if (compute_angle(r2_curr, r1_start, r1_end, r1.is_ccw) > 0):
        dist = compute_dist2(r2_curr, r1_start, r1_end)
        if dist > dmax:
          dmax = dist
          cf1 = 2 * j + 1
      r1_start = r1_end

    # If no positive distance, point is inside polygon
    if dmax == -1:
      print("V-clip: Intersection detected")
      return VCLIP_INTERSECT, None, None
    return VCLIP_CONTINUE, cf1, cf2

  # Compute the point on poly1 closest to v2
  r1_curr = Point.interpolate(r1_start, r1_end, s1)
  # Check if v1 is in v2's Voronoi region
  r2_prev = r2.get_v(v_r2 - 1)
  s2_prev = compute_s(r1_curr, r2_prev, r2_curr)
  if (s2_prev < 1):
    # Go to the previous edge of r2
    return VCLIP_CONTINUE, cf1, (cf2 - 1) % (2 * r2.n)
  r2_next = r2.get_v(v_r2 + 1)
  s2_next = compute_s(r1_curr, r2_curr, r2_next)
  if (s2_next > 0):
    # Go to the next edge of r2
    return VCLIP_CONTINUE, cf1, (cf2 + 1) % (2 * r2.n)
  # We found the closest feature
  return VCLIP_DISJOINT, cf1, cf2


def poly_poly_edge_edge(r1, r2, cf1, cf2):
  v_r1 = cf1 // 2
  v_r2 = cf2 // 2
  # Get edge enpoints
  r1_start = r1.get_v(v_r1)
  r1_end = r1.get_v(v_r1 + 1)
  r2_start = r2.get_v(v_r2)
  r2_end = r2.get_v(v_r2 + 1)

  # Check if the edges intersect
  if (compute_angle(r1_start, r2_start, r2_end) * compute_angle(
      r1_end, r2_start, r2_end) < 0
      and compute_angle(r2_start, r1_start, r1_end) * compute_angle(
          r2_end, r1_start, r1_end) < 0):
    print("V-clip: Intersection detected")
    return VCLIP_INTERSECT, None, None

  # Compute distances of each vertex to opposing edge
  d1_start = compute_point_edge_dist2(r1_start, r2_start, r2_end)
  d1_end = compute_point_edge_dist2(r1_end, r2_start, r2_end)
  d2_start = compute_point_edge_dist2(r2_start, r1_start, r1_end)
  d2_end = compute_point_edge_dist2(r2_end, r1_start, r1_end)

  # Update vertex with smallest distance
  if (d1_start <= d1_end and d1_start <= d2_start and d1_start <= d2_end):
    return VCLIP_CONTINUE, (cf1 - 1) % (2 * r1.n), cf2
  elif (d1_end <= d2_start and d1_end <= d2_end):
    return VCLIP_CONTINUE, (cf1 + 1) % (2 * r1.n), cf2
  elif (d2_start <= d2_end):
    return VCLIP_CONTINUE, cf1, (cf2 - 1) % (2 * r2.n)
  else:
    return VCLIP_CONTINUE, cf1, (cf2 + 1) % (2 * r2.n)


def poly_poly_v_clip(r1, r2, cf1=0, cf2=0):
  loop = 0
  state = VCLIP_CONTINUE
  while state == VCLIP_CONTINUE:
    # Compute next closest feature
    if (cf1 % 2 == 0 and cf2 % 2 == 0):
      state, cf1, cf2 = poly_poly_vertex_vertex(r1, r2, cf1, cf2)
    elif (cf1 % 2 == 0):
      state, cf2, cf1 = poly_poly_edge_vertex(r2, r1, cf2, cf1)
    elif (cf2 % 2 == 0):
      state, cf1, cf2 = poly_poly_edge_vertex(r1, r2, cf1, cf2)
    else:
      state, cf1, cf2 = poly_poly_edge_edge(r1, r2, cf1, cf2)

    loop += 1
    if (loop == VCLIP_MAX_ITERS):
      print(f"V-clip: Cycle detected, current features: {cf1}, {cf2}")
      break
  return cf1, cf2


def mpoint_mpoly_vertex_vertex(mp, mr, t, cf):
  p = mp.at(t)
  v_r = cf // 2
  # Get enpoints of previous end next edge
  r_prev = mr.at_v(v_r - 1, t)
  r_curr = mr.at_v(v_r, t)
  # Check if the point is in r's Voronoi region
  s_prev = compute_s(p, r_prev, r_curr)
  if (s_prev < 1):
    # Go to the previous edge
    return VCLIP_CONTINUE, (cf - 1) % (2 * mr.poly.n)
  r_next = mr.at_v(v_r + 1, t)
  s_next = compute_s(p, r_curr, r_next)
  if (s_next > 0):
    # Go the next edge
    return VCLIP_CONTINUE, (cf + 1) % (2 * mr.poly.n)

  # We found the closest feature
  if (s_prev > 1 or s_next < 0):
    return VCLIP_DISJOINT, cf
  # Point is on the vertex
  print("V-clip: Intersection detected")
  return VCLIP_INTERSECT, None


def mpoint_mpoly_edge_vertex(mp, mr, t, cf):
  p = mp.at(t)
  v_r = cf // 2
  # Get edge enpoints
  r_start = mr.at_v(v_r, t)
  r_end = mr.at_v(v_r + 1, t)
  # Check if the point is in the edge's Voronoin region
  s = compute_s(p, r_start, r_end)
  if (s < 0):
    # Go to the start vertex
    return VCLIP_CONTINUE, (cf - 1) % (2 * mr.poly.n)
  elif (s > 1):
    # Go to the end vertex
    return VCLIP_CONTINUE, (cf + 1) % (2 * mr.poly.n)

  # Check for local minimum
  if (compute_angle(p, r_start, r_end, mr.poly.is_ccw) < EPSILON):
    # Found local minimum
    dmax = -1
    r_start = mr.at_v(0, t)
    for j in range(mr.poly.n):
      # Find the edge with the largest
      # positive distance to the given point
      r_end = mr.at_v(j + 1, t)
      if (compute_angle(p, r_start, r_end, mr.poly.is_ccw) > EPSILON):
        dist = compute_dist2(p, r_start, r_end)
        if dist > dmax:
          dmax = dist
          cf = 2 * j + 1
      r_start = r_end

    # If no positive distance, point is inside polygon
    if dmax == -1:
      print("V-clip: Intersection detected")
      return VCLIP_INTERSECT, None
    return VCLIP_CONTINUE, cf
  # We found the closest feature
  return VCLIP_DISJOINT, cf


def mpoint_mpoly_v_clip(mp, mr, t, cf=0):
  return point_poly_v_clip(mp.at(t), mr.at(t), cf)


def mpoint_mpoly_v_clip_2(mp, mr, t, cf=0):
  loop = 0
  state = VCLIP_CONTINUE
  while state == VCLIP_CONTINUE:
    # Compute next closest feature
    if (cf % 2 == 0):
      state, cf = mpoint_mpoly_vertex_vertex(mp, mr, t, cf)
    else:  # (cf % 2 == 1)
      state, cf = mpoint_mpoly_edge_vertex(mp, mr, t, cf)

    loop += 1
    if (loop == VCLIP_MAX_ITERS):
      print(f"V-clip: Cycle detected, current feature: {cf}")
      break
  return cf


def mpoly_mpoly_vertex_vertex(mr1, mr2, t, cf1, cf2):
  v_r1 = cf1 // 2
  v_r2 = cf2 // 2
  # Get enpoints of previous and next edges
  r1_prev = mr1.at_v(v_r1 - 1, t)
  r1_curr = mr1.at_v(v_r1, t)
  r2_curr = mr2.at_v(v_r2, t)
  # Check if v2 is in v1's Voronoi region
  s1_prev = compute_s(r2_curr, r1_prev, r1_curr)
  if (s1_prev < 1):
    # Go to the previous edge of r1
    return VCLIP_CONTINUE, (cf1 - 1) % (2 * mr1.poly.n), cf2
  r1_next = mr1.at_v(v_r1 + 1, t)
  s1_next = compute_s(r2_curr, r1_curr, r1_next)
  if (s1_next > 0):
    # Go to the previous edge of r1
    return VCLIP_CONTINUE, (cf1 + 1) % (2 * mr1.poly.n), cf2
  r2_prev = mr2.at_v(v_r2 - 1, t)
  r2_next = mr2.at_v(v_r2 + 1, t)
  # Check if v1 is in v2's Voronoi region
  s2_prev = compute_s(r1_curr, r2_prev, r2_curr)
  if (s2_prev < 1):
    # Go to the previous edge of r2
    return VCLIP_CONTINUE, cf1, (cf2 - 1) % (2 * mr2.poly.n)
  s2_next = compute_s(r1_curr, r2_curr, r2_next)
  if (s2_next > 0):
    # Go to the next edge of r2
    return VCLIP_CONTINUE, cf1, (cf2 + 1) % (2 * mr2.poly.n)

  # We found the closest feature
  if ((s1_prev > 1 or s1_next < 0) and (s2_prev > 1 or s2_next < 0)):
    return VCLIP_DISJOINT, cf1, cf2
  # Point is on the vertex
  print("V-clip: Intersection detected")
  return VCLIP_INTERSECT, None, None


def mpoly_mpoly_edge_vertex(mr1, mr2, t, cf1, cf2):
  v_r1 = cf1 // 2
  v_r2 = cf2 // 2
  # Get enpoints of current, previous and next edges
  r1_start = mr1.at_v(v_r1, t)
  r1_end = mr1.at_v(v_r1 + 1, t)
  r2_curr = mr2.at_v(v_r2, t)

  # Check if v2 is in the Voronoi region of the edge of poly1
  s1 = compute_s(r2_curr, r1_start, r1_end)
  if (s1 < 0):
    # Go to the start vertex
    return VCLIP_CONTINUE, (cf1 - 1) % (2 * mr1.poly.n), cf2
  elif (s1 > 1):
    # Go to the end vertex
    return VCLIP_CONTINUE, (cf1 + 1) % (2 * mr1.poly.n), cf2

  # Check for local minimum
  if (compute_angle(r2_curr, r1_start, r1_end, mr1.poly.is_ccw) < 0):
    # Found local minimum
    dmax = -1
    r1_start = mr1.at_v(0, t)
    for j in range(mr1.poly.n):
      # Find the edge with the largest
      # positive distance to the given point
      r1_end = mr1.at_v(j + 1, t)
      if (compute_angle(r2_curr, r1_start, r1_end, mr1.poly.is_ccw) > 0):
        dist = compute_dist2(r2_curr, r1_start, r1_end)
        if dist > dmax:
          dmax = dist
          cf1 = 2 * j + 1
      r1_start = r1_end

    # If no positive distance, point is inside polygon
    if dmax == -1:
      print("V-clip: Intersection detected")
      return VCLIP_INTERSECT, None, None
    return VCLIP_CONTINUE, cf1, cf2

  # Compute the point on poly1 closest to v2
  r1_curr = Point.interpolate(r1_start, r1_end, s1)
  # Check if v1 is in v2's Voronoi region
  r2_prev = mr2.at_v(v_r2 - 1, t)
  s2_prev = compute_s(r1_curr, r2_prev, r2_curr)
  if (s2_prev < 1):
    # Go to the previous edge of r2
    return VCLIP_CONTINUE, cf1, (cf2 - 1) % (2 * mr2.poly.n)
  r2_next = mr2.at_v(v_r2 + 1, t)
  s2_next = compute_s(r1_curr, r2_curr, r2_next)
  if (s2_next > 0):
    # Go to the next edge of r2
    return VCLIP_CONTINUE, cf1, (cf2 + 1) % (2 * mr2.poly.n)
  # We found the closest feature
  return VCLIP_DISJOINT, cf1, cf2


def mpoly_mpoly_edge_edge(mr1, mr2, t, cf1, cf2):
  v_r1 = cf1 // 2
  v_r2 = cf2 // 2
  # Get edge enpoints
  r1_start = mr1.at_v(v_r1, t)
  r1_end = mr1.at_v(v_r1 + 1, t)
  r2_start = mr2.at_v(v_r2, t)
  r2_end = mr2.at_v(v_r2 + 1, t)

  # Check if the edges intersect
  if (compute_angle(r1_start, r2_start, r2_end) * compute_angle(
      r1_end, r2_start, r2_end) < 0
      and compute_angle(r2_start, r1_start, r1_end) * compute_angle(
          r2_end, r1_start, r1_end) < 0):
    print("V-clip: Intersection detected")
    return VCLIP_INTERSECT, None, None

  # Compute distances of each vertex to opposing edge
  d1_start = compute_point_edge_dist2(r1_start, r2_start, r2_end)
  d1_end = compute_point_edge_dist2(r1_end, r2_start, r2_end)
  d2_start = compute_point_edge_dist2(r2_start, r1_start, r1_end)
  d2_end = compute_point_edge_dist2(r2_end, r1_start, r1_end)

  # Update vertex with smallest distance
  if (d1_start <= d1_end and d1_start <= d2_start and d1_start <= d2_end):
    return VCLIP_CONTINUE, (cf1 - 1) % (2 * mr1.poly.n), cf2
  elif (d1_end <= d2_start and d1_end <= d2_end):
    return VCLIP_CONTINUE, (cf1 + 1) % (2 * mr1.poly.n), cf2
  elif (d2_start <= d2_end):
    return VCLIP_CONTINUE, cf1, (cf2 - 1) % (2 * mr2.poly.n)
  else:
    return VCLIP_CONTINUE, cf1, (cf2 + 1) % (2 * mr2.poly.n)


def mpoly_mpoly_v_clip(mr1, mr2, t, cf1=0, cf2=0):
  return poly_poly_v_clip(mr1.at(t), mr2.at(t), cf1, cf2)


def mpoly_mpoly_v_clip_2(mr1, mr2, t, cf1=0, cf2=0):
  loop = 0
  state = VCLIP_CONTINUE
  while state == VCLIP_CONTINUE:
    # Compute next closest feature
    if (cf1 % 2 == 0 and cf2 % 2 == 0):
      state, cf1, cf2 = mpoly_mpoly_vertex_vertex(mr1, mr2, t, cf1, cf2)
    elif (cf1 % 2 == 0):
      state, cf2, cf1 = mpoly_mpoly_edge_vertex(mr2, mr1, t, cf2, cf1)
    elif (cf2 % 2 == 0):
      state, cf1, cf2 = mpoly_mpoly_edge_vertex(mr1, mr2, t, cf1, cf2)
    else:
      state, cf1, cf2 = mpoly_mpoly_edge_edge(mr1, mr2, t, cf1, cf2)

    loop += 1
    if (loop == VCLIP_MAX_ITERS):
      print(f"V-clip: Cycle detected, current features: {cf1}, {cf2}")
      break
  return cf1, cf2