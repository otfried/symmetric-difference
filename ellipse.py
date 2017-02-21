#
# Experiments for finding the best ellipse for a given polygon
#

from geometry import *
import sys
from cvxopt import matrix, solvers
from scipy.optimize import minimize

import numpy as np

assert sys.hexversion >= 0x3000000

solvers.options['show_progress'] = False

# --------------------------------------------------------------------
#   Convex polygon intersection
# --------------------------------------------------------------------

def intersect_polygons(P, Q):
  # set up LP for intersection
  p = len(P.hs)
  q = len(Q.hs)
  G = matrix(0., (p + q, 3))
  h = matrix(0., (p + q, 1))
  for i in range(p):
    h[i] = P.hs[i].b
    G[i, 0] = P.hs[i].a.x
    G[i, 1] = P.hs[i].a.y
    G[i, 2] = 1.
  for j in range(q):
    h[p + j] = Q.hs[j].b
    G[p + j, 0] = Q.hs[j].a.x
    G[p + j, 1] = Q.hs[j].a.y
    G[p + j, 2] = 1.
  c = matrix([0., 0., -1.])  # maximize "distance" to boundary
  sol = solvers.lp(c, G, h)
  assert sol['status'] == 'optimal'
  x = sol['x']
  t = V(x[0], x[1])
  if x[2] <= 0:
    return EmptyPolygon()
  # dualize the problem
  pts = []
  for hs in P.hs + Q.hs:
    n = hs.b - hs.a @ t
    assert n > 0, hs
    pts.append(hs.a / n)
  hull = convex_hull(pts)
  v = []
  for i in range(1, len(hull)):
    # dualize lines between consecutive corners
    v.append(Halfspace(Line(hull[i-1], hull[i])).uniform().a)
  v.append(Halfspace(Line(hull[-1], hull[0])).uniform().a)
  return Polygon(v).transform(1.0, t)

# --------------------------------------------------------------------

def convex_union(P, Q):
  pts = convex_hull(P.v + Q.v)
  return Polygon(pts)

# --------------------------------------------------------------------

def max_inscribed(P, Q):
  m = len(P.hs)
  n = m * len(Q.v)        # number of constraints
  # solve min c @ x s.t. G x <= h
  c = matrix([-1., 0., 0.])
  G = matrix(0., (n, 3))
  h = matrix(0., (n, 1))
  for i in range(len(P.hs)):
    for j in range(len(Q.v)):
      row = j * m + i
      h[row] = P.hs[i].b
      G[row, 0] = Q.v[j] @ P.hs[i].a
      G[row, 1] = P.hs[i].a.x
      G[row, 2] = P.hs[i].a.y
  sol = solvers.lp(c, G, h)
  assert sol['status'] == 'optimal'
  x = sol['x']
  lam = x[0]
  t = V(x[1], x[2])
  return lam, t

# --------------------------------------------------------------------

def min_circumscribed(P, Q):
  m = len(Q.hs)
  n = m * len(P.v)        # number of constraints
  # solve min c @ x s.t. G x <= h
  c = matrix([1., 0., 0.])
  G = matrix(0., (n, 3))
  h = matrix(0., (n, 1))
  for i in range(len(Q.hs)):
    for j in range(len(P.v)):
      row = j * m + i
      h[row] = -(Q.hs[j].a @ P.v[i])
      G[row, 0] = -Q.hs[j].b
      G[row, 1] = -Q.hs[j].a.x
      G[row, 2] = -Q.hs[j].a.y
  sol = solvers.lp(c, G, h)
  assert sol['status'] == 'optimal'
  x = sol['x']
  lam = x[0]
  t = V(x[1], x[2])
  return lam, t

# ====================================================================

def load_poly(fname):
  Ps = parse_polygons(fname)
  print("Area(P) = %g" % Ps[0].area())
  return Ps[0]

def draw_polys(P, Q, R = None):
  out = start_ipe("test.ipe")
  P.draw(out)
  Q.draw(out)
  if R:
    R.draw(out)
  stop_ipe(out)
  
# --------------------------------------------------------------------

def fun(x, P, Q):
  Q1 = Q.transform(x[0], V(x[1], x[2]))
  return Q1.area() - 2 * intersect_polygons(P, Q1).area()

def optimize(P, Q):
  lam1, t1 = max_inscribed(P, Q)
  # print("Max inscribed lambda = %g" % lam1)
  x0 = [ lam1, t1.x, t1.y ]
  res = minimize(fun, x0, args=(P, Q),
                 method='Nelder-Mead',
                 options={'disp': False})
  if not res.success:
    sys.stderr.write("Failed to optimize!\n")
  return res.x[0], V(res.x[1], res.x[2]), res.fun

# --------------------------------------------------------------------

def main():
  if len(sys.argv) != 2:
    sys.stderr.write("Usage: python3 ellipse.py <ipefile>\n")
    sys.exit(9)
  fname = sys.argv[1]
  P = load_poly(fname)
  Q = generate_kgon(37)
  shears = np.linspace(-1.0, 1.0, 11)
  stretches = np.linspace(0.4, 2.0, 9)
  values = np.empty([11, 9])
  best_shear = None
  best_stretch = None
  best_lam = None
  best_t = None
  best_val = 9999999.0
  for i in range(len(shears)):
    shear = shears[i]
    for j in range(len(stretches)):
      stretch = stretches[j]
      Q1 = Q.affine_transform(1.0, shear, 0.0, stretch)
      lam, t, val = optimize(P, Q1)
      print("%d %d %g %g: %g" % (i, j, shear, stretch, val))
      values[i, j] = val
      if val < best_val:
        best_stretch = stretch
        best_shear = shear
        best_lam = lam
        best_t = t
        best_val = val
  print(values)
  print(best_shear, best_stretch, best_lam, best_t)
  Q2 = Q.affine_transform(1.0, best_shear,
                          0.0, best_stretch).transform(best_lam, best_t)
  R = intersect_polygons(P, Q2)
  draw_polys(P, Q2, R)

# --------------------------------------------------------------------

if __name__ == "__main__":
  main()

# --------------------------------------------------------------------
