#
# Some geometric primitives
#

from math import sqrt, atan2
from functools import reduce
import re,sys

assert sys.hexversion >= 0x3000000

class V():
  def __init__(self, x, y):
    self.x = float(x)
    self.y = float(y)

  def __str__(self):
    return "(%g,%g)" % (self.x, self.y)

  def ipe(self, suffix = ""):
    return "%g %g%s" % (self.x, self.y, suffix)

  def __repr__(self):
    return "V(%g, %g)" % (self.x, self.y)

  def __sub__(self, rhs):
    return V(self.x - rhs.x, self.y - rhs.y)

  def __add__(self, rhs):
    return V(self.x + rhs.x, self.y + rhs.y)

  def __neg__(self):
    return V(-self.x, -self.y)

  # vector * scalar
  def __mul__(self, rhs):
    return V(self.x * rhs, self.y * rhs)

  # scalar * vector
  def __rmul__(self, lhs):
    return V(self.x * lhs, self.y * lhs)

  # vector / scalar
  def __truediv__(self, rhs):
    return V(self.x / rhs, self.y / rhs)

  def __eq__(self, rhs):
    return self.x == rhs.x and self.y == rhs.y

  def __lt__(self, rhs):
    return (self.x < rhs.x) or (self.x == rhs.x and self.y < rhs.y)

  def __le__(self, rhs):
    return not(rhs < self)

  def sqlen(self):
    return self.x**2 + self.y**2

  def norm(self):
    return sqrt(self.sqlen())

  def __matmul__(self, rhs):
    return self.x * rhs.x + self.y * rhs.y

  def cross(self, rhs):
    return self.x * rhs.y - self.y * rhs.x

  def normalized(self):
    n = self.norm()
    return V(self.x / n, self.y / n)

  def flipped(self):
    return V(self.x, -self.y)

# --------------------------------------------------------------------

class Line():
  def __init__(self, p, q):
    self.p = p
    self.q = q

  def __str__(self):
    return "%s--%s" % (self.p, self.q)

  def __repr__(self):
    return "Line(%s, %s)" % (repr(self.p), repr(self.q))

  def dir(self):
    return (self.q - self.p).normalized()

  def slope(self):
    d = self.dir()
    if d.x < 0: d = -d
    return atan2(d.y, d.x)

  # comparison is by slope only!
  def __eq__(self, rhs):
    return self.dir() == rhs.dir()

  def __lt__(self, rhs):
    return self.slope() < rhs.slope()

  # which side of this line does x lie on? < 0 ==  left
  def side(self, x):
    return (x - self.p).cross(self.q - self.p)
    
  def intersect(self, rhs):
    denom = rhs.dir().cross(self.dir())
    if denom == 0.0:
      return None
    lam = (self.p - rhs.p).cross(rhs.dir()) / denom
    return self.p + self.dir() * lam

  # is this line occluded by lo and hi (which bracket this line slopewise)
  def occluded(self, lo, hi):
    p = lo.intersect(hi)
    return self.side(p) <= 0

  def flipped(self):
    return Line(self.p.flipped(), self.q.flipped())

# --------------------------------------------------------------------

# Halfspace h represented as a * x <= b
class Halfspace():
  # init with a directed line, hs lies to its left
  def __init__(self, a, b = None):
    if isinstance(a, Line):
      d = a.q - a.p
      self.a = V(d.y, -d.x)  # rotate right 90 degrees
      self.b = self.a @ a.p
    else:
      self.a = a
      self.b = b

  def __repr__(self):
    return "Halfspace(%s, %g)" % (self.a, self.b)

  def normalized(self):
    """Normalized halfspace where |a| == 1"""
    n = self.a.norm()
    return Halfspace(self.a / n, self.b / n)

  def uniform(self):
    """Same halfspace with b == 1"""
    assert self.b > 0, self
    return Halfspace(self.a / self.b, 1.)

  def __contains__(self, p):
    return (self.a @ p <= self.b)

  def side(self, p):
    return (self.a @ p - self.b)
  

# --------------------------------------------------------------------
#   Convex polygons
# --------------------------------------------------------------------

class EmptyPolygon():
  def __init__(self):
    pass
  def __repr__(self):
    return "EmptyPolygon()"
  def area(self):
    return 0.0
  def transform(self, lam, t):
    return self
  def draw(self, ipe, f = 64):
    pass  # not drawn at all

# --------------------------------------------------------------------

class Polygon():
  # v is a list of vertices in ccw order
  def __init__(self, v):
    self.v = v[:]   # make a copy
    # ensure ccw order
    if (self.v[1] - self.v[0]).cross(self.v[2] - self.v[0]) < 0:
      self.v.reverse()
    self.v0 = min(self.v)
    self.vn = max(self.v)
    i = self.v.index(self.v0)
    self.v = self.v[i:] + self.v[:i]
    self.hs = []
    for i in range(1, len(self.v)):
      self.hs.append(Halfspace(Line(self.v[i-1], self.v[i])))
    self.hs.append(Halfspace(Line(self.v[-1], self.v0)))

  def __repr__(self):
    return "Polygon(%s)" % repr(self.v)
    
  def area(self):
    total = 0.0
    for i in range(2, len(self.v)):
      total += (self.v[i-1] - self.v0).cross(self.v[i] - self.v0)
    return total / 2.0
      
  def transform(self, lam=1.0, t=V(0,0)):
    v = list(map(lambda x: lam * x + t, self.v))
    return Polygon(v)

  def chains(self):
    i = self.v.index(self.vn)
    lo = self.v[:i+1]
    hi = self.v[i:] + [self.v0]
    hi.reverse()
    return lo, hi

  def draw(self, ipe, f = 64):
    ipe.write("<path>\n")
    ipe.write((self.v0 * f).ipe(" m\n"))
    for i in range(1, len(self.v)):
      ipe.write((self.v[i] * f).ipe(" l\n"))
    ipe.write("h</path>\n")

# --------------------------------------------------------------------

def start_ipe(fname):
  out = open(fname, "w")
  out.write("""<ipe version="70000">
<ipestyle name="basic">
</ipestyle>
<page>""")
  sys.stderr.write("Writing ipe file '%s' ..." % fname)
  return out

def stop_ipe(out):
  out.write("</page></ipe>\n")
  out.close()
  sys.stderr.write("written\n")

# --------------------------------------------------------------------

def _parse_vertex(s):
  f = s.split(maxsplit=2)
  return V(f[0], f[1])
  
def parse_polygons(fname, factor = 64):
  data = open(fname, "r").read()
  i = data.find("<page")
  data = data[i:]
  polygons = []
  path_objs = re.findall("<path[^>]*>\\s*(\\S[^<]+\\S)\\s*l\\s+h\\s*</path>",
                         data)
  for path in path_objs:
    coords = re.split("\\s*[mlh]\\s*", path.strip())
    v = map(lambda s: _parse_vertex(s) / factor, coords)
    polygons.append(Polygon(list(v)))
  return polygons

# --------------------------------------------------------------------

def _turn_left(p, q, r):
  return (q-p).cross(r-p) > 0

def _keep_left(hull, r):
  while len(hull) > 1 and not _turn_left(hull[-2], hull[-1], r):
    hull.pop()
  hull.append(r)
  return hull

def convex_hull(pts):
  pts.sort()
  l = reduce(_keep_left, pts, [])
  pts.reverse()
  u = reduce(_keep_left, pts, [])
  l.extend(u[1:-1])
  return l

# --------------------------------------------------------------------
