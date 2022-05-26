from math import atan2
from random import randint
from random import random
from random import shuffle
from random import uniform

from .point import Point


class Poly():
    """
    Simple class representing a 2d polygon

    The vertices of the polygon must form a closed loop
    (i.e.: pts[0] == pts[-1])
    """

    __slots__ = ['pts', 'c', 'n', 'is_ccw']

    def __init__(self, pts, c=Point(0, 0)):
        self.pts = pts
        self.c = c
        self.n = len(pts) - 1
        self.is_ccw = self.pts_ccw(pts)

    def get_v(self, i):
        return self.pts[i % self.n]

    def apply(self, pose):
        pts = [p.apply(pose, self.c) for p in self.pts]
        c = self.c.apply(pose, self.c)
        return Poly(pts, c)

    # Test if the polygon is counter-clockwise
    # Checks the angle between three consecutive points
    # We can chose any 3 consecutive points,
    # since we are working with convex polygons.
    @staticmethod
    def pts_ccw(pts):
        tot = 0
        for i in range(len(pts) - 2):
            a = pts[i]
            b = pts[i + 1]
            c = pts[i + 2]
            tot += (a.x - b.x) * (c.y - b.y) - (a.y - b.y) * (c.x - b.x)
        return tot < 0

    # Computes the centroid of a list of points representing a polygon
    @staticmethod
    def pts_centroid(pts):
        n = len(pts) - 1

        area = 0
        for i in range(n):
            p1 = pts[i]
            p2 = pts[i + 1]
            area += p1.x * p2.y - p2.x * p1.y
        area /= 2

        cx = cy = 0
        for i in range(n):
            p1 = pts[i]
            p2 = pts[i + 1]
            cx += (p1.x + p2.x) * (p1.x * p2.y - p2.x * p1.y)
            cy += (p1.y + p2.y) * (p1.x * p2.y - p2.x * p1.y)
        cx /= 6 * area
        cy /= 6 * area

        return Point(cx, cy)

    # Generate a random n-gon inside the given box
    # If center=True, the polygon will have its centroid (approximately) at (0, 0)
    # Round decimal values for readability
    @staticmethod
    def random(xmin, ymin, xmax, ymax, n, center=True, digits=2):
        def gen_random_poly(xmin, ymin, xmax, ymax, n):
            xpool = []
            ypool = []

            for i in range(n):
                xpool.append(uniform(xmin, xmax))
                ypool.append(uniform(ymin, ymax))

            xpool.sort()
            ypool.sort()

            xmin = xpool[0]
            xmax = xpool[-1]
            ymin = ypool[0]
            ymax = ypool[-1]

            xvec = []
            yvec = []

            lasttop = lastbot = xmin
            lastleft = lastright = ymin

            for i in range(1, n - 1):
                x = xpool[i]
                if random() < 0.5:
                    xvec.append(x - lasttop)
                    lasttop = x
                else:
                    xvec.append(lastbot - x)
                    lastbot = x

                y = ypool[i]
                if random() < 0.5:
                    yvec.append(y - lastleft)
                    lastleft = y
                else:
                    yvec.append(lastright - y)
                    lastright = y

            xvec.append(xmax - lasttop)
            xvec.append(lastbot - xmax)
            yvec.append(ymax - lastleft)
            yvec.append(lastright - ymax)

            shuffle(yvec)

            vec = list(zip(xvec, yvec))
            vec.sort(key=lambda v: atan2(v[1], v[0]))

            x = y = 0
            minPolyX = minPolyY = 0
            points = []

            for i in range(n):
                points.append((x, y))

                x += vec[i][0]
                y += vec[i][1]

                minPolyX = min(x, minPolyX)
                minPolyY = min(y, minPolyY)

            xshift = xmin - minPolyX
            yshift = ymin - minPolyY

            for i in range(n):
                points[i] = (round(points[i][0] + xshift, digits),
                             round(points[i][1] + yshift, digits))

            points.append(points[0])

            duplicates = False
            for i in range(n - 1):
                p1 = points[i]
                p2 = points[i + 1]
                if p1[0] == p2[0] and p1[1] == p2[1]:
                    duplicates = True
                    break

            collinear = False
            for i in range(n - 2):
                p1 = points[i]
                p2 = points[i + 1]
                p3 = points[i + 2]
                if (p1[0] - p2[0]) * (p2[1] - p3[1]) == (p1[1] - p2[1]) * (
                        p2[0] - p3[0]):
                    collinear = True
                    break

            convex = True
            a = points[i]
            b = points[i + 1]
            c = points[i + 2]
            ccw = (a[0] - b[0]) * (c[1] - b[1]) - (a[1] - b[1]) * (
                c[0] - b[0]) < 0
            for i in range(1, n - 2):
                a = points[i]
                b = points[i + 1]
                c = points[i + 2]
                if ccw and (a[0] - b[0]) * (c[1] - b[1]) - (a[1] - b[1]) * (
                        c[0] - b[0]) > 0:
                    convex = False
                    break
                elif not ccw and (a[0] - b[0]) * (c[1] - b[1]) - (
                        a[1] - b[1]) * (c[0] - b[0]) < 0:
                    convex = False
                    break
            return points, duplicates, collinear, convex

        duplicates = True
        collinear = True
        convex = False
        while duplicates or collinear or not convex:
            points, duplicates, collinear, convex = gen_random_poly(
                xmin, ymin, xmax, ymax, n)

        pts = [Point(p[0], p[1]) for p in points]
        c = Poly.pts_centroid(pts)
        if center:
            for i in range(len(pts)):
                pts[i] = Point(
                    round(pts[i].x - round(c.x, digits), digits),
                    round(pts[i].y - round(c.y, digits), digits))
            return Poly(pts)
        else:
            return Poly(pts, c)

    def __repr__(self):
        return f"[c={repr(self.c)}, pts={self.pts}]"

    def __str__(self):
        return f"Poly(c={self.c}, pts={self.pts})"