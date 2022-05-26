from math import cos
from math import sin
from math import sqrt
from random import uniform


class Point():
    """
    Simple class representing a 2d point
    """

    __slots__ = ['x', 'y']

    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y

    @staticmethod
    def interpolate(p_1, p_2, ratio):
        x = p_1.x * (1 - ratio) + p_2.x * ratio
        y = p_1.y * (1 - ratio) + p_2.y * ratio
        return Point(x, y)

    def apply(self, pose, c):
        co = cos(pose.theta)
        si = sin(pose.theta)
        x = (self.x - c.x) * co - (self.y - c.y) * si + c.x + pose.pos.x
        y = (self.x - c.x) * si + (self.y - c.y) * co + c.y + pose.pos.y
        return Point(x, y)

    def dist(self, other):
        return sqrt(self.dist2(other))

    def dist2(self, other):
        return pow(self.x - other.x, 2) + pow(self.y - other.y, 2)

    # Generate a random point inside the given box
    # Round decimal values for readability
    @staticmethod
    def random(xmin, ymin, xmax, ymax, digits=2):
        x = uniform(xmin, xmax)
        y = uniform(ymin, ymax)
        return Point(round(x, digits), round(y, digits))

    def __repr__(self):
        return f"({self.x}, {self.y})"

    def __str__(self):
        return f"Point({self.x}, {self.y})"