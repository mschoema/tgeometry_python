from .constants import EPSILON
from .poly import Poly
from .pose import Pose


class MPoly():
    """
    Simple class representing a 2d moving polygon
    """

    __slots__ = ['poly', 'start_pose', 'end_pose', 'linear']

    def __init__(self, poly, start_pose, end_pose):
        self.poly = poly
        self.start_pose = start_pose
        self.end_pose = end_pose
        self.linear = abs(start_pose.theta - end_pose.theta) < EPSILON

    def at(self, t):
        pose = Pose.interpolate(self.start_pose, self.end_pose, t)
        return self.poly.apply(pose)

    def at_v(self, i, t):
        pose = Pose.interpolate(self.start_pose, self.end_pose, t)
        return self.poly.get_v(i).apply(pose, self.poly.c)

    # Generate a random moving n-gon
    # The polyogn is created inside the box (-wmax, -hmax, wmax, hmax)
    # The start and end poses are created inside the box (xmin, ymin, xmax, ymax)
    # The polygon can thus intersect the (xmin, ymin, xmax, ymax) box during its movement
    # If linear = True, the rotation angle theta will be 0
    # Round decimal values for readability
    @staticmethod
    def random(xmin,
               ymin,
               xmax,
               ymax,
               wmax,
               hmax,
               n,
               linear=False,
               digits=2):
        poly = Poly.random(-wmax, -hmax, wmax, hmax, n, digits)
        p1 = Pose.random(xmin, ymin, xmax, ymax, linear, digits)
        p2 = Pose.random(xmin, ymin, xmax, ymax, linear, digits)
        return MPoly(poly, p1, p2)

    def __repr__(self):
        return f"({repr(self.poly)}, {repr(self.start_pose)}, {repr(self.end_pose)})"

    def __str__(self):
        return f"MPoly({self.poly}, {self.start_pose}, {self.end_pose})"