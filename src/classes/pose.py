from .angle import Angle
from .point import Point


class Pose():
  """
    Simple class representing a 2d pose
    """

  __slots__ = ['pos', 'theta']

  def __init__(self, pos=Point(), theta=0):
    self.pos = pos
    self.theta = theta

  @staticmethod
  def interpolate(pose_1, pose_2, ratio):
    pos = Point.interpolate(pose_1.pos, pose_2.pos, ratio)
    theta = Angle.interpolate(pose_1.theta, pose_2.theta, ratio)
    return Pose(pos, theta)

  # Generate a random pose inside the given box
  # If linear = True, the rotation angle theta will be 0
  # Round decimal values for readability
  @staticmethod
  def random(xmin, ymin, xmax, ymax, linear=False, digits=2):
    pos = Point.random(xmin, ymin, xmax, ymax, digits)
    if linear:
      theta = 0
    else:
      theta = Angle.random(digits)
    return Pose(pos, theta)

  def __repr__(self):
    return f"(pos={repr(self.pos)}, theta={repr(self.theta)})"

  def __str__(self):
    return f"Pose(pos={repr(self.pos)}, theta={repr(self.theta)})"