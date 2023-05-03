from .point import Point


class MPoint():
  """
    Simple class representing a 2d moving point
    """

  __slots__ = ['start_point', 'end_point']

  def __init__(self, start_point, end_point):
    self.start_point = start_point
    self.end_point = end_point

  def at(self, t):
    return Point.interpolate(self.start_point, self.end_point, t)

  # Generate a random moving point inside the given box
  # Round decimal values for readability
  @staticmethod
  def random(xmin, ymin, xmax, ymax, digits=2):
    p1 = Point.random(xmin, ymin, xmax, ymax, digits)
    p2 = Point.random(xmin, ymin, xmax, ymax, digits)
    return MPoint(p1, p2)

  def __repr__(self):
    return f"({repr(self.start_point)}, {repr(self.end_point)})"

  def __str__(self):
    return f"MPoint({self.start_point}, {self.end_point})"