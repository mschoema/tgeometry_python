from math import pi
from random import uniform

from .constants import EPSILON


class Angle():
  """
    Simple class representing an angle theta in ]-pi, pi]
    """

  @staticmethod
  def interpolate(theta_1, theta_2, ratio):
    theta_delta = theta_2 - theta_1
    # If abs(theta_delta) == pi: Always turn counter-clockwise
    if (abs(theta_delta) < EPSILON):
      theta = theta_1
    elif (theta_delta > 0 and abs(theta_delta) <= pi):
      theta = theta_1 + theta_delta * ratio
    elif (theta_delta > 0 and abs(theta_delta) > pi):
      theta = theta_2 + (2 * pi - theta_delta) * (1 - ratio)
    elif (theta_delta < 0 and abs(theta_delta) < pi):
      theta = theta_1 + theta_delta * ratio
    else:  # (theta_delta < 0 and abs(theta_delta) >= pi)
      theta = theta_1 + (2 * pi + theta_delta) * ratio
    if theta > pi:
      theta = theta - 2 * pi
    return theta

  # Generate a random angle in ]-pi, pi]
  # Round decimal values for readability
  @staticmethod
  def random(digits=2):
    theta = uniform(-pi, pi)
    if theta == -pi:
      theta = pi
    return round(theta, digits)
