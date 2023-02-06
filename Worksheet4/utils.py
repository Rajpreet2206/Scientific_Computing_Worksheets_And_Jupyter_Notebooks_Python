import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams, patches
import time
from particle import *


def potential_function(x, y):
    return np.log(np.sqrt(x**2 + y**2))


def direct_outside_potential(coordinates, list_particles):
    potential = 0.
    for particle in list_particles:
        x, y = np.abs(particle.x - coordinates[0]), np.abs(particle.y - coordinates[1])
        potential += potential_function(x, y)
    return potential


def get_quad_points(cell_center, a, K):
    quad_points = []
    delta_theta = 2 * np.pi / K
    for i in range(K):
        x, y = a*np.cos(i*delta_theta) + cell_center[0], a*np.sin(i*delta_theta) + cell_center[1]
        quad_points.append([x,y])
    return quad_points
