import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams, patches
import utils
import math


class Cell():
    
    """The class for a cell.
    
    Arguments:
        n_crit: maximum number of particles in a leaf cell.
    
    Attributes:
        nparticles (int): number of particles in the cell
        particle_index (array of int): array of particles
        nchild (int):  an integer whose last 8 bits is used to keep track 
        of the empty child cells
        child_index (array of int): array of child index
        parent (int): index of parent cell
        x, y (float): coordinates of the cell's center
        r (float): radius of the cell (half of the side length for cubic cell)      
    """
    
    def __init__(self, n_crit, center=[0.5, 0.5], side_length=1., M=3):
        self.nparticles = 0        # number of particels
        self.particle_index = []     # array of particle indices
        self.nchild = 0       # binary counter to keep track of empty cells
        self.child_index = np.zeros(4, dtype=np.int)         # array of child index
        self.parent = 0       # index of parent cell
        self.x = center[0]
        self.y = center[1]                    # center of the cell
        self.r = side_length / 2           # radius of the cell
        self.n_crit = n_crit
        K = 2*M +1
        self.K = K
        self.a = side_length
        self.quad_points = utils.get_quad_points(center, self.a, K)
        self.quad_angles = np.array([i*2*np.pi/K for i in range(K)])
        self.inner_ring_expansion_terms = np.zeros(K, dtype=np.float) # potential at quadrature points due to particles outside the box
        self.outer_ring_expansion_terms = np.zeros(K, dtype=np.float) # potential at quadrature points due to particles inside the box
        self.kappa = None
        self.h = 2 * np.pi * self.a/ K

    def distance(self, other):
        return np.sqrt((self.x-other.x)**2 + (self.y-other.y)**2)
    
    def calculate_total_mass(self, list_cells, list_particles):
        mass = 0.
        if self.nparticles < self.n_crit:
            for p_idx in self.particle_index:
                mass += list_particles[p_idx].m
        else:
            for c_idx in self.child_index:
                if c_idx != 0:
                    mass += list_cells[c_idx].calculate_total_mass(list_cells, list_particles)
        return mass
    
    def eval_total_mass(self, list_cells, list_particles):
        self.kappa = self.calculate_total_mass(list_cells, list_particles)
    
    def calculate_outer_ring_terms_leaf(self, list_particles):
        assert self.nparticles < self.n_crit, "This function needs to be called only for leaf cells"
        section_list = np.array(list_particles)[self.particle_index]
        for i in range(self.K):
            self.outer_ring_expansion_terms[i] = utils.direct_outside_potential(self.quad_points[i], section_list) 
    
    def calculate_outer_ring_from_child(self, list_cells):
        self.outer_ring_expansion_terms = np.zeros(self.K, dtype=np.float)
        for k in range(self.K):
            for c_idx in self.child_index:
                child = list_cells[c_idx]
                dx, dy = self.quad_points[k][0] - child.x, self.quad_points[k][1] - child.y
                r = np.sqrt(dx**2 + dy**2)
                theta = math.atan2(dy,dx)
                self.outer_ring_expansion_terms[k] += child.outer_ring_expansion(r, theta)
    
    def outer_ring_expansion(self, r, theta):
        assert self.kappa is not None, "Runtime Error: Total mass of individual cell is not yet calculate"
        ################# Task 1 begins #########################

        ################## Task 1 begins ########################
        return potential
    
    def inner_ring_expansion(self, r, theta):
        ################ Task 2 begins #######################

        ################ Task 2 ends ########################
        return potential

    def plot(self, ax, linewidth=1, edgecolor='r', facecolor='none'):
        rect = patches.Rectangle((self.x-self.r, self.y-self.r), 2*self.r, 2*self.r, linewidth=linewidth,
                                 edgecolor=edgecolor, facecolor=facecolor)
        ax.add_patch(rect)
        
    def plot_quad_points(self, color='b', marker='o'):
        for i in range(self.K):
            plt.plot(self.quad_points[i][0], self.quad_points[i][1], color=color, marker=marker)
    
    def plot_ring(self, ax, linewidth=1, edgecolor='b', facecolor='none'):
        circle = patches.Circle((self.x, self.y), radius=self.a, linewidth=linewidth, edgecolor=edgecolor, facecolor=facecolor)
        ax.add_patch(circle)
        
    def plot_cell_and_ring(self, ax):
        self.plot_ring(ax)
        self.plot(ax)
        self.plot_quad_points()
        
    def plot1(self, n_crit, linewidth=1, color='k'):
        if self.nparticles >= n_crit:
            plt.hlines(self.y, self.x - self.r, self.x + self.r, linewidth=linewidth, color=color)
            plt.vlines(self.x, self.y - self.r, self.y + self.r, linewidth=linewidth, color=color)
