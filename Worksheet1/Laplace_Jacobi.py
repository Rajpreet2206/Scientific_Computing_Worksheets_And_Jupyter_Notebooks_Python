%matplotlib inline
import numpy as np
# Print with smaller accurary to make comparisons easier
np.set_printoptions(precision=7, suppress=True)

import matplotlib.pyplot as plt
import scipy.sparse as sp
class Grid(object):
    def __init__(self, num_points):
        # Num points does not include boundary!
        self.num_points = num_points
        self.spacing = 1.0/(num_points + 1)
        
    def get_nth_point(self, i):
        return self.spacing * (i + 1)
        

def make_grid(level):
    # Level is refinement
    # Consider interval [0,1]
    # Level one has one point at 1/2
    # Level two has points at 1/4, 1/2, 3/4
    num_points = 2**(level + 1) - 1
    return Grid(num_points)


grid = make_grid(level=2)
grid.get_nth_point(1)

def create_matrix(grid):
    n = grid.num_points
    h = grid.spacing
    
    def val(i,j):
        if i == j:
            return 2/h**2
        elif np.abs(i-j) == 1:
            return -1/h**2
        else:
            return 0.0
        
    return np.fromfunction(np.vectorize(val), (n, n), dtype=float)


discretization = create_matrix(grid)
discretization

def analytical_eigenvalue(grid, k):
    h = grid.spacing
    c = k * np.pi
    return (4/h**2) * np.sin(c * grid.spacing * 0.5)**2

grid = make_grid(2)
discretization = create_matrix(grid)

numerical_eigenvals = np.linalg.eig(discretization)[0]

analytical_eigenvals = np.array([analytical_eigenvalue(grid,k) for k in range(1,grid.num_points+1)])
print("Numerical:", np.sort(numerical_eigenvals))
print("Analytical:", np.sort(analytical_eigenvals))
# They are in different order
diff = np.sort(numerical_eigenvals) - np.sort(analytical_eigenvals)
print("Difference:", diff)

# The numerical ones are:
print(np.linalg.eig(discretization)[1])

def analytical_eigenvector(grid, k):
    n = grid.num_points
    h = grid.spacing
    c = k * np.pi
    
    eigenvec = np.array([np.sin(c * i * h) for i in range(1, n+1)])
    # Note: numpy always normalizes eigenvectors to unit length
    # So we also have to normalize the analytical ones
    return eigenvec/np.linalg.norm(eigenvec)

for k in range(0, grid.num_points):
    print(analytical_eigenvector(grid, k+1))

"""
Returns a tuple of (iteration matrix, inverse diagonal part)
"""
def build_jacobi_matrices(grid, matrix, omega):
    n = grid.num_points
    h = grid.spacing
    A = np.zeros(matrix.shape)
    I = np.eye(grid.num_points,grid.num_points)
    inv_diag = I * (1.0/np.diag(matrix))
    return (I - omega * inv_diag @ matrix, inv_diag)
  
M, N = build_jacobi_matrices(grid, discretization, omega=2/3)
M

def analytical_eigenvalue_jacobi(grid, k):
    h = grid.spacing
    c = k * np.pi
    return 1 - 2*2/3 * np.sin(c*h/2)**2#

analytical_eigenvals = np.array([analytical_eigenvalue_jacobi(grid,k) for k in range(1,grid.num_points+1)])
analytical_eigenvals

for k in range(0, grid.num_points):
    print(analytical_eigenvector(grid, k+1))


fig, ax = plt.subplots()
omegas = np.linspace(0, 1, 100)
opt_1 = np.abs(1 - omegas)
opt_2 = np.abs(1 - 2 * omegas)
ax.plot(omegas, opt_1, label="$| 1 - \omega | $", linestyle='dotted')
ax.plot(omegas, opt_2, label="$| 1 - 2\omega | $", linestyle='dotted')
ax.axvline(2/3, label="$\omega = 2/3$", c='black')

ax.plot(omegas, np.maximum(opt_1, opt_2),
        label="$\max(| 1 - \omega |, | 2 - \omega |)$",
        linestyle='solid', c='green')

ax.legend()

