def outer_ring_expansion(cell, r, theta):
    assert cell.kappa is not None, "Runtime Error: Total mass of individual cell is not yet calculated!"
    ################# Task 1 begins #########################
    term1 = cell.kappa * np.log(r)

    ratio = cell.a / r
    num = 1 - (ratio**2)
    den = 1 - (2 * ratio * np.cos(theta - cell.quad_angles)) + (ratio**2)
    f_s = cell.outer_ring_expansion_terms - (cell.kappa * np.log(cell.a))
    term2 = np.sum(f_s * num / den) * cell.h / (2 * np.pi)

    potential =  term1 + term2
    ################## Task 1 begins ########################
    return potential
