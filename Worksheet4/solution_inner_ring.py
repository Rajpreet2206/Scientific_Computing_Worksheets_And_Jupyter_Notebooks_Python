def inner_ring_expansion(cell, r, theta):
    ################ Task 2 begins #######################
    ratio = r / cell.a
    num = 1 - (ratio**2)
    den = 1 - (2 * ratio * np.cos(theta - cell.quad_angles)) + (ratio**2)

    potential = (np.sum(cell.inner_ring_expansion_terms * num / den) * cell.h / (2 * np.pi))
    ################ Task 2 ends ########################
    return potential
