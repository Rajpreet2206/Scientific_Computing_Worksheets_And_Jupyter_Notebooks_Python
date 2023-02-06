def stormer_verlet(init_condition, time_steps=10, dt=.1, g=9.81, r=1):
    output = []
    output.append(init_condition)
    # TODO: implement the Velocity Störmer Verlet method
    # In each time step, the new solution for the position and velocity should be appended to output.    ##########################TODO BEGINS##########################
    for i in range(time_steps):
        force_t = -1.0 * g * output[-1][0] 
        solution_pos = output[-1][0] + (dt * output[-1][1] / r) + (0.5 * dt**2 * force_t)
        force_t_dt = -1.0 * g * solution_pos
        solution_vel = output[-1][1] + ((0.5 * dt) * (force_t + force_t_dt))
        output.append(np.array([solution_pos, solution_vel]))
    ##########################TODO ENDS##########################
    return output
