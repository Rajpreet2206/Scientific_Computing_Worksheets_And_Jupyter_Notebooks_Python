def explicit_euler(A, init_condition, time_steps=1, dt=0.1):
    output = []
    I = np.eye(A.shape[0])
    output.append(init_condition)
    # TODO: implement the explicit euler method
    # In each time step, the new solution for the position should be appended to output.
    ##########################TODO BEGINS##########################
    for i in range(time_steps):
        solution = np.matmul((I + A*dt ), output[-1])
        output.append(solution)
    ##########################TODO ENDS##########################
    return output
