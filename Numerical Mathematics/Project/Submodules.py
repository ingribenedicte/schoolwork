def Newton(f, df, x_0, err_tol):
    x = x_0
    error = abs(f(x))
    step_error_list = []
    while error > err_tol:
        x_next = x - f(x) / df(x)
        step_error = x_next - x
        step_error_list.append(step_error)
        x = x_next
        error = abs(f(x))
    return x, step_error_list, len(step_error_list)

def Improved_Newton(f, df, ddf, x_0, err_tol):
    x = x_0
    error = abs(f(x))
    step_error_list = []
    while error > err_tol:
        x_next = x - (f(x) * df(x)) / (df(x)** 2 - f(x) * ddf(x))
        step_error = x_next - x
        step_error_list.append(step_error)
        x = x_next
        error = abs(f(x))
    return x, step_error_list, len(step_error_list)

def Olver(f, df, ddf, x_0, err_tol):
    x = x_0
    error = abs(f(x))
    step_error_list = []
    while error > err_tol:
        x_next = x - f(x) / df(x) - (ddf(x) * (f(x)) ** 2) / (2 * (df(x)) ** 3)
        step_error = x_next - x
        step_error = x_next - x
        step_error_list.append(step_error)
        x = x_next
        error = abs(f(x))
    return x, step_error_list, len(step_error_list)

def Improved_Olver(f, df, ddf, x_0, err_tol):
    x = x_0
    error = abs(f(x))
    step_error_list = []
    while error > err_tol:
        x_next = x  # +noe
        step_error = x_next - x
        step_error_list.append(step_error)
        x = x_next
        error = abs(f(x))
    return x, step_error_list, len(step_error_list)