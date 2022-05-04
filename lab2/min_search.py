import numpy as np
from scipy.misc import derivative


# def const():
#     return 0.3


def fibonacci(f, x_start, x_end, eps):
    fibonacciSeries = [1, 1]
    function_calcs = 0

    while 1:
        newItem = fibonacciSeries[len(fibonacciSeries) - 1] + fibonacciSeries[len(fibonacciSeries) - 2]
        fibonacciSeries.append(newItem)
        if newItem > abs(x_end - x_start) / eps:
            break

    n = len(fibonacciSeries) - 1
    x1 = x_start + (fibonacciSeries[n - 2] / fibonacciSeries[n]) * (x_end - x_start)
    x2 = x_start + (fibonacciSeries[n - 1] / fibonacciSeries[n]) * (x_end - x_start)
    fx1 = f(x1)
    fx2 = f(x2)
    function_calcs += 2

    iters = 0
    while n > 2:
        print("Iter:", iters, ", Curr interval: (", x_start, ", ", x_end, "), diff between prev interval: ", (abs(x_end - x_start)), ", xmin: ", (x_start + x_end)/2.0, ", func calculations: ", function_calcs)
        n -= 1
        iters += 1

        if fx1 < fx2:
            x_end = x2
            x2 = x1
            fx2 = fx1
            x1 = x_start + ((fibonacciSeries[n - 2] / fibonacciSeries[n]) * (x_end - x_start))
            fx1 = f(x1)
            function_calcs += 1

        else:
            x_start = x1
            x1 = x2
            fx1 = fx2
            x2 = x_start + (fibonacciSeries[n - 1] / fibonacciSeries[n]) * (x_end - x_start)
            fx2 = f(x2)
            function_calcs += 1

    return (x_start + x_end) / 2.0


def linear_search(f, eps):
    alpha = 0
    der = derivative(f, alpha)
    der_test = derivative(f, alpha + eps / 2)
    if (der >= 0 and der - der_test > 0) or (der < 0 and der - der_test < 0):
        while der * derivative(f, alpha) > 0:
            alpha += eps / 2
    else:
        while der * derivative(f, alpha) > 0:
            alpha -= eps / 2

    return alpha


def step_drill(f, step):
    v = 0.4
    der = 0
    i = 0
    while f(der) <= f(step):
        step = step * v
        i += 1
    return step


def gold_section(f, eps):
    a = 0
    b = 1
    phi = (1 + np.sqrt(5)) / 2
    rphi = 2 - phi
    x1 = a + rphi * (b - a)
    x2 = b - rphi * (b - a)
    f1 = f(x1)
    f2 = f(x2)
    while np.fabs(a - b) > eps:
        if f1 < f2:
            b = x2
            x2 = x1
            f2 = f1
            x1 = a + rphi * (b - a)
            f1 = f(x1)
        else:
            a = x1
            x1 = x2
            f1 = f2
            x2 = b - rphi * (b - a)
            f2 = f(x2)
    return x1 + x2 / 2
