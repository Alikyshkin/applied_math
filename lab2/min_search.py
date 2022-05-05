import numpy as np
from scipy.misc import derivative


# Постоянная величина шага
def const(f, eps, step):
    return 0.1

def step_decreaser(f, e, step):
    v = 0.4
    der = 0
    i = 0
    while f(der) <= f(step):
        step = step * v
        i += 1
    return step


# Метод золотого сечения
def gold_section(f, eps, step):
    a = 0
    b = 10
    phi = (1 + np.sqrt(5)) / 2
    rphi = 2 - phi
    x1 = a + rphi * abs(b - a)
    x2 = b - rphi * abs(b - a)
    f1 = f(x1)
    f2 = f(x2)
    while abs(a - b) > eps:
        if f1 < f2:
            b = x2
            x2 = x1
            f2 = f1
            x1 = a + rphi * abs(b - a)
            f1 = f(x1)
        else:
            a = x1
            x1 = x2
            f1 = f2
            x2 = b - rphi * abs(b - a)
            f2 = f(x2)
    step = x1 + x2 / 2
    return step


def fib(n):
    fib1 = fib2 = 1
    n = int(n) - 2
    while n > 0:
        fib1, fib2 = fib2, fib1 + fib2
        n -= 1
    return fib2


def fibonacci(f, eps, step):
    a = 0
    b = 10
    n = 0
    while fib(n) < (b - a) / (2 * eps):
        n += 1
    l = a + (fib(n - 2) / fib(n)) * (b - a)
    m = a + (fib(n - 1) / fib(n)) * (b - a)
    yL = f(l)
    yM = f(m)
    for k in range(1, n - 1):
        if yL > yM:
            a = l
            l = m
            m = a + (fib(n - k - 1) / fib(n - k)) * (b - a)
            yM = f(m)
            yL = f(l)
        else:
            b = m
            m = l
            l = a + (fib(n - k - 2) / fib(n - k)) * (b - a)
            yL = f(l)
            yM = f(m)
    m = l + eps
    if f(l) <= f(m):
        a = l
    else:
        b = m
    return a + b / 2

def linear_search(f, e, step):
    alpha = 0
    der = derivative(f, alpha)
    der_test = derivative(f, alpha + e / 2)
    if (der >= 0 and der - der_test > 0) or (der < 0 and der - der_test < 0):
        while der * derivative(f, alpha) > 0:
            alpha += e / 2
    else:
        while der * derivative(f, alpha) > 0:
            alpha -= e / 2

    return alpha
