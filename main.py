import math
from math import *

point1 = -4
point2 = 2
epsil = 1e-6


def f(x):
    return exp(sin(x)) * (x ** 2)


def sign(x):
    if (x > 0):
        return 1
    if (x < 0):
        return -1
    if (x == 0):
        return 0


# * Метод дихотомии =============

def DichotomyMethod(x_start, x_end, eps):
    delta = eps / 2
    interval_length = x_end - x_start
    iters = 0
    function_calcs = 0

    while interval_length > eps:
        x1 = (x_start + x_end - delta) / 2
        x2 = (x_start + x_end + delta) / 2

        if f(x1) > f(x2):
            x_start = x1

        else:
            x_end = x2

        interval_length = x_end - x_start
        function_calcs += 2
        iters += 1
        # print("Iter:", iters, ", Curr interval: (", x_start, ", ", x_end, "), diff between prev interval: ", abs(x_end - x_start), ", xmin: ", (x_start + x_end)/2.0, ", func calculations: ", function_calcs)
    return (x_start + x_end) / 2.0, f((x_start + x_end) / 2.0), function_calcs, iters


result = DichotomyMethod(point1, point2, epsil)
print('===Метод дихотомии===\n\n', result, "\n")


# * Золотого сечения =============

def GoldenSectionMethod(x_start, x_end, eps):
    interval_length = x_end - x_start

    x1 = x_start + (3 - math.sqrt(5)) / 2 * (x_end - x_start)
    x2 = x_start + (math.sqrt(5) - 1) / 2 * (x_end - x_start)
    y1 = f(x1)
    y2 = f(x2)
    iters = 0
    function_calcs = 2

    while interval_length > eps:
        if y1 > y2:
            x_start = x1
            x1 = x2
            x2 = x_start + (math.sqrt(5) - 1) / 2 * (x_end - x_start)
            y1 = y2
            y2 = f(x2)
            function_calcs += 1

        else:
            x_end = x2
            x2 = x1
            x1 = x_start + (3 - math.sqrt(5)) / 2 * (x_end - x_start)
            y2 = y1
            y1 = f(x1)
            function_calcs += 1

        interval_length = x_end - x_start
        iters += 1
        # print("Iter:", iters, ", Curr interval: (", x_start, ", ", x_end, "), diff between prev interval: ", abs(x_end - x_start), ", xmin: ", (x_start + x_end)/2.0, ", func calculations: ", function_calcs)

    return (x_start + x_end) / 2.0, f((x_start + x_end) / 2.0), function_calcs, iters


result = GoldenSectionMethod(point1, point2, epsil)
print('===Метод золотого сечения===\n\n', result, "\n")


# * Фиббоначи =============

def FibonacciMethod(x_start, x_end, eps):
    fibonacciSeries = [1, 1]
    function_calcs = 0

    while (1):
        newItem = fibonacciSeries[len(fibonacciSeries) - 1] + fibonacciSeries[len(fibonacciSeries) - 2]
        fibonacciSeries.append(newItem)
        if (newItem > abs(x_end - x_start) / eps):
            break

    n = len(fibonacciSeries) - 1
    x1 = x_start + (fibonacciSeries[n - 2] / fibonacciSeries[n]) * (x_end - x_start)
    x2 = x_start + (fibonacciSeries[n - 1] / fibonacciSeries[n]) * (x_end - x_start)
    fx1 = f(x1)
    fx2 = f(x2)
    function_calcs += 2

    iters = 0
    while (n > 2):
        # print("Iter:", iters, ", Curr interval: (", x_start, ", ", x_end, "), diff between prev interval: ", (abs(x_end - x_start)), ", xmin: ", (x_start + x_end)/2.0, ", func calculations: ", function_calcs)
        n -= 1
        iters += 1

        if (fx1 < fx2):
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

    return (x_start + x_end) / 2.0, f((x_start + x_end) / 2.0), function_calcs, iters


result = FibonacciMethod(point1, point2, epsil)
print('===Метод фиббоначи===\n\n', result, "\n")


# * Комбинированный метод Брента

def parabolaMin(x1, x2, x3, y1, y2, y3):
    return x2 - (pow(x2 - x1, 2) * (y2 - y3) - pow(x2 - x3, 2) * (y2 - y1)) / (
                2.0 * ((x2 - x1) * (y2 - y3) - (x2 - x3) * (y2 - y1)))


def isSame(val1, val2, val3, eps):
    return abs(val1 - val2) < eps and abs(val1 - val3) < eps and abs(val2 - val3) < eps;


def brent(x_start, x_end, eps):
    x = w = v = float((x_start + x_end) / 2)
    fx = fw = fv = float(f(x))
    function_calcs = 1
    d = float(x_end - x_start)
    e = float(x_end - x_start)
    K = float((3 - math.sqrt(5)) / 2)
    u = 0
    iters = 1

    while (abs(d) > eps):
        g = e
        e = d

        if (not (isSame(x, w, v, eps) and isSame(fx, fw, fv, eps))):
            u = parabolaMin(x_start, w, x_end, fx, fw, fv)

        if (x_start + eps <= u and x_end - eps >= u and abs(u - x) < 0.5 * g):
            d = abs(u - x)
        else:
            if (x < (x_start + x_end) / 2):
                u = x + K * (x_end - x)
                d = x_end - x
            else:
                u = x - K * (x - x_start)
                d = x - x_start
        if (abs(u - x) < eps):
            u = x + sign(u - x) * eps

        fu = f(u)
        function_calcs += 1

        if (fu <= fx):
            if (u >= x):
                x_start = x
            else:
                x_end = x

            v = w
            w = x
            x = u
            fv = fw
            fw = fx
            fx = fu

        else:
            if (u >= x):
                x_end = u

            else:
                x_start = u

            if (fu <= fw or w == x):
                w = u
                fv = fw
                fw = fu

            else:
                if (fu <= fu or v == x or v == w):
                    v = u
                    fv = fu

        iters += 1
    # print("Iter:", iters, ", Curr interval: (", x_start, ", ", x_end, "), diff between prev interval: ", (abs(x_end - x_start)), ", xmin: ", (x_start + x_end)/2.0, ", func calculations: ", function_calcs)

    return (x_start + x_end) / 2.0, f((x_start + x_end) / 2.0), function_calcs, iters


result = brent(point1, point2, epsil)
print('===Комбинированный метод Брента===\n\n', result, "\n")


# * Парабол =============

def ParabolaMethod(x_start, x_end, eps):
    iters = 0
    function_calcs = 0
    prevMin = 0
    middle = x_end - x_start
    f_middle = f(middle)

    if (x_start == 0):
        x_start += eps

    if (x_end == 0):
        x_end -= eps

    f_start = f(x_start)
    function_calcs += 1
    f_end = f(x_end)
    function_calcs += 1

    curMin = middle - (
                pow(middle - x_start, 2) * (f_middle - f_end) - pow(middle - x_end, 2) * (f_middle - f_start)) / (
                         2 * ((middle - x_start) * (f_middle - f_end) - (middle - x_end) * (f_middle - f_start)))
    f_curMin = f(curMin)
    function_calcs += 1

    while abs(curMin - prevMin) > eps:
        prevMin = curMin

        if (x_start < curMin and curMin < middle):
            if (f_curMin >= f_middle):
                x_start = curMin
                f_start = f_curMin
            else:
                if (f_curMin < f_middle):
                    x_end = middle
                    f_end = f_middle
                    middle = curMin
                    f_middle = f_curMin
        else:
            if middle < curMin and curMin < x_end:
                if f_middle >= f_curMin:
                    x_start = middle
                    f_start = f_middle
                    middle = curMin
                    f_middle = f_curMin
            else:
                if f_middle < f_curMin:
                    x_end = curMin
                    f_end = f_curMin

        iters += 1
        curMin = middle - (
                    pow(middle - x_start, 2) * (f_middle - f_end) - pow(middle - x_end, 2) * (f_middle - f_start)) / (
                             2 * ((middle - x_start) * (f_middle - f_end) - (middle - x_end) * (f_middle - f_start)))
        f_curMin = f(curMin)
        function_calcs += 1
    # print("x_start: ", x_start,  "x_end: ", x_end)

    # print("func calculations: ", function_calcs, ", iter: ", iters, "\ncurmin: ", curMin)

    return (x_start + x_end) / 2.0, f((x_start + x_end) / 2.0), curMin, function_calcs, iters


result = ParabolaMethod(point1, point2, epsil)
print('===Метод парабол===\n\n', result, "\n")