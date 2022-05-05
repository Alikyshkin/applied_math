from random import *
import numpy as np
import matplotlib.pyplot as plt
import min_search as ms

# Метод градиентного спуска
def descent_gradient(f, grad, x, eps, min_f, n=100):
    steps = 0
    points = [np.copy(x)]
    step = 1
    for _ in range(n):
        steps += 1
        gr = grad(x)
        g = lambda l_: f(x - gr * l_)
        l = min_f(g, eps, step)
        step = l
        x = x - gr * l
        points.append(np.copy(x))
        if np.linalg.norm(gr) < eps:
            break
    return x, steps, points


# Метод сопряженных градиентов
def conjugate_gradient(f, grad, x0, eps, min_f, n=100):
    k = 0
    xk = x0
    pk = -grad(x0)
    steps = 0
    points = [np.copy(x0)]
    step = 0.7
    while True:
        steps += 1
        alpha = min_f(lambda a: f(xk + a * pk), eps, step)
        xk1 = xk + alpha * pk
        step = alpha
        points.append(np.copy(xk1))
        if np.linalg.norm(grad(xk1)) < eps:
            return xk1, steps, points
        if k + 1 == n:
            k = 0
            x0 = xk
        else:
            beta = (np.linalg.norm(grad(xk1)) ** 2) / (np.linalg.norm(grad(xk)) ** 2)
            pk1 = -grad(xk1) + beta * pk
            pk = pk1
            k += 1
        xk = xk1


# Для графика поле
def make_field(f):
    x = np.linspace(-15, 15, 100)
    y = np.linspace(-15, 15, 100)

    X, Y = np.meshgrid(x, y)  # сетка
    Z = f([X, Y])
    return X, Y, Z


# Рисуем график
def draw_level_lines(func, points: list):
    x, y, z = make_field(func)
    points_x, points_y = [], []
    for p in points:
        points_x.append(p[0])
        points_y.append(p[1])
    fig, ax = plt.subplots()
    ax.contour(x, y, z)
    ax.scatter(points_x, points_y,
               c=[random() for _ in range(len(points_x))])
    ax.plot(points_x, points_y, c='red')
    plt.show()


def do_report(func, t):
    ans, steps, points = t
    print(f"{ans} steps = {steps}")
    draw_level_lines(func, points)


def test_on(func, grad, x0, eps, min_f):
    print(f"Все ответы для x0 = {x0}, eps = {eps}")
    print(f"Метод сопряженных градиентов = ", end="")
    do_report(func, conjugate_gradient(func, grad, np.copy(x0), eps, min_f))
    print(f"Метод градиентного спуска = ", end="")
    do_report(func, descent_gradient(func, grad, np.copy(x0), eps, min_f))


def main():
    # f = x^3 * sin(x)
    # f = x^2
    # функции которые мы придумали
    # x[0] = x, x[1] = y
    def func(x):
        return 4 * x[0] ** 2 + 2 * x[1] ** 2 + 3 * x[0] * x[1] + 2 * x[0] + 2 * x[1]
        # return x[0] ** 2 + x[1] ** 2 + x[0] * x[1] + 10 * x[0] - 10 * x[1] + 10
        # return 3 * x[0] ** 2 + x[1] ** 2 - x[0] * x[1] - 4 * x[0]

    def grad(x):
        return np.array([8 * x[0] + 3 * x[1] + 2, 3 * x[0] + 4 * x[1] + 2])
        # return np.array([2 * x[0] + x[1] + 10, x[0] + 2 * (x[1] - 5)])
        # return np.array([6 * x[0] - x[1] - 4, 2 * x[1] - x[0]])

    print(f"\n==Метод с постоянным шагом==\n")
    test_on(func, grad, np.array([15.0, 10.0]), 0.001, ms.const)
    print(f"\n==Метод дробления==\n")
    test_on(func, grad, np.array([15.0, 10.0]), 0.001, ms.step_decreaser)
    print(f"\n==Метод линейного поиска==\n")
    test_on(func, grad, np.array([15.0, 10.0]), 0.001, ms.linear_search)
    print(f"\n==Метод золотого сечения==\n")
    test_on(func, grad, np.array([15.0, 10.0]), 0.1, ms.gold_section)
    print(f"\n==Метод Фибоначчи==\n")
    test_on(func, grad, np.array([15.0, 10.0]), 0.001, ms.fibonacci)


main()
