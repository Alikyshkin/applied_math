from random import *
from numdifftools import Gradient
import numpy as np
import matplotlib.pyplot as plt
import min_search as ms


# метод наискорейшего спуска
def steepest_descent_gradient(f, grad, x, eps, min_f, n=100):
    steps = 0
    points = [np.copy(x)]
    for _ in range(n):
        steps += 1
        gr = grad(x)
        g = lambda l_: f(x - gr * l_)
        l = min_f(g, -10, 10, eps)
        x = x - gr * l
        points.append(np.copy(x))
        if np.linalg.norm(gr) < eps:
            break
    return x, steps, points


# метод градиентного спуска
def descent_gradient(grad, v, k, eps, n=100):
    steps = 0
    points = [np.copy(v)]
    for i in range(1, n + 1):
        steps += 1
        l = 1 / min(i + 3, k + 3)
        v = v - grad(v) * l
        points.append(np.copy(v))
        if np.linalg.norm(grad(v)) < eps:
            break
    return v, steps, points


# метод ньютона
def newton(grad, hesse, v, eps, n=100):
    steps = 0
    points = [np.copy(v)]
    for _ in range(n):
        steps += 1
        g = -1 * grad(v)
        h = hesse(v)
        s = np.linalg.solve(h, g)
        v += s
        points.append(np.copy(v))
        if np.linalg.norm(s) <= eps:
            return v, steps, points


# метод сопряженных градиентов
def conjugate_gradient(f, grad, x0, eps, min_f, n=100):
    k = 0
    xk = x0
    pk = -grad(x0)
    steps = 0
    points = [np.copy(x0)]
    while True:
        steps += 1
        alpha = min_f(lambda a: f(xk + a * pk), 0, 1, eps)
        xk1 = xk + alpha * pk
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


# метод сопряженных направлений
def conugate_directions(f, grad, x, eps, min_f):  # Пауэлл
    n = len(x)
    D = np.eye(n)
    steps = 0
    points = [np.copy(x)]
    while True:
        steps += 1
        r = min_f(lambda l: f(x + l * D[:, n - 1]), -10, 10, eps)
        x = x + r * D[:, n - 1]
        points.append(np.copy(x))
        y = x
        for i in range(n):
            r = min_f(lambda l: f(x + l * D[:, i]), -10, 10, eps)
            x = x + r * D[:, i]
        for i in range(n - 1):
            D[:, i] = D[:, i + 1]
        D[:, n - 1] = y - x

        if np.linalg.norm(D[:, n - 1]) < eps:
            return x, steps, points


# для графика поле
def make_field(f):
    x = np.linspace(-15, 15, 100)
    y = np.linspace(-15, 15, 100)

    X, Y = np.meshgrid(x, y)  # сетка
    Z = f([X, Y])
    return X, Y, Z


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


def generate_matrix(n, k):
    matrix = [[0] * n for _ in range(n)]
    for i in range(n):
        matrix[i][i] = randint(1, k)
    i = randint(0, n - 1)
    matrix[i][i] = k
    j = i
    while j == i:
        j = randint(0, n - 1)
    matrix[j][j] = 1
    return matrix


def generate_function(n, k):
    matrix = generate_matrix(n, k)
    b = [randint(0, 100) for _ in range(n)]
    func = lambda x: sum(x[i] ** 2 * matrix[i][i] - b[i] * x[i] for i in range(n))
    return func


# траектория реализованных методов (линия шагов оптимизации)
def n_by_k(n):
    fig, ax = plt.subplots()
    s = [0] * 10000
    for j in range(2, 5):
        func = generate_function(j * 10, n)
        grad = Gradient(func)
        start = np.array([randint(-10000, 10000) for _ in range(j * 10)])
        ans, steps, points = steepest_descent_gradient(func, grad, start, 0.001, ms.fibonacci, 10000)
        s[j] = steps
        print(f"{j * 10} = {steps}")

    y = [0]
    x = [0]
    ax.plot(x, y, color='red')
    for i in range(2, 4):
        y = [s[i], s[i + 1]]
        x = [i * 10, ((i + 1) * 10)]
        ax.plot(x, y, color='red')
    plt.show()


def test_on(func, grad, hesse, x0, eps, min_f, ans):
    print(f"test all functions on x0 = {x0} eps = {eps}\n  ans = {ans}")
    print(f"steepest descent method = ", end="")
    do_report(func, steepest_descent_gradient(func, grad, np.copy(x0), eps, min_f))
    print(f"conjugate_gradient = ", end="")
    do_report(func, conjugate_gradient(func, grad, np.copy(x0), eps, min_f))
    print(f"pavel_method = ", end="")
    do_report(func, conugate_directions(func, grad, np.copy(x0), eps, min_f))
    print(f"descent_gradient = ", end="")
    do_report(func, descent_gradient(grad, np.copy(x0), 20, eps))
    print(f"newton method = ", end="")
    do_report(func, newton(grad, hesse, np.copy(x0), eps))
    n_by_k(10)


def main():
    # f = x^3 * sin(x)
    # f = x^2
    # функции которые мы придумали
    def func(x):
        return 4 * x[0] ** 2 + 2 * x[1] ** 2 + 3 * x[0] * x[1] + 2 * x[0] + 2 * x[1]
        return x[0] ** 2 + x[1] ** 2
        return (x[0] + 1) ** 2 + x[1] ** 2
        return 6 * x[0] ** 2 - 4 * x[0] * x[1] + 3 * x[1] ** 2 + 4 * sqrt(5) * (x[0] + 2 * x[1]) + 22
        return (x[0]) ** 2 + x[1] ** 2 + x[0] * x[1]
        return (x[0] - 3) ** 2
        return x ** 3 * sin(x)

    # d/dx(x^3 sin(x)) = x^2 (3 sin(x) + x cos(x))
    def grad(x):
        return np.array([8 * x[0] + 3 * x[1] + 2, 3 * x[0] + 4 * x[1] + 2])
        return np.array([2 * x[0], 2 * x[1]])
        return np.array([12 * x[0] - 4 * x[1] + 4 * sqrt(5), -4 * x[0] + 6 * x[1] - 8 * sqrt(5)])
        return np.array([2 * x[0] + x[1], 2 * x[1] + x[0]])
        return x ** 2 * 3 * sin(x) + x ** 3 * cos(x)

    # d^2/dx^2(x^3 sin(x)) = x (6 x cos(x) - (x^2 - 6) sin(x))
    def hesse(x):
        return np.array([[8, 3], [3, 4]])
        return np.array([[2, 0], [0, 2]])
        return np.array([[2, 1], [1, 2]])
        return np.array([[2]])
        return np.array([[12, -4], [-4, 6]])  # only one variable

    test_on(func, grad, hesse, np.array([15.0, 10.0]), 0.001, ms.fibonacci, [-0.08695585, -0.43478307])


if __name__ == '__main__':
    main()
