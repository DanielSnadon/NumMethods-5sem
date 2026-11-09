import json
import math
from pathlib import Path


def function(x):
    return math.sqrt(x) / (4 + 3*x)


def pervoobraznaya(x): # Аналитическое
    t = math.sqrt(x)
    return 2*t/3 - 4/(3*math.sqrt(3)) * math.atan(math.sqrt(3)*t/2)


def integrate(f, a, b, h):
    n = round((b - a) / h)
    if n < 1 or not math.isclose(n*h, b - a, rel_tol=1e-12):
        raise ValueError("Ошибка integrate: шаг должен разбивать отрезок на целое число участков.")
    if n % 2 != 0:
        raise ValueError("Ошибка integrate: для Симпсона число участков должно быть чётным.")

    nodes = [a + i*h for i in range(n)] + [b]
    values = [f(x) for x in nodes]

    left = 0.0
    right = 0.0
    middle = 0.0
    trapezoid = 0.0

    for i in range(n):
        left += h * values[i]
        right += h * values[i + 1]
        middle += h * f((nodes[i] + nodes[i + 1]) / 2)
        trapezoid += h * (values[i] + values[i + 1]) / 2

    simpsonSum = values[0] + values[n]
    for i in range(1, n):
        koeff = 4 if i % 2 == 1 else 2
        simpsonSum += koeff * values[i]
    simpson = h * simpsonSum / 3

    return left, right, middle, trapezoid, simpson


def rungeRomberg(h1, h2, ratio, p):
    correction = (h2 - h1) / (ratio**p - 1)
    ipp = h2 + correction
    return ipp, abs(correction)


def main():
    with Path(__file__).with_name("5.json").open(encoding="utf-8") as file:
        data = json.load(file)

    a, b = data["a"], data["b"]
    h1, h2 = data["h1"], data["h2"]
    if not (0 <= a < b and 0 < h2 < h1):
        raise ValueError("Необходимо, чтобы 0 <= a < b и 0 < h2 < h1.")

    fh1 = integrate(function, a, b, h1)
    fh2 = integrate(function, a, b, h2)
    ratio = h1 / h2
    exact = pervoobraznaya(b) - pervoobraznaya(a)

    methods = [
        ("Левые прямоугольники", 1),
        ("Правые прямоугольники", 1),
        ("Средние прямоугольники", 2),
        ("Трапеции", 2),
        ("Симпсон", 4),
    ]

    print(f"Интеграл sqrt(x)/(4+3*x) на [{a:g}; {b:g}]")
    print(f"h1 = {h1:g}, h2 = {h2:g}, h1/h2 = {ratio:g}")
    print(f"Точное значение F(b)-F(a) ≈ {exact:.12f}")

    print("\nРезультаты и уточнение Рунге–Ромберга:")
    print(f"{'Метод':24} {'I(h1)':>15} {'I(h2)':>15} {'p':>2} {'I(РР)':>15} {'Оценка для h2':>15}")
    refinedValues = []
    for (name, p), h1Value, h2Value in zip(methods, fh1, fh2):
        refined, estimate = rungeRomberg(h1Value, h2Value, ratio, p)
        refinedValues.append(refined)
        print(f"{name:24} {h1Value:15.12f} {h2Value:15.12f} {p:2d} {refined:15.12f} {estimate:15.6e}")

    print("\nАбсолютные погрешности относительно F(b)-F(a):")
    print(f"{'Метод':24} {'Для h1':>15} {'Для h2':>15} {'После РР':>15}")
    for i, (name, _) in enumerate(methods):
        print(f"{name:24} {abs(fh1[i]-exact):15.6e} "
              f"{abs(fh2[i]-exact):15.6e} {abs(refinedValues[i]-exact):15.6e}")

    for label, results in [("для h1", fh1), ("для h2", fh2), ("после РР", refinedValues)]:
        best = min(range(len(results)), key=lambda i: abs(results[i] - exact))
        print(f"Лучшая точность {label}: {methods[best][0]}.")


if __name__ == "__main__":
    main()
