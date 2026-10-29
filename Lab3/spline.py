import json
import matplotlib.pyplot as plt


def progonka(lower, diagonal, upper, b):
    n = len(diagonal)
    p = [0.0] * n
    q = [0.0] * n
    p[0] = -upper[0] / diagonal[0]
    q[0] = b[0] / diagonal[0]
    
    for i in range(1, n):
        denominator = diagonal[i] + lower[i] * p[i - 1]
        if denominator == 0:
            raise ValueError("Ошибка метода прогонки: нулевой знаменатель.\n")
        p[i] = -upper[i] / denominator
        q[i] = (b[i] - lower[i] * q[i - 1]) / denominator

    c = [0.0] * n
    c[-1] = q[-1]
    for i in range(n - 2, -1, -1):
        c[i] = p[i] * c[i + 1] + q[i]

    return c


def buildSpline(nodes, values):
    n = len(nodes)
    if n < 2 or len(values) != n:
        raise ValueError("Ошибка: нужно хотя бы два узла и столько же значений.\n")

    # h[i] и m[i] -> заставляет задуматься.
    h = []
    m = []
    for i in range(n - 1):
        h.append(nodes[i + 1] - nodes[i])
        m.append((values[i + 1] - values[i]) / h[i])

    lower = [0.0] * n
    diagonal = [0.0] * n
    upper = [0.0] * n
    b = [0.0] * n
    diagonal[0] = 1.0
    diagonal[-1] = 1.0

    for i in range(1, n - 1):
        lower[i] = h[i - 1]
        diagonal[i] = 2 * (h[i - 1] + h[i])
        upper[i] = h[i]
        b[i] = 3 * (m[i] - m[i - 1])

    c = progonka(lower, diagonal, upper, b)

    a = []
    b = []
    d = []
    for i in range(n - 1):
        a.append(values[i])
        b.append(m[i] - h[i] * (2 * c[i] + c[i + 1]) / 3)
        d.append((c[i + 1] - c[i]) / (3 * h[i]))

    return a, b, c, d


def getValue(nodes, a, b, c, d, x):
    if x < nodes[0] or x > nodes[-1]:
        raise ValueError("Ошибка: точка должна находиться между крайними узлами.\n")
    
    interval = len(nodes) - 2
    for i in range(len(nodes) - 1):
        if x < nodes[i + 1]:
            interval = i
            break

    t = x - nodes[interval]
    result = a[interval] + b[interval]*t + c[interval]*t**2 + d[interval]*t**3
    return result, interval


def main():
    with open("2.json", encoding="utf-8") as file:
        data = json.load(file)

    nodes = data["nodes"]
    values = data["values"]
    xStar = data["xStar"]
    a, b, c, d = buildSpline(nodes, values)

    print("Исходная таблица:")
    for i in range(len(nodes)):
        print(f"x[{i}] = {nodes[i]}, y[{i}] = {values[i]}")

    print("\nКоэффициенты и многочлены участков:")
    for i in range(len(nodes) - 1):
        print(f"\nУчасток {i+1}: [{nodes[i]}; {nodes[i+1]}]")
        print(f"a = {a[i]:.12g}, b = {b[i]:.12g}, c = {c[i]:.12g}, d = {d[i]:.12g}")
        print(f"t = x - ({nodes[i]})")
        print(f"S(x) = ({a[i]:.12g}) + ({b[i]:.12g})*t + ({c[i]:.12g})*t^2 + ({d[i]:.12g})*t^3")

    result, interval = getValue(nodes, a, b, c, d, xStar)
    print(f"\nx* = {xStar}, S(x*) = {result:.15f}")
    print(f"Номер участка: {interval+1}, индекс Python: {interval}")
    print(f"Интервал: [{nodes[interval]}; {nodes[interval+1]}]")

    print("\nПроверка прохождения каждого участка через оба его конца:")
    for i in range(len(nodes) - 1):
        for j in (i, i + 1):
            t = nodes[j] - nodes[i]
            value = a[i] + b[i]*t + c[i]*t**2 + d[i]*t**3
            print(f"Участок {i+1}, x={nodes[j]}: S={value:.12f}, y={values[j]}, ошибка={abs(value-values[j]):.3e}")

    last = len(nodes) - 2
    secondLeft = 2 * c[0]
    secondRight = 2*c[last] + 6*d[last]*(nodes[-1] - nodes[-2])
    print("\nВторые производные на концах:")
    print(f"S''(x_0) = {secondLeft:.12g}")
    print(f"S''(x_n) = {secondRight:.12g}")

    # График
    for i in range(len(nodes) - 1):
        plotX = [nodes[i] + (nodes[i+1]-nodes[i])*j/100 for j in range(101)]
        plotY = []
        for x in plotX:
            t = x - nodes[i]
            plotY.append(a[i] + b[i]*t + c[i]*t**2 + d[i]*t**3)
        plt.plot(plotX, plotY, label=f"Участок {i+1}")
    plt.scatter(nodes, values, color="red", label="Узлы")
    plt.scatter([xStar], [result], color="black", label="(x*, S(x*))")
    plt.title("Кубический сплайн")
    plt.xlabel("x")
    plt.ylabel("S(x)")
    plt.grid()
    plt.legend()
    plt.show()
    plt.close()


if __name__ == "__main__":
    main()
