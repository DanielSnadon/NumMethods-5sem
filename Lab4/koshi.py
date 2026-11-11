import argparse
import json
import math
from pathlib import Path
import matplotlib.pyplot as plt

xStart, xEnd = 2.0, 3.0
yStart = math.sqrt(2)
zStart = 1.5 * math.sqrt(2)


def acceleration(x, y, z): # z'
    return (0.75 * y - 0.5 * z) / (x * (x - 1))


def perfect(x):
    if x == xStart:
        return yStart
    
    return -2 * x**1.5 + math.sqrt(2) * (2*x + 1) * math.sqrt(x - 1)


def eulerStep(x, y, z, h):
    return y + h*z, z + h*acceleration(x, y, z)


def rkStep(x, y, z, h):
    # 1
    dy1 = z
    dz1 = acceleration(x, y, z)

    # 2
    yMid = y + h/2 * dy1
    zMid = z + h/2 * dz1
    dy2 = zMid
    dz2 = acceleration(x + h/2, yMid, zMid)

    # 3
    yMid = y + h/2 * dy2
    zMid = z + h/2 * dz2
    dy3 = zMid
    dz3 = acceleration(x + h/2, yMid, zMid)

    # 4
    yEnd = y + h * dy3
    zEnd = z + h * dz3
    dy4 = zEnd
    dz4 = acceleration(x + h, yEnd, zEnd)

    nextY = y + h/6 * (dy1 + 2*dy2 + 2*dy3 + dy4)
    nextZ = z + h/6 * (dz1 + 2*dz2 + 2*dz3 + dz4)
    
    return nextY, nextZ


def oneStep(nodes, h, step):
    y = [yStart]
    z = [zStart]

    for i in range(len(nodes) - 1):
        nextY, nextZ = step(nodes[i], y[i], z[i], h)
        y.append(nextY)
        z.append(nextZ)

    return y, z


def adamsWeight(values, i):
    return (55*values[i] - 59*values[i-1] + 37*values[i-2] - 9*values[i-3]) / 24


def solveAdams(nodes, h):
    y, z = oneStep(nodes[:4], h, rkStep)
    dz = [acceleration(nodes[i], y[i], z[i]) for i in range(4)]

    for i in range(3, len(nodes) - 1):
        nextY = y[i] + h * adamsWeight(z, i)
        nextZ = z[i] + h * adamsWeight(dz, i)
        y.append(nextY)
        z.append(nextZ)
        dz.append(acceleration(nodes[i + 1], nextY, nextZ))

    return y, z


def makeGrid(h):
    if not math.isfinite(h) or h <= 0:
        raise ValueError("Ошибка makeGrid: шаг h должен быть конечным положительным числом.")
    
    n = round((xEnd - xStart) / h)

    if not math.isclose(n*h, xEnd - xStart, rel_tol=1e-12, abs_tol=1e-12):
        raise ValueError("Ошибка makeGrid: шаг h должен разбивать [2; 3] на целое число участков.")
    
    if n < 4:
        raise ValueError("Ошибка makeGrid: нужно минимум четыре шага: три для старта и один для Адамса.")
    
    return [xStart + i*h for i in range(n)] + [xEnd]


def solve(nodes, h):
    return [
        ("Эйлер", 1, oneStep(nodes, h, eulerStep)),
        ("Рунге–Кутта", 4, oneStep(nodes, h, rkStep)),
        ("Адамс", 4, solveAdams(nodes, h)),
    ]


def printTable(nodes, badGrid, goodGrid):
    for (name, order, (y, z)), (_, _, (goodY, _)) in zip(badGrid, goodGrid):
        print(f"\n{name}, порядок p = {order}")
        print(f"{'x':>6} {'y(h)':>13} {'z(h)':>13} {'точное y':>13} "f"{'y(h/2)':>13} {'ошибка h':>12} {'ошибка h/2':>12} {'РР для h/2':>12}")

        badGridErrors, goodGridErrors, delta = [], [], []
        for i, x in enumerate(nodes):
            exact = perfect(x)
            goodGridValue = goodY[2*i]
            badGridErrors.append(abs(y[i] - exact))
            goodGridErrors.append(abs(goodGridValue - exact))
            delta.append(abs(goodGridValue - y[i]) / (2**order - 1))

            print(f"{x:6.3f} {y[i]:13.9f} {z[i]:13.9f} {exact:13.9f} "f"{goodGridValue:13.9f} {badGridErrors[-1]:12.3e} "f"{goodGridErrors[-1]:12.3e} {delta[-1]:12.3e}")

        print(f"Максимумы на общих узлах: ошибка h = {max(badGridErrors):.3e}; "f"ошибка h/2 = {max(goodGridErrors):.3e}; оценка РР = {max(delta):.3e}")


def drawGraphs(nodes, solutions, h, folder, show):
    plotX = [xStart + (xEnd - xStart)*i/200 for i in range(201)]
    plt.figure(figsize=(8, 5))
    plt.plot(plotX, [perfect(x) for x in plotX], "k--", label="Точное решение")
    for (name, _, (y, _)), marker in zip(solutions, ["o", "s", "^"]):
        plt.plot(nodes, y, marker=marker, markersize=4, label=name)
    plt.title(f"4.1 - Решение задачи Коши, h = {h:g}")
    plt.xlabel("x")
    plt.ylabel("y(x)")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(folder / "koshi_solution.png", dpi=150)

    plt.figure(figsize=(8, 5))
    for (name, _, (y, _)), marker in zip(solutions, ["o", "s", "^"]):
        points = [(x, abs(value - perfect(x))) for x, value in zip(nodes, y)]
        points = [(x, error) for x, error in points if error > 0]
        plt.semilogy([p[0] for p in points], [p[1] for p in points],
                     marker=marker, markersize=4, label=name)
    plt.title(f"4.1 - Абсолютная ошибка, h = {h:g}")
    plt.xlabel("x")
    plt.ylabel("|y_численное − y_точное|")
    plt.grid(alpha=0.3, which="both")
    plt.legend()
    plt.tight_layout()
    plt.savefig(folder / "koshi_error.png", dpi=150)
    if show:
        plt.show()
    plt.close("all")


def main():
    parser = argparse.ArgumentParser(description="ЛР 4.1, вариант 20")
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    folder = Path(__file__).resolve().parent
    with (folder / "1.json").open(encoding="utf-8") as file:
        h = float(json.load(file)["h"])

    nodes, goodNodes = makeGrid(h), makeGrid(h/2)
    badGrid, goodGrid = solve(nodes, h), solve(goodNodes, h/2)

    print(f"Вариант 20: h = {h:g}, h/2 = {h/2:g}")
    print("Начальная задача: y(2) = sqrt(2), y'(2) = 3*sqrt(2)/2.")
    
    printTable(nodes, badGrid, goodGrid)
    drawGraphs(nodes, badGrid, h, folder, args.show)

    print("\nГрафики сохранены: koshi_solution.png, koshi_error.png.")


if __name__ == "__main__":
    main()
