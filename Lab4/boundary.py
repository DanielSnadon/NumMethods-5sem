import argparse
import json
import math
from pathlib import Path
import matplotlib.pyplot as plt


def exactY(x):
    return 2*x + 1 + math.exp(2*x)


def makeGrid(h):
    if not math.isfinite(h) or h <= 0:
        raise ValueError("Ошибка makeGrid: шаг должен быть конечным положительным числом.")
    
    n = round(1 / h)

    if n < 2 or not math.isclose(n*h, 1, rel_tol=1e-12):
        raise ValueError("Ошибка makeGrid: шаг должен делить [0; 1] на целое число частей, не меньше двух.")
    
    return [i*h for i in range(n)] + [1.0]


def rkStep(y, z, w, h):
    dy1, dz1, dw1 = z, w, 2*w

    # 2
    zMid = z + h/2 * dz1
    wMid = w + h/2 * dw1
    dy2, dz2, dw2 = zMid, wMid, 2*wMid

    # 3
    zMid = z + h/2 * dz2
    wMid = w + h/2 * dw2
    dy3, dz3, dw3 = zMid, wMid, 2*wMid

    # 4
    zEnd = z + h * dz3
    wEnd = w + h * dw3
    dy4, dz4, dw4 = zEnd, wEnd, 2*wEnd

    nextY = y + h/6 * (dy1 + 2*dy2 + 2*dy3 + dy4)
    nextZ = z + h/6 * (dz1 + 2*dz2 + 2*dz3 + dz4)
    nextW = w + h/6 * (dw1 + 2*dw2 + 2*dw3 + dw4)

    return nextY, nextZ, nextW


def shoot(nodes, h, s):
    y, z, w = 2.0, 4.0, s
    values = [y]
    for _ in range(len(nodes) - 1):
        y, z, w = rkStep(y, z, w, h)
        values.append(y)

    promax = z - 2*y + 4

    return values, promax


def solveShooting(nodes, h):
    s0, s1 = 0.0, 8.0

    _, promax0 = shoot(nodes, h, s0)
    _, promax1 = shoot(nodes, h, s1)

    difference = promax1 - promax0

    if abs(difference) < 1e-12:
        raise ValueError("Ошибка solveShooting: две невязки совпали: секущую построить нельзя.")
    
    s = s1 - promax1 * (s1 - s0) / difference

    values, promax = shoot(nodes, h, s)

    return values, [(s0, promax0), (s1, promax1), (s, promax)]


def solveLinear(matrix, rightSide): # Гаусс
    n = len(rightSide)
    for k in range(n):
        pivot = max(range(k, n), key=lambda i: abs(matrix[i][k]))

        if matrix[pivot][k] == 0:
            raise ValueError("Ошибка solveLinear: разностная система вырождена для этой сетки.")
        
        matrix[k], matrix[pivot] = matrix[pivot], matrix[k]
        rightSide[k], rightSide[pivot] = rightSide[pivot], rightSide[k]

        for i in range(k + 1, n):
            if matrix[i][k] == 0:
                continue

            factor = matrix[i][k] / matrix[k][k]
            matrix[i][k] = 0.0
            for j in range(k + 1, n):
                matrix[i][j] -= factor * matrix[k][j]
            rightSide[i] -= factor * rightSide[k]

    values = [0.0] * n
    for i in range(n - 1, -1, -1):
        known = sum(matrix[i][j] * values[j] for j in range(i + 1, n))
        values[i] = (rightSide[i] - known) / matrix[i][i]

    return values


def solveDifferences(nodes, h):
    n = len(nodes) - 1
    matrix = [[0.0] * (n + 1) for _ in nodes]
    rightSide = [0.0] * (n + 1)

    # first
    matrix[0][:3] = [-3.0, 4.0, -1.0]
    rightSide[0] = 8*h

    # mid
    for i in range(1, n):
        x = nodes[i]
        halfSlope = (2*x + 1)*h/2
        matrix[i][i-1] = x + halfSlope
        matrix[i][i] = -2*x + 2*h*h
        matrix[i][i+1] = x - halfSlope

    # last
    matrix[n][-3:] = [1.0, -4.0, 3 - 4*h]
    rightSide[n] = -8*h

    return solveLinear(matrix, rightSide)


def printTable(nodes, name, order, badGrid, goodGrid):
    print(f"\n{name}, порядок p = {order}")
    print(f"{'x':>8} {'y(h)':>13} {'точное y':>13} {'y(h/2)':>13} "f"{'ошибка h':>12} {'ошибка h/2':>12} {'РР для h/2':>12}")

    badGridErrors, goodGridErrors, estimates = [], [], []
    for i, x in enumerate(nodes):
        exact = exactY(x)
        goodGridValue = goodGrid[2*i]
        badGridErrors.append(abs(badGrid[i] - exact))
        goodGridErrors.append(abs(goodGridValue - exact))
        estimates.append(abs(goodGridValue - badGrid[i]) / (2**order - 1))
        print(f"{x:8.4f} {badGrid[i]:13.9f} {exact:13.9f} {goodGridValue:13.9f} "f"{badGridErrors[-1]:12.3e} {goodGridErrors[-1]:12.3e} {estimates[-1]:12.3e}")
    
    print(f"Максимумы на общих узлах: ошибка h = {max(badGridErrors):.3e}; "f"ошибка h/2 = {max(goodGridErrors):.3e}; оценка РР = {max(estimates):.3e}")


def drawGraphs(nodes, solutions, h, folder, show):
    plotX = [i/200 for i in range(201)]
    plt.figure(figsize=(8, 5))
    plt.plot(plotX, [exactY(x) for x in plotX], "k--", label="Точное решение")
    for (name, values), marker in zip(solutions, ["o", "s"]):
        plt.plot(nodes, values, marker=marker, markersize=4, label=name)
    plt.title(f"4.2 · Решение краевой задачи, h = {h:g}")
    plt.xlabel("x")
    plt.ylabel("y(x)")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(folder / "boundary_solution.png", dpi=150)
    plt.figure(figsize=(8, 5))
    for (name, values), marker in zip(solutions, ["o", "s"]):
        points = [(x, abs(y - exactY(x))) for x, y in zip(nodes, values)]
        points = [(x, e) for x, e in points if e > 0]
        plt.semilogy([x for x, _ in points], [e for _, e in points],
                     marker=marker, markersize=4, label=name)
    plt.title(f"4.2 · Абсолютная ошибка, h = {h:g}")
    plt.xlabel("x")
    plt.ylabel("|y численное − y точное|")
    plt.grid(alpha=0.3, which="both")
    plt.legend()
    plt.tight_layout()
    plt.savefig(folder / "boundary_error.png", dpi=150)
    if show:
        plt.show()
    plt.close("all")


def main():
    parser = argparse.ArgumentParser(description="ЛР 4.2, вариант 20")
    parser.add_argument("--show", action="store_true", help="Также открыть окна графиков")
    args = parser.parse_args()

    folder = Path(__file__).resolve().parent
    with (folder / "2.json").open(encoding="utf-8") as file:
        h = float(json.load(file)["h"])
    nodes, goodGridNodes = makeGrid(h), makeGrid(h/2)
    shooting, trials = solveShooting(nodes, h)
    goodGridShooting, goodGridTrials = solveShooting(goodGridNodes, h/2)
    differences = solveDifferences(nodes, h)
    goodGridDifferences = solveDifferences(goodGridNodes, h/2)

    print(f"Вариант 20: h = {h:g}, h/2 = {h/2:g}")
    print("Точное решение: y(x) = 2*x + 1 + exp(2*x).")

    for step, shots in [(h, trials), (h/2, goodGridTrials)]:
        print(f"\nСтрельба: h = {step:g}, s = y''(0)")
        for i, (s, promax) in enumerate(shots, 1):
            print(f"  Выстрел {i}: s = {s:.10g}, R(s) = {promax:+.3e}")

    printTable(nodes, "Метод стрельбы (RK4)", 4, shooting, goodGridShooting)
    printTable(nodes, "Конечные разности", 2, differences, goodGridDifferences)
    drawGraphs(nodes, [("Стрельба (RK4)", shooting), ("Конечные разности", differences)], h, folder, args.show)
    
    print("\nГрафики сохранены: boundary_solution.png, boundary_error.png.")


if __name__ == "__main__":
    main()
