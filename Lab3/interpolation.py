import json
import math
import matplotlib.pyplot as plt


def func(x):
    return math.acos(x) + x


def lagrangeValue(nodes, values, x):
    result = 0.0

    for i in range(len(nodes)):
        temp = values[i]

        for j in range(len(nodes)):
            if j != i:
                temp *= (x - nodes[j]) / (nodes[i] - nodes[j])
        
        result += temp

    return result


def divDiffs(nodes, values):
    n = len(nodes)
    table = [[0.0] * n for _ in range(n)]

    for i in range(n):
        table[i][0] = values[i]

    for k in range(1, n):
        for i in range(n - k):
            table[i][k] = ((table[i + 1][k - 1] - table[i][k - 1]) / (nodes[i + k] - nodes[i]))
    return table


def newtonValue(nodes, table, x):
    result = table[0][0]
    skobki = 1.0

    for k in range(1, len(nodes)):
        skobki *= x - nodes[k - 1]
        result += table[0][k] * skobki

    return result


def main():
    with open("1.json", encoding="utf-8") as file:
        data = json.load(file)
    xStar = data["xStar"]

    for nodeSet in data["nodeSets"]:
        nodes = nodeSet["nodes"]
        n = len(nodes)
        if n == 0 or len(set(nodes)) != n:
            raise ValueError("Ошибка json-файла: нужен непустой набор различных узлов.\n")

        values = [func(x) for x in nodes]
        table = divDiffs(nodes, values)

        print("\n >>>>>>>> Набор", nodeSet["name"], " <<<<<<<<")
        print("Исходные узлы и значения:")
        for i in range(n):
            print(f"x[{i}] = {nodes[i]:.6f}, y[{i}] = {values[i]:.12f}")

        print("\nМногочлен Лагранжа:")
        print("P_L(x) =")
        for i in range(n):
            temp = f"({values[i]:.12g})"
            for j in range(n):
                if j != i:
                    temp += f" * (x - ({nodes[j]})) / ({nodes[i] - nodes[j]:.12g})"
            if i > 0:
                print("+", end=" ")
            print(temp)

        print("\nТаблица разделённых разностей:")
        print("x_i\t\tf_0\t\tf_1\t\tf_2\t\tf_3")
        for i in range(n):
            print(f"{nodes[i]:10.6f}", end=" ")
            for k in range(n - i):
                print(f"{table[i][k]:16.10f}", end=" ")
            print()

        print("\nМногочлен Ньютона:")
        print("P_N(x) =")
        skobkiText = ""
        for k in range(n):
            if k > 0:
                skobkiText += f" * (x - ({nodes[k - 1]}))"
                print("+", end=" ")
            print(f"({table[0][k]:.12g})" + skobkiText)

        lagrange = lagrangeValue(nodes, values, xStar)
        newton = newtonValue(nodes, table, xStar)
        exact = func(xStar)

        print(f"\nx* = {xStar}")
        print(f"Лагранж: P_L(x*) = {lagrange:.15f}")
        print(f"Ньютон:  P_N(x*) = {newton:.15f}")
        print(f"Исходная функция: f(x*) = {exact:.15f}")
        print(f"Погрешность Лагранжа: {abs(exact - lagrange):.12g}")
        print(f"Погрешность Ньютона:  {abs(exact - newton):.12g}")
        print(f"Разность Лагранж - Ньютон: {lagrange - newton:.12g}")
        print(f"Модуль разности методов:  {abs(lagrange - newton):.12g}")

        # График
        left = min(min(nodes), xStar)
        right = max(max(nodes), xStar)
        plotX = [left + (right - left) * i / 200 for i in range(201)]
        funcY = [func(x) for x in plotX]
        lagrangeY = [lagrangeValue(nodes, values, x) for x in plotX]
        newtonY = [newtonValue(nodes, table, x) for x in plotX]
        plt.figure()
        plt.plot(plotX, funcY, label="f(x)")
        plt.plot(plotX, lagrangeY, label="Лагранж")
        plt.plot(plotX, newtonY, "--", label="Ньютон")
        plt.scatter(nodes, values, color="red", label="Узлы")
        plt.axvline(xStar, color="gray", linestyle=":", label="x*")
        plt.scatter([xStar], [lagrange], color="green", label="P(x*)")
        plt.title("Набор " + nodeSet["name"])
        plt.xlabel("x")
        plt.ylabel("y")
        plt.grid()
        plt.legend()
        plt.show()
        plt.close()


if __name__ == "__main__":
    main()
