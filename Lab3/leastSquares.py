import json
import matplotlib.pyplot as plt


def solveLU(matrix, rightSide):
    n = len(matrix)
    lu = [row.copy() for row in matrix]
    permutedRight = rightSide.copy()

    for k in range(n):
        pivot = k
        for i in range(k + 1, n):
            if abs(lu[i][k]) > abs(lu[pivot][k]):
                pivot = i

        if lu[pivot][k] == 0:
            raise ValueError("Ошибка: система вырождена: коэффициенты нельзя найти однозначно.\n")

        lu[k], lu[pivot] = lu[pivot], lu[k]
        permutedRight[k], permutedRight[pivot] = permutedRight[pivot], permutedRight[k]

        for i in range(k + 1, n):
            lu[i][k] /= lu[k][k]
            for j in range(k + 1, n):
                lu[i][j] -= lu[i][k] * lu[k][j]

    z = [0.0] * n
    for i in range(n):
        z[i] = permutedRight[i]
        for j in range(i):
            z[i] -= lu[i][j] * z[j]

    coefficients = [0.0] * n
    for i in range(n - 1, -1, -1):
        coefficients[i] = z[i]
        for j in range(i + 1, n):
            coefficients[i] -= lu[i][j] * coefficients[j]
        coefficients[i] /= lu[i][i]

    return coefficients


def leastSquares(nodes, values, degree):
    size = degree + 1
    if len(nodes) != len(values) or len(set(nodes)) < size:
        raise ValueError("Ошибка: нужно одинаковое число x и y и не менее degree+1 различных x.\n")

    matrix = [[0.0] * size for _ in range(size)]
    rightSide = [0.0] * size

    for i in range(size):
        for j in range(size):
            for k in range(len(nodes)):
                matrix[i][j] += nodes[k] ** (i + j)

        for k in range(len(nodes)):
            rightSide[i] += values[k] * nodes[k] ** i

    coefficients = solveLU(matrix, rightSide)
    return coefficients, matrix, rightSide


def polyValue(coefficients, x):
    result = 0.0
    for i in range(len(coefficients)):
        result += coefficients[i] * x**i
    return result


def main():
    with open("3.json", encoding="utf-8") as file:
        data = json.load(file)

    nodes = data["nodes"]
    values = data["values"]
    if len(nodes) != len(values) or len(set(nodes)) < 4:
        raise ValueError("Ошибка: для степеней 1, 2, 3 нужно хотя бы 4 различных x и столько же y.\n")

    print("Исходные точки:")
    for i in range(len(nodes)):
        print(f"x[{i}] = {nodes[i]}, y[{i}] = {values[i]}")

    allCoefficients = []
    errors = []
    for degree in (1, 2, 3):
        coefficients, matrix, rightSide = leastSquares(nodes, values, degree)
        allCoefficients.append(coefficients)

        print(f"\nСтепень {degree}. Нормальная система: A | b")
        for i in range(degree + 1):
            for value in matrix[i]:
                print(f"{value:12.6f}", end=" ")
            print(f"| {rightSide[i]:12.6f}")

        print("Коэффициенты:")
        for i in range(len(coefficients)):
            print(f"a[{i}] = {coefficients[i]:.12g}")

        text = f"({coefficients[0]:.12g})"
        for i in range(1, len(coefficients)):
            text += f" + ({coefficients[i]:.12g})*x^{i}"
        print(f"P_{degree}(x) = " + text)

        sum = 0.0
        for i in range(len(nodes)):
            residual = polyValue(coefficients, nodes[i]) - values[i]
            sum += residual**2
        mse = sum / len(nodes)
        errors.append(mse)
        print(f"Сумма квадратов ошибок sum = {sum:.12g}")
        print(f"MSE{degree} = {mse:.12g}")

    print("\nСравнение ошибок на исходных точках:")
    for i in range(3):
        print(f"Степень {i+1}: MSE = {errors[i]:.12g}")
    best = errors.index(min(errors)) + 1
    print(f"Наименьшая MSE у многочлена степени {best}.")

    # График
    left = min(nodes)
    right = max(nodes)
    plotX = [left + (right-left)*i/200 for i in range(201)]
    plt.scatter(nodes, values, color="black", label="Исходные точки")
    for i in range(3):
        plotY = [polyValue(allCoefficients[i], x) for x in plotX]
        plt.plot(plotX, plotY, label=f"Степень {i+1}")
    plt.title("Метод наименьших квадратов")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid()
    plt.legend()
    plt.show()
    plt.close()


if __name__ == "__main__":
    main()
