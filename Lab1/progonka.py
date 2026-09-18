import sys

from toolbox import *

# That's it.
def progonka(lowerDiagonal, mainDiagonal, upperDiagonal, rightSide, epsilon=EpsilonConst):
    size = len(mainDiagonal)

    if size < 2:
        raise ValueError("Ошибка прогонки: для метода прогонки требуется хотя бы два уравнения.")

    if (len(lowerDiagonal) != size - 1 or len(upperDiagonal) != size - 1 or len(rightSide) != size):
        raise ValueError("Ошибка прогонки: размеры диагоналей и правой части не согласованы")

    alpha = [0.0 for i in range(size)]
    beta = [0.0 for i in range(size)]

    if abs(mainDiagonal[0]) < epsilon:
        raise ValueError("Ошибка прогонки: нулевой знаменатель в первой строке прогонки.")

    alpha[0] = -upperDiagonal[0] / mainDiagonal[0]
    beta[0] = rightSide[0] / mainDiagonal[0]

    for row in range(1, size - 1):
        temp = (mainDiagonal[row] + lowerDiagonal[row - 1] * alpha[row - 1])

        if abs(temp) < epsilon:
            raise ValueError(f"Ошибка прогонки: нулевой знаменатель на шаге {row + 1} прогонки.")

        alpha[row] = -upperDiagonal[row] / temp
        beta[row] = (rightSide[row] - lowerDiagonal[row - 1] * beta[row - 1]) / temp

    lastTemp = mainDiagonal[-1] + lowerDiagonal[-1] * alpha[-2]

    if abs(lastTemp) < epsilon:
        raise ValueError("Ошибка прогонки: нулевой знаменатель в последней строке прогонки.")

    solution = [0.0 for i in range(size)]
    solution[-1] = (rightSide[-1] - lowerDiagonal[-1] * beta[-2]) / lastTemp

    alpha[-1] = 0.0
    beta[-1] = solution[-1]

    for row in range(size - 2, -1, -1):
        solution[row] = alpha[row] * solution[row + 1] + beta[row]

    return solution, alpha, beta

def main():
    if len(sys.argv) != 2:
        print("Ошибка ввода: верный ввод: python progonka.py 2.json")
        return

    data = loadInputData(sys.argv[1])
    solution, alpha, beta = progonka(data["lowerDiagonal"], data["mainDiagonal"], data["upperDiagonal"], data["rightSide"])

    print("Прогоночные коэффициенты:")
    for index in range(len(alpha)):
        print(f"alpha[{index + 1}] = {formatNumber(alpha[index])}, " f"beta[{index + 1}] = {formatNumber(beta[index])}")
    printVector("Решение:", solution) 

if __name__ == "__main__":
    main()