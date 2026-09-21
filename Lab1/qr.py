import math
import sys

import numpy as np
from toolbox import *

# That's it?..
def qr(matrix, epsilon=EpsilonConst):
    checkSquareMatrix(matrix)
    size = len(matrix)
    vectorMatrix = createOneMatrix(size)
    upperMatrix = copyMatrix(matrix)

    for column in range(size - 1):
        for row in range(size - 1, column, -1):
            firstValue = upperMatrix[column][column]
            secondValue = upperMatrix[row][column]

            if abs(secondValue) < epsilon:
                continue

            length = math.hypot(firstValue, secondValue) # ~sqrt(.^2 + .^2)
            cosugla = firstValue / length
            sinugla = secondValue / length

            for currentColumn in range(column, size):
                oldTopValue = upperMatrix[column][currentColumn]
                oldBottomValue = upperMatrix[row][currentColumn]
                upperMatrix[column][currentColumn] = cosugla * oldTopValue + sinugla * oldBottomValue
                upperMatrix[row][currentColumn] = -sinugla * oldTopValue + cosugla * oldBottomValue
            upperMatrix[row][column] = 0.0

            for currentRow in range(size):
                oldLeftValue = vectorMatrix[currentRow][column]
                oldRightValue = vectorMatrix[currentRow][row]
                vectorMatrix[currentRow][column] = cosugla * oldLeftValue + sinugla * oldRightValue
                vectorMatrix[currentRow][row] = -sinugla * oldLeftValue + cosugla * oldRightValue

    return vectorMatrix, upperMatrix

# 2x2 is enough :p
def quasiCheck(matrix, epsilon):
    size = len(matrix)

    for row in range(size):
        for column in range(row - 1):
            if abs(matrix[row][column]) >= epsilon:
                return False

    prevWasLarge = False
    for row in range(1, size):
        currentLargeSubDiagonal = abs(matrix[row][row - 1]) >= epsilon
        if prevWasLarge and currentLargeSubDiagonal:
            return False
        prevWasLarge = currentLargeSubDiagonal

    return True


def SobstQuasiZnach(matrix, epsilon):
    sobstZnach = []
    index = 0
    size = len(matrix)

    while index < size:
        if index < size - 1 and abs(matrix[index + 1][index]) >= epsilon:
            firstDiagonal = matrix[index][index]
            upperValue = matrix[index][index + 1]
            lowerValue = matrix[index + 1][index]
            secondDiagonal = matrix[index + 1][index + 1]

            # l^2 - (a+d)l + (ad - bc) = 0
            t1 = firstDiagonal + secondDiagonal
            t2 = firstDiagonal * secondDiagonal - upperValue * lowerValue
            discriminant = t1 * t1 - 4.0 * t2

            if discriminant >= 0:
                root = math.sqrt(discriminant)
                sobstZnach.append((t1 + root) / 2.0)
                sobstZnach.append((t1 - root) / 2.0)
            else:
                root = math.sqrt(-discriminant)
                sobstZnach.append(complex(t1 / 2.0, root / 2.0))
                sobstZnach.append(complex(t1 / 2.0, -root / 2.0))

            index += 2
        else:
            sobstZnach.append(matrix[index][index])
            index += 1

    return sobstZnach


def qrWrap(matrix, epsilon, maxIterations):
    checkSquareMatrix(matrix)
    currentMatrix = copyMatrix(matrix)

    for iteration in range(maxIterations + 1):
        if quasiCheck(currentMatrix, epsilon):
            sobstZnach = SobstQuasiZnach(currentMatrix, epsilon)
            return sobstZnach, iteration, currentMatrix

        vectorMatrix, upperMatrix = qr(currentMatrix)
        currentMatrix = multiplyMatrix(upperMatrix, vectorMatrix)

    raise ValueError("Ошибка ожиданий: QR-алгоритм не сошёлся за заданное число итераций.")

def useNumPy(matrix):
    return np.linalg.eigvals(np.array(matrix, dtype=float))

def main():
    if len(sys.argv) != 2:
        print("Использование: python qr.py 5.json")
        return

    data = loadInputData(sys.argv[1])

    originalMatrix = data["matrix"]
    epsilon = data["epsilon"]
    maxIterations = data["maxIterations"]

    vectorMatrix, upperMatrix = qr(originalMatrix)
    temp = multiplyMatrix(vectorMatrix, upperMatrix)

    sobstZnach, iterations, finalMatrix = qrWrap(originalMatrix, epsilon, maxIterations)

    printMatrix("Исходная матрица A:", originalMatrix)
    printMatrix("Матрица Q:", vectorMatrix)
    printMatrix("Матрица R:", upperMatrix)
    printMatrix("Произведение Q * R:", temp)

    print(f"Максимальная ошибка Q * R = A: " f"{matrixMaxDiff(temp, originalMatrix):.3e}")

    print(f"Количество QR-итераций: {iterations}")
    printMatrix("Итоговая матрица:", finalMatrix)
    printVector("Собственные значения, найденные QR-алгоритмом:", sobstZnach)

    hasComplexValues = any(isinstance(value, complex) for value in sobstZnach)

    if hasComplexValues:
        numpySobstZnach = useNumPy(originalMatrix)
        printVector("Проверка комплексных собственных значений:", list(numpySobstZnach))

if __name__ == "__main__":
    main()