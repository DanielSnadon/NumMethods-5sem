import math
import sys

from toolbox import *

def upperValue(matrix):
    value = 0.0

    for row in range(len(matrix)):
        for column in range(row + 1, len(matrix)):
            value += matrix[row][column] ** 2

    return math.sqrt(value)

def largestOffDiagonal(matrix):
    size = len(matrix)
    pivotRow = 0
    pivotColumn = 1

    for row in range(size):
        for column in range(row + 1, size):
            if abs(matrix[row][column]) > abs(matrix[pivotRow][pivotColumn]):
                pivotRow = row
                pivotColumn = column

    return pivotRow, pivotColumn

# That's it.
def jacobi(matrix, epsilon, maxIterations):
    checkSquareMatrix(matrix)
    if len(matrix) == 1:
        return [matrix[0][0]], [[1.0]], 0, copyMatrix(matrix)
    if not isSymmetric(matrix):
        raise ValueError("Ошибка метода вращений Якоби: метод требует симметричную матрицу.")

    size = len(matrix)
    currentMatrix = copyMatrix(matrix)
    vectorMatrix = createOneMatrix(size)
    iteration = 0

    while upperValue(currentMatrix) >= epsilon:
        if iteration >= maxIterations:
            raise ValueError("Ошибка метода вращений Якоби: метод Якоби не сошёлся за заданное число итераций")

        pivotRow, pivotColumn = largestOffDiagonal(currentMatrix)
        angle = 0.5 * math.atan2(2.0 * currentMatrix[pivotRow][pivotColumn], currentMatrix[pivotRow][pivotRow] - currentMatrix[pivotColumn][pivotColumn])
        cosugla = math.cos(angle)
        sinugla = math.sin(angle)

        oldFirstDiagonal = currentMatrix[pivotRow][pivotRow]
        oldSecondDiagonal = currentMatrix[pivotColumn][pivotColumn]
        oldPivot = currentMatrix[pivotRow][pivotColumn]

        for index in range(size):
            if index == pivotRow or index == pivotColumn:
                continue

            oldFirstValue = currentMatrix[index][pivotRow]
            oldSecondValue = currentMatrix[index][pivotColumn]
            newFirstValue = cosugla * oldFirstValue + sinugla * oldSecondValue
            newSecondValue = -sinugla * oldFirstValue + cosugla * oldSecondValue

            currentMatrix[index][pivotRow] = newFirstValue
            currentMatrix[pivotRow][index] = newFirstValue
            currentMatrix[index][pivotColumn] = newSecondValue
            currentMatrix[pivotColumn][index] = newSecondValue

        currentMatrix[pivotRow][pivotRow] = cosugla * cosugla * oldFirstDiagonal + 2.0 * cosugla * sinugla * oldPivot + sinugla * sinugla * oldSecondDiagonal
        currentMatrix[pivotColumn][pivotColumn] = sinugla * sinugla * oldFirstDiagonal - 2.0 * cosugla * sinugla * oldPivot + cosugla * cosugla * oldSecondDiagonal

        currentMatrix[pivotRow][pivotColumn] = 0.0
        currentMatrix[pivotColumn][pivotRow] = 0.0

        for row in range(size):
            oldFirstValue = vectorMatrix[row][pivotRow]
            oldSecondValue = vectorMatrix[row][pivotColumn]
            vectorMatrix[row][pivotRow] = cosugla * oldFirstValue + sinugla * oldSecondValue
            vectorMatrix[row][pivotColumn] = -sinugla * oldFirstValue + cosugla * oldSecondValue
        iteration += 1

    sobstZnach = [currentMatrix[index][index] for index in range(size)]
    return sobstZnach, vectorMatrix, iteration, currentMatrix

def main():
    if len(sys.argv) != 2:
        print("Использование: python jacobi.py 4.json")
        return

    data = loadInputData(sys.argv[1])
    originalMatrix = data["matrix"]
    epsilon = data["epsilon"]

    sobstZnach, vectorMatrix, iterations, diagonalMatrix = jacobi(originalMatrix, epsilon, data["maxIterations"])

    lambdaMatrix = createMatrix(len(sobstZnach), len(sobstZnach))
    for index in range(len(sobstZnach)):
        lambdaMatrix[index][index] = sobstZnach[index]

    leftSide = multiplyMatrix(originalMatrix, vectorMatrix)
    rightSide = multiplyMatrix(vectorMatrix, lambdaMatrix)

    print(f"Epsilon: {epsilon}")
    print(f"Количество вращений: {iterations}")
    printVector("Собственные значения:", sobstZnach)
    printMatrix("Матрица собственных векторов V:", vectorMatrix)

    for column in range(len(sobstZnach)):
        eigenvector = [
            vectorMatrix[row][column]
            for row in range(len(sobstZnach))
        ]
        printVector(f"Собственный вектор v{column + 1}:", eigenvector)

    printMatrix("Почти диагональная матрица:", diagonalMatrix)
    printMatrix("Произведение A * V:", leftSide)
    printMatrix("Произведение V * Lambda:", rightSide)
    print(f"Максимальная ошибка A * V = V * Lambda: {matrixMaxDiff(leftSide, rightSide):.3e}")

if __name__ == "__main__":
    main()
