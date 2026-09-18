import sys
from toolbox import *

# That's it.
def lu(matrix, epsilon=EpsilonConst):
    size = len(matrix)
    lowerMatrix = createOneMatrix(size)
    upperMatrix = copyMatrix(matrix)
    permutationMatrix = createOneMatrix(size)
    swapCount = 0

    for column in range(size):
        pivotRow = column

        for row in range(column + 1, size):
            if abs(upperMatrix[row][column]) > abs(upperMatrix[pivotRow][column]):
                pivotRow = row

        if abs(upperMatrix[pivotRow][column]) < epsilon:
            raise ValueError("Ошибка lu: матрица вырождена, ибо не найден ненулевой опорный элемент.")

        if pivotRow != column:
            upperMatrix[column], upperMatrix[pivotRow] = upperMatrix[pivotRow], upperMatrix[column]
            permutationMatrix[column], permutationMatrix[pivotRow] = permutationMatrix[pivotRow], permutationMatrix[column]

            for previousColumn in range(column):
                lowerMatrix[column][previousColumn], lowerMatrix[pivotRow][previousColumn] = (lowerMatrix[pivotRow][previousColumn], lowerMatrix[column][previousColumn])

                swapCount += 1
        
        for row in range(column + 1, size):
            multiplier = upperMatrix[row][column] / upperMatrix[column][column]
            lowerMatrix[row][column] = multiplier
            upperMatrix[row][column] = 0.0

            for currentColumn in range(column + 1, size):
                upperMatrix[row][currentColumn] -= multiplier * upperMatrix[column][currentColumn]
    
    return lowerMatrix, upperMatrix, permutationMatrix, swapCount

def forwardPodst(lowerMatrix, rightSide):
    size = len(lowerMatrix)
    result = [0.0 for i in range(size)]

    for row in range(size):
        knownSum = 0.0
        for column in range(row):
            knownSum += lowerMatrix[row][column] * result[column]

        result[row] = (rightSide[row] - knownSum) / lowerMatrix[row][row]

    return result

def backwardPodst(upperMatrix, rightSide):
    size = len(upperMatrix)
    result = [0.0 for i in range(size)]

    for row in range(size - 1, -1, -1):
        knownSum = 0.0
        for column in range(row + 1, size):
            knownSum += upperMatrix[row][column] * result[column]

        if abs(upperMatrix[row][row]) < EpsilonConst:
            raise ValueError("Ошибка обратного хода: нулевой диагональный элемент при обратной подстановке.")

        result[row] = (rightSide[row] - knownSum) / upperMatrix[row][row]

    return result

def solve(lowerMatrix, upperMatrix, permutationMatrix, rightSide):
    permutedRightSide = multiplyMatrixAndVector(permutationMatrix, rightSide)
    temp = forwardPodst(lowerMatrix, permutedRightSide)
    return backwardPodst(upperMatrix, temp)

def findInverse(lowerMatrix, upperMatrix, permutationMatrix):
    size = len(lowerMatrix)
    inverseMatrix = createMatrix(size, size)

    for column in range(size):
        unoVector = [0.0 for i in range(size)]
        unoVector[column] = 1.0
        solution = solve(lowerMatrix, upperMatrix, permutationMatrix, unoVector)

        for row in range(size):
            inverseMatrix[row][column] = solution[row]

    return inverseMatrix

def findDeterminant(upperMatrix, swapCount):
    if swapCount % 2 == 0:
        determinant = 1.0
    else:
        determinant = -1.0

    for index in range(len(upperMatrix)):
        determinant *= upperMatrix[index][index]

    return determinant

def main():
    if len(sys.argv) != 2:
        print("Ошибка ввода: верный ввод: python lu.py 1.json")
        return

    data = loadInputData(sys.argv[1])
    matrix = data["matrix"]
    rightSide = data["rightSide"]
    checkLinearSystem(matrix, rightSide)

    lowerMatrix, upperMatrix, permutationMatrix, swapCount = lu(matrix)
    solution = solve(lowerMatrix, upperMatrix, permutationMatrix, rightSide)
    inverseMatrix = findInverse(lowerMatrix, upperMatrix, permutationMatrix)
    determinant = findDeterminant(upperMatrix, swapCount)

    MultResult = multiplyMatrix(lowerMatrix, upperMatrix)
    permutedMatrix = multiplyMatrix(permutationMatrix, matrix)
    InverseChecked = multiplyMatrix(matrix, inverseMatrix)

    printMatrix("Матрица L:", lowerMatrix)
    printMatrix("Матрица U:", upperMatrix)
    printMatrix("Произведение L * U:", MultResult)
    printMatrix("Произведение P * A:", permutedMatrix)
    printVector("Решение Ax = b:", solution)
    printMatrix("Обратная матрица A^(-1):", inverseMatrix)
    print(f"Определитель A: {formatNumber(determinant)}")
    printMatrix("Произведение A * A^(-1):", InverseChecked)

if __name__ == "__main__":
    main()