import json

EpsilonConst = 1e-12

def loadInputData(filePath):
    with open(filePath, "r", encoding="utf-8") as inputFile:
        return json.load(inputFile)

def copyMatrix(matrix):
    return [row.copy() for row in matrix]

def createMatrix(rowCount, columnCount, value=0.0):
    return [[value for i in range(columnCount)] for j in range(rowCount)]

def createOneMatrix(size):
    matrix = createMatrix(size, size)

    for index in range(size):
        matrix[index][index] = 1.0

    return matrix

def transpMatrix(matrix):
    rowCount = len(matrix)
    columnCount = len(matrix[0])
    result = createMatrix(columnCount, rowCount)

    for row in range(rowCount):
        for column in range(columnCount):
            result[column][row] = matrix[row][column]

    return result

def multiplyMatrix(firstMatrix, secondMatrix):
    firstRowCount = len(firstMatrix)
    firstColumnCount = len(firstMatrix[0])
    secondRowCount = len(secondMatrix)
    secondColumnCount = len(secondMatrix[0])

    if firstColumnCount != secondRowCount:
        raise ValueError("Ошибка умножения: невозможно умножить матрицы таких размеров.")

    result = createMatrix(firstRowCount, secondColumnCount)

    for row in range(firstRowCount):
        for column in range(secondColumnCount):
            for i in range(firstColumnCount):
                result[row][column] += (firstMatrix[row][i] * secondMatrix[i][column])

    return result

def multiplyMatrixAndVector(matrix, vector):
    if len(matrix[0]) != len(vector):
        raise ValueError("Ошибка умножения: размеры матрицы и вектора не согласованы.")

    result = [0.0 for i in range(len(matrix))]

    for row in range(len(matrix)):
        for column in range(len(vector)):
            result[row] += matrix[row][column] * vector[column]

    return result

def subtractVectors(firstVector, secondVector):
    if len(firstVector) != len(secondVector):
        raise ValueError("Ошибка вычитания: размеры векторов не совпадают.")

    return [firstVector[index] - secondVector[index] for index in range(len(firstVector))]

def vectorInfNorm(vector):
    return max(abs(value) for value in vector)

def matrixInfNorm(matrix):
    maxRowSum = 0.0

    for row in matrix:
        rowSum = 0.0
        for value in row:
            rowSum += abs(value)
        if rowSum > maxRowSum:
            maxRowSum = rowSum

    return maxRowSum

def matrixMaxDiff(firstMatrix, secondMatrix):
    if len(firstMatrix) != len(secondMatrix):
        raise ValueError("Ошибка сравнения: размеры матриц не совпадают.")

    maxDifference = 0.0

    for row in range(len(firstMatrix)):
        if len(firstMatrix[row]) != len(secondMatrix[row]):
            raise ValueError("Ошибка сравнения: размеры матриц не совпадают.")

        for column in range(len(firstMatrix[row])):
            difference = abs(firstMatrix[row][column] - secondMatrix[row][column])
            maxDifference = max(maxDifference, difference)

    return maxDifference

def checkSquareMatrix(matrix):
    if not matrix:
        raise ValueError("Ошибка квадратности: матрица должна быть квадратной.")

    size = len(matrix)

    for row in matrix:
        if len(row) != size:
            raise ValueError("Ошибка квадратности: матрица должна быть квадратной.")

def checkLinearSystem(matrix, rightSide):
    checkSquareMatrix(matrix)
    if len(matrix) != len(rightSide):
        raise ValueError("Ошибка системности: размер матрицы не совпадает с размером правой части")

def isSymmetric(matrix, epsilon=EpsilonConst):
    checkSquareMatrix(matrix)

    for row in range(len(matrix)):
        for column in range(row + 1, len(matrix)):
            if abs(matrix[row][column] - matrix[column][row]) > epsilon:
                return False

    return True

def formatNumber(value, digits=6):
    if isinstance(value, complex):
        realPart = round(value.real, digits)
        imaginaryPart = round(value.imag, digits)

        sign = "+"
        absImaginaryPart = imaginaryPart
        if imaginaryPart < 0:
            sign = "-"
            absImaginaryPart = -imaginaryPart

        realText = f"{realPart:.{digits}f}"
        imaginaryText = f"{absImaginaryPart:.{digits}f}"

        return realText + " " + sign + " " + imaginaryText + "i"

    if abs(value) < 10 ** (-digits):
        value = 0.0
    return f"{value:.{digits}f}"

def printMatrix(title, matrix, digits=6):
    print(title)
    for row in matrix:
        print("  [" + "  ".join(formatNumber(value, digits) for value in row) + "]")

def printVector(title, vector, digits=6):
    print(title + " [" + ", ".join(formatNumber(value, digits) for value in vector) + "]")