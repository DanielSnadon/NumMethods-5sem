import math
import sys

from lu import lu, solve
from toolbox import *

def rebuild(matrix, rightSide):
    size = len(matrix)
    iterationMatrix = createMatrix(size, size)
    constVector = [0.0 for i in range(size)]

    for row in range(size):
        diagonal = matrix[row][row]
        if abs(diagonal) < EpsilonConst:
            raise ValueError("Ошибка построения итерационной формы: на диагонали есть ноль.")

        constVector[row] = rightSide[row] / diagonal

        for column in range(size):
            if row != column:
                iterationMatrix[row][column] = -matrix[row][column] / diagonal

    return iterationMatrix, constVector

# That's it.
def mpi(iterationMatrix, constVector, epsilon, maxIterations):
    size = len(constVector)
    previousSolution = [0.0 for i in range(size)]

    for iteration in range(1, maxIterations + 1):
        currentSolution = [0.0 for i in range(size)]

        for row in range(size):
            currentSolution[row] = constVector[row]
            for column in range(size):
                currentSolution[row] += iterationMatrix[row][column] * previousSolution[column]

        hasBadValue = False
        for value in currentSolution:
            if not math.isfinite(value):
                hasBadValue = True
                break
        if hasBadValue:
            return previousSolution, iteration, False

        difference = vectorInfNorm(subtractVectors(currentSolution, previousSolution))

        if difference < epsilon:
            return currentSolution, iteration, True

        previousSolution = currentSolution

    return previousSolution, maxIterations, False

# That's it.
def zeidelMethod(iterationMatrix, constVector, epsilon, maxIterations):
    size = len(constVector)
    previousSolution = [0.0 for _ in range(size)]

    for iteration in range(1, maxIterations + 1):
        currentSolution = previousSolution.copy()

        for row in range(size):
            currentSolution[row] = constVector[row]

            for column in range(row):
                currentSolution[row] += iterationMatrix[row][column] * currentSolution[column]

            for column in range(row + 1, size):
                currentSolution[row] += iterationMatrix[row][column] * previousSolution[column]

        hasBadValue = False
        for value in currentSolution:
            if not math.isfinite(value):
                hasBadValue = True
                break
        if hasBadValue:
            return previousSolution, iteration, False

        difference = vectorInfNorm(subtractVectors(currentSolution, previousSolution))
            
        if difference < epsilon:
            return currentSolution, iteration, True

        previousSolution = currentSolution

    return previousSolution, maxIterations, False

def main():
    if len(sys.argv) != 2:
        print("Ошибка ввода: верный ввод: python mpizeidel.py 2.json")
        return

    data = loadInputData(sys.argv[1])
    matrix = data["matrix"]
    rightSide = data["rightSide"]
    epsilon = data["epsilon"]
    maxIterations = data["maxIterations"]
    checkLinearSystem(matrix, rightSide)

    iterationMatrix, constVector = rebuild(matrix, rightSide)

    mpiSolution, mpiIterations, mpiConverged = mpi(iterationMatrix, constVector, epsilon, maxIterations)
    zeidelSolution, zeidelIterations, zeidelConverged = zeidelMethod(iterationMatrix, constVector, epsilon, maxIterations)

    lowerMatrix, upperMatrix, permutationMatrix, useless = lu(matrix)
    ref = solve(lowerMatrix, upperMatrix, permutationMatrix, rightSide)

    print(f"Epsilon: {epsilon}")
    printMatrix("Матрица B:", iterationMatrix)
    printVector("Вектор c:", constVector)
    printVector("Эталонное решение через LU:", ref)

    print("\nМетод простых итераций:")
    print(f"Сошёлся: {'да' if mpiConverged else 'нет'}")
    print(f"Количество итераций: {mpiIterations}")
    printVector("Решение:", mpiSolution)

    print("\nМетод Зейделя:")
    print(f"Сошёлся: {'да' if zeidelConverged else 'нет'}")
    print(f"Количество итераций: {zeidelIterations}")
    printVector("Решение:", zeidelSolution)
    zeidelResidual = subtractVectors(
        multiplyMatrixAndVector(matrix, zeidelSolution),
        rightSide,
    )
    
    if mpiConverged and zeidelConverged:
        if zeidelIterations < mpiIterations:
            fasterMethod = "метод Зейделя"
        elif mpiIterations < zeidelIterations:
            fasterMethod = "метод простых итераций"
        else:
            fasterMethod = "одинаково"
        print(f"\nБыстрее по количеству итераций: {fasterMethod}")

if __name__ == "__main__":
    main()
