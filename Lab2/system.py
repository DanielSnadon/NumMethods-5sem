import math
import sys

from toolbox import *

def functions(point, a):
    x1, x2 = point

    if x1 <= 0 or x2 <= 0:
        raise ValueError("Ошибка functions: некорректные параметры.")
    
    return [x1 * x1 - a * math.log10(x2) - 1, x1 * x1 - a * x1 * x2 + a]

def simpleStep(point, a):
    x1, x2 = point
    result = 1 + a * math.log10(x2)

    if result <= 0:
        raise ValueError("Ошибка simpleStep: итерация вышла из области положительных решений.")
    
    return [math.sqrt(result), (x1 * x1 + a) / (a * x1)]

def jacobian(point, a):
    x1, x2 = point

    return [[2 * x1, -a / (x2 * math.log(10))], [2 * x1 - a * x2, -a * x1]]

def newtonStep(point, a):
    values = functions(point, a)
    angle = jacobian(point, a)
    first, second = angle
    determinant = first[0] * second[1] - first[1] * second[0]

    if abs(determinant) < 1e-14:
        raise ValueError("Ошибка newtonStep: матрица Якоби близка к вырожденной.")

    delta1 = (-values[0] * second[1] + first[1] * values[1]) / determinant
    delta2 = (-first[0] * values[1] + values[0] * second[0]) / determinant

    return [point[0] + delta1, point[1] + delta2]

def solve(start, a, epsilon, maxIterations, method, q=None):
    point = start.copy()
    history = [(0, point, None, vectorInfNorm(functions(point, a)))]

    for iteration in range(1, maxIterations + 1):
        newPoint = simpleStep(point, a) if method == "simple" else newtonStep(point, a)
        difference = vectorInfNorm([newPoint[i] - point[i] for i in range(2)])
        nevyazka = vectorInfNorm(functions(newPoint, a))
        history.append((iteration, newPoint, difference, nevyazka))

        error = q / (1 - q) * difference if method == "simple" else difference

        if error <= epsilon and nevyazka <= epsilon:
            return newPoint, history
        
        point = newPoint

    raise ValueError(f"Метод {method} не достиг заданной точности.")


def printHistory(title, history, reference):
    print("\n" + title)
    print("  k           x1           x2       изменение       невязка     ошибка к эталону")
    for iteration, point, difference, nevyazka in history:
        if difference is None:
            differenceStr = "—"
        else:
            differenceStr = f"{difference:.3e}"

        error = vectorInfNorm([point[i] - reference[i] for i in range(2)])

        iterationStr = f"{iteration:3d}"
        xStr = f"{point[0]:11.8f}"
        yStr = f"{point[1]:11.8f}"
        nevyazkaStr = f"{nevyazka:12.3e}"
        errorStr = f"{error:16.3e}"

        line = "  ".join([
            iterationStr,
            xStr,
            yStr,
            differenceStr.rjust(14),
            nevyazkaStr,
            errorStr,
        ])
        print(line)

def main():
    if len(sys.argv) != 2:
        print("Использование: python system.py 2.json")
        return

    data = loadInputData(sys.argv[1])
    a = data["a"]
    start = data["start"]
    epsilon = data["epsilon"]
    q = data["q"]
    maxIterations = data["maxIterations"]

    if not all(math.isfinite(value) for value in (a, *start, epsilon, q)):
        raise ValueError("Ошибка: все числовые параметры должны быть конечными.")
    
    if type(maxIterations) is not int or maxIterations <= 0:
        raise ValueError("Ошибка: maxIterations должно быть положительным целым числом.")
    
    if len(start) != 2 or a <= 0 or not 0 < q < 1 or not 0 < epsilon < 1:
        raise ValueError("Ошибка: проверьте a, start, epsilon и q в JSON.")
    
    functions(start, a)

    reference, _ = solve(start, a, 1e-13, maxIterations, "newton")
    
    simpleRoot, simpleHistory = solve(start, a, epsilon, maxIterations, "simple", q)
    
    newtonRoot, newtonHistory = solve(start, a, epsilon, maxIterations, "newton")

    print(f"Эталонная пара: ({reference[0]:.13f}, {reference[1]:.13f})")
    printHistory("Метод простых итераций", simpleHistory, reference)
    printHistory("Метод Ньютона", newtonHistory, reference)
    print(f"\nЧисло итераций: простые — {len(simpleHistory)-1}, Ньютон — {len(newtonHistory)-1}.")
    print(f"Полученные решения: {simpleRoot} и {newtonRoot}")

if __name__ == "__main__":
    main()
