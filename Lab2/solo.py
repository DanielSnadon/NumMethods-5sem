import math
import sys
from toolbox import *

def function(x, coefficient, constant):
    return math.tan(x) - coefficient * x * x + constant

def derivative(x, coefficient):
    return 1.0 / math.cos(x) ** 2 - 2.0 * coefficient * x

def fi(x, coefficient, constant):
    result = (math.tan(x) + constant) / coefficient

    if result < 0:
        raise ValueError("Ошибка Fi: подкоренное выражение стало отрицательным.")
    
    return math.sqrt(result)

def simpleIterations(start, coefficient, constant, epsilon, q, maxIterations):
    x = start
    history = [(0, x, None, abs(function(x, coefficient, constant)))]

    for iteration in range(1, maxIterations + 1):
        newX = fi(x, coefficient, constant)
        difference = abs(newX - x)
        nevyazka = abs(function(newX, coefficient, constant))
        history.append((iteration, newX, difference, nevyazka))

        if q / (1.0 - q) * difference <= epsilon and nevyazka <= epsilon:
            return newX, history
        
        x = newX

    raise ValueError("Ошибка SimpleIterations: простые итерации не достигли заданной точности.")

def newton(start, coefficient, constant, epsilon, maxIterations):
    x = start
    history = [(0, x, None, abs(function(x, coefficient, constant)))]

    for iteration in range(1, maxIterations + 1):
        angle = derivative(x, coefficient)
        if abs(angle) < 1e-14:
            raise ValueError("Ошибка Newton: производная близка к нулю.")

        newX = x - function(x, coefficient, constant) / angle
        difference = abs(newX - x)
        nevyazka = abs(function(newX, coefficient, constant))
        history.append((iteration, newX, difference, nevyazka))

        if difference <= epsilon and nevyazka <= epsilon:
            return newX, history
        
        x = newX

    raise ValueError("Ошибка Newton: метод не достиг заданной точности.")

def printHistory(title, history, reference):
    print("\n" + title)
    print("  k              x           изменение         |f(x)|       ошибка к эталону")
    
    for iteration, x, difference, nevyazka in history:
        if difference is None:
            differenceStr = "—"
        else:
            differenceStr = f"{difference:.3e}"

        iterationStr = f"{iteration:3d}"
        xStr = f"{x:14.10f}"
        nevyazkaStr = f"{nevyazka:12.3e}"
        errorStr = f"{abs(x - reference):16.3e}"

        line = "  ".join([iterationStr, xStr, differenceStr.rjust(16), nevyazkaStr, errorStr])
        print(line)

def main():
    if len(sys.argv) != 2:
        print("Использование: python solo.py 1.json")
        return

    data = loadInputData(sys.argv[1])
    coefficient = data["coefficient"]
    constant = data["constant"]
    left, right = data["interval"]
    start = data["start"]
    epsilon = data["epsilon"]
    q = data["q"]
    maxIterations = data["maxIterations"]

    if not all(math.isfinite(value) for value in (coefficient, constant, left, right, start, epsilon, q)):
        raise ValueError("Ошибка: все числовые параметры должны быть конечными.")
    
    if type(maxIterations) is not int or maxIterations <= 0:
        raise ValueError("Ошибка: maxIterations должно быть положительным целым числом.")
    
    if not (coefficient > 0 and 0 < left < right < math.pi / 2):
        raise ValueError("Ошибка: нужны c > 0 и интервал внутри (0, π/2).")
    
    if not (left <= start <= right and 0 < epsilon < 1 and 0 < q < 1):
        raise ValueError("Ошибка: проверьте start, epsilon и q в JSON.")
    
    if function(left, coefficient, constant) * function(right, coefficient, constant) > 0:
        raise ValueError("Ошибка: на концах интервала функция имеет одинаковые знаки.")

    reference, _ = newton(start, coefficient, constant, 1e-13, maxIterations)

    simpleRoot, simpleHistory = simpleIterations(start, coefficient, constant, epsilon, q, maxIterations)

    newtonRoot, newtonHistory = newton(start, coefficient, constant, epsilon, maxIterations)

    print(f"Эталонный корень = {reference:.13f}")
    printHistory("Метод простых итераций", simpleHistory, reference)
    printHistory("Метод Ньютона", newtonHistory, reference)
    print(f"\nЧисло итераций: простые — {len(simpleHistory)-1}, Ньютон — {len(newtonHistory)-1}.")
    print(f"Полученные корни: {simpleRoot:.12f} и {newtonRoot:.12f}")

if __name__ == "__main__":
    main()
