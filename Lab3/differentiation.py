import json
from pathlib import Path


def chooseNodes(nodes, xStar):
    if len(nodes) < 3:
        raise ValueError("Ошибка chooseNodes: нужно не менее трёх узлов.")
    if any(nodes[i] >= nodes[i + 1] for i in range(len(nodes) - 1)):
        raise ValueError("Ошибка chooseNodes: узлы должны идти в строго возрастающем порядке.")
    if not nodes[0] <= xStar <= nodes[-1]:
        raise ValueError("Ошибка chooseNodes: точка вычисления должна лежать внутри таблицы.")

    middle = 1
    for i in range(2, len(nodes) - 1):
        if abs(nodes[i] - xStar) < abs(nodes[middle] - xStar):
            middle = i 

    return middle - 1, middle, middle + 1


def derivatives(x0, x1, x2, y0, y1, y2, xStar):
    f1 = (y1 - y0) / (x1 - x0)
    f2 = (y2 - y1) / (x2 - x1)
    f3 = (f2 - f1) / (x2 - x0)

    # zeroDerivative = y0 + f1*(x-x0) + f3*(x-x0)*(x-x1).
    firstDerivative = f1 + f3 * (2 * xStar - x0 - x1)
    secondDerivative = 2 * f3

    return firstDerivative, secondDerivative


def main():
    inputPath = Path(__file__).with_name("4.json")
    with inputPath.open(encoding="utf-8") as file:
        data = json.load(file)

    nodes = data["nodes"]
    values = data["values"]
    xStar = data["xStar"]
    if len(nodes) != len(values):
        raise ValueError("Количество x и y должно совпадать.")

    print("Исходная таблица:")
    print(f"{'i':>3} {'x_i':>10} {'y_i':>10}")
    for i in range(len(nodes)):
        print(f"{i:3d} {nodes[i]:10.6g} {values[i]:10.6g}")

    print(f"\nТочка вычисления: x* = {xStar:g}")
    i, j, k = chooseNodes(nodes, xStar)
    print("Использованные узлы:")
    for index in (i, j, k):
        print(f"  i={index}: x={nodes[index]:g}, y={values[index]:g}")

    first, second = derivatives(
        nodes[i], nodes[j], nodes[k],
        values[i], values[j], values[k], xStar,
    )
    print(f"\nf'(x*)  ≈ {first:.10g}")
    print(f"f''(x*) ≈ {second:.10g}")


if __name__ == "__main__":
    main()
