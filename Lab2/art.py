import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt

from toolbox import loadInputData


def makePoints(left, right, count=400):
    return [left + (right - left) * i / count for i in range(count + 1)]


def drawSolo(data, outputPath):
    coefficient = data["coefficient"]
    constant = data["constant"]
    left, right = data["interval"]
    start = data["start"]

    width = right - left
    xValues = makePoints(left - 3 * width, right + 3 * width)
    tangentValues = [math.tan(x) for x in xValues]
    parabolaValues = [coefficient * x * x - constant for x in xValues]

    plt.figure(figsize=(8, 5))
    plt.plot(xValues, tangentValues, label="tan(x)")
    plt.plot(xValues, parabolaValues, label=f"{coefficient}x² − {constant}")
    plt.axvspan(left, right, color="green", alpha=0.15, label="Интервал с корнем")
    plt.axvline(start, color="black", linestyle="--", label=f"Начало x₀ = {start}")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outputPath)
    plt.close()


def drawSystem(data, outputPath):
    a = data["a"]
    startX1, startX2 = data["start"]

    xValues = makePoints(startX1 * 0.7, startX1 * 1.3)
    firstValues = [10 ** ((x * x - 1) / a) for x in xValues]
    secondValues = [(x * x + a) / (a * x) for x in xValues]

    plt.figure(figsize=(8, 5))
    plt.plot(xValues, firstValues, label="x₂ = 10^((x₁² − 1)/a)")
    plt.plot(xValues, secondValues, label="x₂ = (x₁² + a)/(a·x₁)")
    plt.scatter(startX1, startX2, color="black", label="Начальная точка")
    plt.xlabel("x₁")
    plt.ylabel("x₂")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outputPath)
    plt.close()


def main():
    if len(sys.argv) != 3:
        print("Использование: python art.py 1.json 2.json")
        return

    soloPath = Path(sys.argv[1])
    systemPath = Path(sys.argv[2])
    drawSolo(loadInputData(soloPath), soloPath.with_name("solo.png"))
    drawSystem(loadInputData(systemPath), systemPath.with_name("system.png"))
    print("Графики сохранены как .png.")


if __name__ == "__main__":
    main()
