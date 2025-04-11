import numpy as np
import matplotlib.pyplot as plt


def plot_from_file(file_name):
    data = []

    try:
        with open(file_name, "r") as f:
            file_data = np.loadtxt(f).transpose()
            t = file_data[0]
            data = file_data[1:,:]
    except FileNotFoundError:
        print(f"Файл {file_name} не найден. Пропускаем.")
    except Exception as e:
        print(f"Ошибка при чтении файла {file_name}: {e}")

    # Отрисовка графиков
    plt.figure(figsize=(9, 5))
    for i, var in enumerate(data):
    	plt.plot(t, var, linewidth=2, label=f"x[{i+1}]")

    plt.xlabel(r"Ось $r$", labelpad=-4, fontsize=12)
    plt.ylabel(r"$u(r)$", fontsize=12)
    plt.title(r"Численное решение задачи Коши")
    plt.grid(True)
    plt.legend(loc='best')
    plt.grid(which="major", linewidth=1)
    plt.grid(which="minor", linewidth=0.4)
    plt.minorticks_on()
    plt.tight_layout()
    plt.show()

plot_from_file("output.txt")