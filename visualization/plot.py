import matplotlib.pyplot as plt
import numpy as np

def plot_2d():
    file_name = "output.txt"

    with open(file_name, "r") as f:
        file_data = np.loadtxt(f).transpose()
        t = file_data[0]
        u = file_data[1]
        v = file_data[2]

    plt.plot(u, v, linewidth=2, label=f"c")

    # Настраиваем первый график
    plt.xlabel(r"Ось $r$", labelpad=-4, fontsize=12)
    plt.ylabel(r"$u(r)$", fontsize=12)
    plt.title("Сравнение аналитического и численного решений")
    plt.grid(True)
    plt.legend(loc='best')

    plt.tight_layout()

    plt.show()

def plot_1d():
    print("enter theta")
    theta = input()
    file_name = "beta_" + theta + ".txt"

    with open(file_name, "r") as f:
        file_data = np.loadtxt(f).transpose()
        M = file_data[0]
        beta = file_data[1]

    plt.plot(M, beta / np.pi * 180.0, linewidth=2, label=f"$\\theta_0 = {theta}\degree$")

    # Настраиваем первый график
    plt.xlabel(r"$M_\infty$", fontsize=12)
    plt.ylabel(r"$\beta$", fontsize=12)
    # plt.title("Сравнение аналитического и численного решений")
    plt.grid(True)
    plt.xlim((1, 5))
    plt.ylim((0, 80))
    plt.legend(loc='best')

    plt.tight_layout()

    plt.show()

plot_1d()