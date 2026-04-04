import numpy as np
import matplotlib.pyplot as plt

data_root_path = '../output_results'
a = np.loadtxt(data_root_path + "r_s_out.txt")
x = np.linspace(np.min(a[:, 0]), np.max(a[:, 0]))

plt.plot(a[:, 0], a[:, 1])
plt.plot(x, x * np.tan(0.440851))
plt.grid()
plt.show()