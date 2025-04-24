import numpy as np
import matplotlib.pyplot as plt

a = np.loadtxt("r_s_out.txt");

plt.plot(a[:, 0], a[:, 1])
plt.grid()
plt.show()