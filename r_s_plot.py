import numpy as np
import matplotlib.pyplot as plt

a = np.loadtxt("r_s_out.txt")
x = np.linspace(np.min(a[:, 0]), np.max(a[:, 0]))

plt.plot(a[:, 0], a[:, 1])
plt.plot(x, x * np.tan(0.440851))
plt.grid()
plt.show()