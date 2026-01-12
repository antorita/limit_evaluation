import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("lim_output.txt")
sigma = data[:,1]
masses = data[:,0]

plt.xlabel(r"$m_{\chi}\ \mathrm{[GeV/c^2]}$")
plt.ylabel(r"$\sigma\ \mathrm{[cm^2]}$")

plt.loglog(masses, sigma, marker='o', linestyle='-')
plt.grid(True, which="both", ls="--")
plt.show()
