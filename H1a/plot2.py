import numpy as np
import matplotlib.pyplot as plt

kinetic_energies = np.genfromtxt("ex2_kin_energy.csv",delimiter=",")
potential_energies = np.genfromtxt("ex2_pot_energy.csv",delimiter=",")
total_energies = np.genfromtxt("ex2_tot_energy.csv",delimiter=",")

kb = 8.6171e-5

t = np.linspace(0,1,len(kinetic_energies))

plt.plot(t,kinetic_energies, label = "kinetic")
plt.plot(t,potential_energies, label = "potential")
plt.plot(t,total_energies, label = "total")
plt.legend()
#plt.savefig("plots/plot2.png")

print(f"average temperature: {2*kinetic_energies.mean()/(3*kb*256)} K")
