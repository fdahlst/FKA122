import matplotlib.pyplot as plt
import numpy as np

pot_energy = np.genfromtxt("ex1_pot_energy.csv", delimiter = ",")
lattice_const = np.array([4.0,4.01,4.02,4.03,4.04,4.05,4.06,4.07,4.08])

#Fit polynomial
fit = np.polyfit(lattice_const**3,pot_energy, 2)

x = np.linspace(lattice_const[0]**3,lattice_const[-1]**3)
y = fit[0]*x**2+fit[1]*x+fit[2]

plt.plot(x,y,label="quadratic fit")
plt.scatter(lattice_const**3,pot_energy, label = "data")
plt.xlabel("Unit cell volume(Å^3)")
plt.ylabel("Energy(Ev/unit cell)")
plt.legend()
plt.tight_layout()
plt.savefig("plots/1.png")

min_energy = min(pot_energy);
print(f"Mininum energy{min_energy}eV/Unit Cell with lattice constant {lattice_const[np.where(pot_energy == min_energy)][0]}Å")




