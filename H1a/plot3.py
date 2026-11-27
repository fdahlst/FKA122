import numpy as np
import matplotlib.pyplot as plt

temperatures = np.genfromtxt("data/task3/ex3_temperatures.csv",delimiter=",")
pressures = np.genfromtxt("data/task3/ex3_pressures.csv",delimiter=",")
end_positions = np.genfromtxt("data/task3/end_positions.csv", delimiter = ",")
positions = np.genfromtxt("data/task3/positions.csv", delimiter = ",")


t = np.linspace(0,5,len(temperatures))

plt.plot(t,temperatures)
plt.xlabel("time (ps)")
plt.ylabel("temperature (K)")
plt.savefig("plots/task3/temp.png")
plt.cla()

plt.plot(t,pressures)
plt.xlabel("time (ps)")
plt.ylabel("pressures (eV/Å^3)")
plt.savefig("plots/task3/pressures.png")
plt.cla()


plt.figure(figsize=(12, 4))
plt.subplot(1,3,1)
plt.plot(t, positions[:,0]-positions[0,0])
plt.plot(t, positions[:,3]-positions[0,3])
plt.plot(t, positions[:,6]-positions[0,6])
plt.ylabel("displacement (Å)")
plt.xlabel("time (ps)")
plt.title("x")

plt.subplot(1,3,2)
plt.plot(t, positions[:,1]-positions[0,1])
plt.plot(t, positions[:,4]-positions[0,4])
plt.plot(t, positions[:,7]-positions[0,7])
plt.ylabel("displacement (Å)")
plt.xlabel("time (ps)")
plt.title("y")

plt.subplot(1,3,3)
plt.plot(t, positions[:,2]-positions[0,2])
plt.plot(t, positions[:,5]-positions[0,5])
plt.plot(t, positions[:,8]-positions[0,8])
plt.ylabel("displacement (Å)")
plt.xlabel("time (ps)")
plt.title("z")

plt.tight_layout()
plt.savefig("plots/task3/displacement.png")

