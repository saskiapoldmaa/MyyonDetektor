import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import random
from scipy.integrate import odeint

g=9.80665
R=8.31447
T0=288.15
P0=101325
M=0.0289644
L = 0.0065
E0=105.6584
h = np.linspace(15, 0, 200)


def dEdh(E, h):
    T=T0-1000*L*h
    P=P0*(T/T0)**(g*M/R/L)
    dens=M*P/T/R
    diff= dens*15.3289*(8.386+math.log(1-E0**2/E**2)-2*math.log(E0/E)+E0**2/E**2)*(E**2)/(E**2-E0**2)
    return diff

plt.figure(figsize=(8, 6))

ylist=np.random.normal(6000, 1000, 50)
for i in range(1,50):
    y0=ylist[i]
    E_solution = odeint(dEdh, y0, h)
    random_color = (random.random(),random.random(),random.random())
    plt.plot(h, E_solution, color=random_color)

plt.gca().invert_xaxis()
plt.xlabel('Kõrgus [km]')
plt.gca().set_ylim(0,1e4)
plt.ylabel('Energia [MeV]')
plt.title('Müüonite energia ning kõrgus')
plt.legend()
plt.grid(True)
plt.show()


