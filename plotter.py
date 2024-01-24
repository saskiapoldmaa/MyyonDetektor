import matplotlib.pyplot as plt
import numpy as np

# Read data from the file
foon_path = '/Users/saskiapoldmaa/Documents/UT/foon200'
with open(foon_path, 'r') as foon_file:
    for _ in range(6):
        next(foon_file)
    foon = [float(line.split()[5]) for line in foon_file]

gamma_path = '/Users/saskiapoldmaa/Documents/UT/gammafoon200'
with open(gamma_path, 'r') as gamma_file:
    for _ in range(6):
        next(gamma_file)
    gamma = [float(line.split()[5]) for line in gamma_file]

# Plot histogram
plt.hist(gamma, weights=np.ones(len(gamma))/len(gamma), bins=np.linspace(20, 70, 100), color='red', edgecolor='black', alpha=0.7)
plt.hist(foon, weights=np.ones(len(foon))/len(foon), bins=np.linspace(20, 70, 100), color='blue', edgecolor='black', alpha=0.7)
plt.title('Spekter piirväärtusel 200ADC')
plt.xlabel('Pingeamplituud [mV]')
plt.ylabel('Osakaal spektrist')
plt.legend(['Th-232 gamma-allikaga','Ilma gamma-allikata'])
plt.grid(True)
plt.show()
