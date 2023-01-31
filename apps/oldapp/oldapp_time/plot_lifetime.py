import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import statsmodels.api as sm

filename = "/home/zichuan/openFLY/data/V1H1_lifetime1.txt"

temperature = np.loadtxt(filename, delimiter = " ", usecols = 0)
meanlifetime = np.loadtxt(filename, delimiter = " ", usecols = 1)
stdev = np.loadtxt(filename, delimiter = " ", usecols = 2)

length = len(temperature)
print(length)   

inversetemp = np.zeros(length)
lnmean = np.zeros(length)
errlifetime = np.zeros(length)

for i in range(length):
    inversetemp[i] = 1 / temperature[i]
    lnmean[i] = np.log(meanlifetime[i])
    errlifetime[i] = stdev[i] / (np.sqrt(length) * meanlifetime[i])
    
print("temperature: ", temperature)
print("meanlifetime: ", meanlifetime)
print("stdev: ", stdev)
print("inversetemp: ", inversetemp)
print("lnmean", lnmean)
print("errlifetime", errlifetime)

# Create a plot
fig, axs =  plt.subplots()
axs.scatter(inversetemp, lnmean, color = "blue", marker = ".")

# Linear regression calculation and plotting regression line
m, b = np.polyfit(inversetemp, lnmean, 1)
m1 = int(m)
b1 = round(b, 2)
plt.plot(inversetemp, m*inversetemp + b, ":", label = rf"ln $\tau$ = {m1}/T - {-b1}")
plt.legend(loc="lower right", fontsize = 13)

# Adding error bars
plt.errorbar(inversetemp, lnmean, yerr=errlifetime, linestyle = "None", color = "blue")

# Set x and y ranges
plt.xlim([0, 0.006])

axs.set_ylabel(r"ln( $\tau$ / s )", fontsize = 15)
axs.set_xlabel("(1/T) / (1/K)", fontsize = 15)

# Saving the image
fig.set_size_inches(10,7)
plt.savefig("/home/zichuan/openFLY/data/V1H1lifetime.eps", format = "eps")

model = sm.OLS(lnmean, sm.add_constant(inversetemp))
results = model.fit()

print(results.params)
print(results.summary())