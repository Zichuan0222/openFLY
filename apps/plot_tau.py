import numpy as np
import math
from time import gmtime, strftime
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import statsmodels.api as sm


# ----------------------- Min and Max temperatures ---------------------- #
minT = 3
maxT = 8                        # maxT = maximum temp / 100 + 1
Tsize = maxT - minT

temperature = np.zeros(Tsize)
inverse_temp = np.zeros(Tsize)
taumean = np.zeros(Tsize)
ln_taumean = np.zeros(Tsize)
stdev_tau = np.zeros(Tsize)
stdev_lntau = np.zeros(Tsize)

for i in range (minT, maxT):
    temperature[i - minT] = (i * 100)


for i in range(Tsize):

    # ----------------------------------- Cluster type --------------------------------- #
    clustertype = "V3H_Vdis"


    datadirectory = "/home/zichuan/openFLY/build/data/"       # data txt file directory
    txtname = clustertype + "/" + str(int(temperature[i]))         # txt file name
    filename = datadirectory + txtname + "K.txt"

    tau = np.loadtxt(filename, usecols = 0)
    tausize = len(tau)                                        # no. of data point
    
    inverse_temp[i] = 1 / temperature[i]
    taumean[i] = sum(tau) / tausize
    ln_taumean[i] = np.log(taumean[i])
    squareterm = np.zeros(tausize)

    # Calculate std dev
    for j in range (tausize):
        squareterm[j] = (tau[j] - taumean[i]) * (tau[j] - taumean[i])

    stdev_tau[i] = np.sqrt(sum(squareterm) / tausize) 
    stdev_lntau[i] = stdev_tau[i] / (np.sqrt(tausize) * taumean[i])

    print(temperature[i], taumean[i], stdev_tau[i], ln_taumean[i], stdev_lntau[i])
    
countfilename = datadirectory + clustertype + "/dissociation_count.txt"
iteration = np.loadtxt(countfilename, usecols = 1)
no_iteration = int(np.max(iteration))


# # temperature = np.loadtxt(filename, delimiter = " ", usecols = 0)
# # meanlifetime = np.loadtxt(filename, delimiter = " ", usecols = 1)
# # stdev = np.loadtxt(filename, delimiter = " ", usecols = 2)



######## Create the plot #########
fig, axs =  plt.subplots()
axs.scatter(inverse_temp, ln_taumean, color = "blue", marker = ".")

# Linear regression calculation and plotting regression line
m, b = np.polyfit(inverse_temp, ln_taumean, 1)
m1 = int(m)
b1 = round(b, 2)
plt.plot(inverse_temp, m*inverse_temp + b, ":", label = rf"ln $\tau$ = {m1}/T - {-b1}")
plt.legend(loc="lower right", fontsize = 13)

# Adding error bars 
plt.errorbar(inverse_temp, ln_taumean, yerr=stdev_lntau, linestyle = "None", color = "blue")

# Set x and y ranges
plt.xlim([0.0005, 0.0045])

# Set axis labels and title
axs.set_ylabel(r"ln($\tau$ / s)", fontsize = 15)
axs.set_xlabel("(1/T) / (1/K)", fontsize = 15)
    
plt.title(clustertype + f"\n{no_iteration} iterations at each temperature", fontsize = 16   )

# Saving the image and name after current time
figdirectory = datadirectory +  clustertype + "/"
figname = clustertype + "_" + str(no_iteration) + ".png"                               # Change extension here, e.g. .png or .eps

fig.set_size_inches(10,7)
plt.savefig(figdirectory + figname)
#plt.savefig("/home/zichuan/openFLY/data/V2H.eps", format = "eps")

model = sm.OLS(ln_taumean, sm.add_constant(inverse_temp))
results = model.fit()
 
print(results.params)
print(results.summary())