import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
import numpy as np

for nvac in range(2,6):

    # Import file
    filename = f"/home/zichuan/openFLY/data/formation_energy_{nvac}V.txt"

    cellsize = np.loadtxt(filename, delimiter = " ", usecols = 0)
    Ef_cluster = np.loadtxt(filename, delimiter = " ", usecols = 1)
    Ef_clusterperfect = np.loadtxt(filename, delimiter = " ", usecols = 2)

    print(cellsize)
    print(Ef_cluster)
    print(Ef_clusterperfect)


    # Create plot
    fig, axs = plt.subplots(2, sharex=True)
    plt.suptitle(f"Formation energy of {nvac}-vacancy cluster", fontsize = 20)
    plt.xlabel("number of atoms", fontsize = 15)

    # Plot no.1 for Ef_cluster
    labelpoint = 10
    axs[0].plot(cellsize, Ef_cluster, marker='.')
    axs[0].plot(cellsize[labelpoint], Ef_cluster[labelpoint], marker = "o", color = "r")
    axs[0].set_title(f"from {nvac} vacancies", fontsize = 17)
    axs[0].set_ylabel("Ef / eV", fontsize = 15)

    # Plot no.2 for Ef_clusterperfect
    axs[1].plot(cellsize, Ef_clusterperfect, marker='.')
    axs[1].plot(cellsize[labelpoint], Ef_clusterperfect[labelpoint], marker = "o", color = "r")
    axs[1].set_title("from perfect lattice", fontsize = 17)
    axs[1].set_ylabel("Ef / eV", fontsize = 15)

    # Label a data point
    axs[0].annotate("   future simulation (8x8x8 cell)", (cellsize[labelpoint], Ef_cluster[labelpoint]), xycoords = 'data', fontsize = 13)
    axs[1].annotate("   future simulation (8x8x8 cell)", (cellsize[labelpoint], Ef_clusterperfect[labelpoint]), xycoords = 'data', fontsize = 13)

    # Save images
    fig.set_size_inches(9,9)
    plt.savefig(f"/home/zichuan/openFLY/data/figures/{nvac}V_formation_energy.png", dpi = 80)