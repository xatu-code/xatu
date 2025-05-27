import matplotlib.pyplot as plt
import numpy as np
import sys

# Egap = sys.argv[1]

# -----------------------------------------------------------
# Configure plot
# -----------------------------------------------------------

exciton_file = sys.argv[1]
sp_file      = sys.argv[2]
plotType     = int(sys.argv[3])
# states_file  = "in2se3.eigval"
fontsize     = 15
markersize   = 110  # Adjust size of points in plot
latex        = True # Set to True if latex is installed
if latex:
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif"
    })
plt.rcParams["axes.labelsize"] = fontsize
plt.rcParams["axes.titlesize"] = fontsize
plt.rcParams["xtick.labelsize"] = fontsize
plt.rcParams["ytick.labelsize"] = fontsize
plt.rcParams["legend.fontsize"] = fontsize
plt.rcParams["font.size"] = fontsize

# -----------------------------------------------------------
# Extract conductivies from files
# -----------------------------------------------------------

file = open(exciton_file, "r")
energy, sigma_sp_xx, sigma_sp_yy, sigma_sp_zz = [], [], [], []
sigma_xx, sigma_yy, sigma_zz = [], [], []
E, Vx, Vy, Vz = [], [], [], []

flag = 1
for line in file.readlines():
    line = [float(value) for value in line.split()]
    if not line: # Skip reading exciton velocity matrix elements
        flag = 2
        continue
    if flag == 1:
        energy.append(line[0])
        sigma_xx.append(line[1])
        sigma_yy.append(line[5])
        sigma_zz.append(line[9])
    else:
        E.append(line[0])
        Vx.append(np.abs(line[1]+line[2]*1j)**2)
        Vy.append(np.abs(line[3]+line[4]*1j)**2)
        Vz.append(np.abs(line[5]+line[6]*1j)**2)

E  = np.array(E)
Vx = np.array(Vx)
Vy = np.array(Vy)
Vz = np.array(Vz)

file = open(sp_file, "r")
for line in file.readlines():
    line = [float(value) for value in line.split()]
    sigma_sp_xx.append(line[1])
    sigma_sp_yy.append(line[5])
    sigma_sp_zz.append(line[9])

# -----------------------------------------------------------
# Plot states
# -----------------------------------------------------------

# plot conductivity
if plotType == 1:
    fig, ax = plt.subplots(figsize=(6.4,4.8))
    ax.plot(energy, sigma_xx, "r-", label="BSE_xx",alpha=0.6)
    ax.plot(energy, sigma_yy, "orange", label="BSE_yy",alpha=0.6)
    ax.plot(energy, sigma_zz, "g-", label="BSE_zz",alpha=0.6)
    ax.plot(energy, sigma_sp_xx, "b--", label="IPA_xx")
    ax.plot(energy, sigma_sp_yy, "c--", label="IPA_yy")
    ax.plot(energy, sigma_sp_zz, c="navy", ls='--', label="IPA_zz")

    ax.set_xlabel("E (eV)", fontsize=18)
    ax.set_ylabel(r"$\sigma$ ($e^2/\hbar$)", fontsize=18)
    ax.legend(fontsize=fontsize/1.5)
    ax.set_xlim(energy[0],energy[-1])
    plt.tight_layout()
    plt.savefig("conductivity.png", dpi = 400)

# plot conductivity + oscillator strengths
elif plotType == 2:
    fig, ax = plt.subplots(2,1, figsize=(6.4,4.8))
    ax[0].plot(energy, sigma_xx, "r-", label="BSE_xx")
    ax[0].plot(energy, sigma_yy, "orange", label="BSE_yy")
    ax[0].plot(energy, sigma_zz, "g-", label="BSE_zz")
    ax[0].plot(energy, sigma_sp_xx, "b--", label="IPA_xx")
    ax[0].plot(energy, sigma_sp_yy, "c--", label="IPA_yy")
    ax[0].plot(energy, sigma_sp_zz, "navy", ls='--', label="IPA_zz")

    # ax[0].set_xlabel("E (eV)", fontsize=18)
    ax[0].set_ylabel(r"$\sigma$ ($e^2/\hbar$)", fontsize=18)

    # Create masks for your conditions
    mask_vx = Vx > 0.1
    mask_vy = Vy > 0.1
    mask_vz = Vz > 0.5

    # Plot all Vx bars with a single label
    if np.any(mask_vx):
        ax[1].bar(E[mask_vx], E[mask_vx]*Vx[mask_vx], 
                 width=0.03, color="r", alpha=0.3, label='xx')

    # Plot all Vy bars with a single label
    if np.any(mask_vy):
        ax[1].bar(E[mask_vy], E[mask_vy]*Vy[mask_vy], 
                 width=0.03, color="orange", alpha=0.3, label='yy')

    # Plot all Vy bars with a single label
    if np.any(mask_vz):
        ax[1].bar(E[mask_vz], E[mask_vz]*Vz[mask_vz], 
                 width=0.03, color="green", alpha=0.3, label='zz')

    ax[1].set_ylabel(f'$|V^\\alpha|^2$', fontsize=18) # $A_{vc}(\mathbf{k})v_{vc}^{\alpha}(\mathbf{k})$
    ax[1].set_xlabel("E (eV)", fontsize=18)

    for axs in ax:
        axs.set_xlim(energy[0],energy[-1])
        axs.tick_params(labelsize=12)

    plt.subplots_adjust(wspace=0, hspace=0.1)

    ax[0].legend(fontsize=fontsize/1.5)
    ax[1].legend(fontsize=fontsize/1.5)
    plt.tight_layout()
    plt.savefig("conductivity_oscillators.png", dpi = 400)
plt.show()