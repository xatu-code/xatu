import numpy as np
import matplotlib.pyplot as plt
import sys
# -----------------------------------------------------------
# Configure plot
# -----------------------------------------------------------

fontsize   = 15
latex      = True # Set to True if latex is installed
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
# First extract states from file
# -----------------------------------------------------------
filename = sys.argv[1]
file = open(filename, "r")
lines = file.readlines()
kpoints = []
ks = []
coefs = []
states = []
for line in lines:
    if line[0] == "k":
        continue
    if line[0] == "#":
        states.append(np.array(coefs))
        ks.append(kpoints)
        kpoints = []
        coefs = []
        continue
    line = line.split()
    kpoint = [float(num) for num in line[0:3]]
    kpoints.append(kpoint)
    coefs.append(float(line[-1]))

kpoints = np.array(ks[0])
file.close()

# -----------------------------------------------------------
# Sum densities over degeneracies
# -----------------------------------------------------------

degeneracies = [1, 1, 1, 1, 1, 1, 1, 1]
wfs = []
counter = 0
for deg in degeneracies:
    state = np.zeros(states[0].shape)
    for i in range(deg):
        state += states[counter + i]
    counter += deg
    wfs.append(state)

# -----------------------------------------------------------
# Plot states
# -----------------------------------------------------------
fig, ax = plt.subplots(4, 2, figsize=(6, 8), sharex=True, sharey=True)

ax[0, 0].tripcolor(kpoints[:, 0], kpoints[:, 1], wfs[0], cmap="viridis")
ax[0, 1].tripcolor(kpoints[:, 0], kpoints[:, 1], wfs[1], cmap="viridis")
ax[1, 0].tripcolor(kpoints[:, 0], kpoints[:, 1], wfs[2], cmap="viridis")
ax[1, 1].tripcolor(kpoints[:, 0], kpoints[:, 1], wfs[3], cmap="viridis")
ax[2, 0].tripcolor(kpoints[:, 0], kpoints[:, 1], wfs[4], cmap="viridis")
ax[2, 1].tripcolor(kpoints[:, 0], kpoints[:, 1], wfs[5], cmap="viridis")
ax[3, 0].tripcolor(kpoints[:, 0], kpoints[:, 1], wfs[6], cmap="viridis")
ax[3, 1].tripcolor(kpoints[:, 0], kpoints[:, 1], wfs[7], cmap="viridis")

for axis_row in ax:
    for axis in axis_row:
        axis.axis("square")
# ax[2, 1].axis("off")

ylim = max(kpoints[:, 1])    
ax[0, 0].set_ylim([-ylim, ylim])
xlim = max(kpoints[:, 0])
ax[0, 0].set_xlim([-xlim, xlim])
# ax[0].set_xlabel(r"$k_x$ (\AA$^{-1}$)")

ax[0, 0].set_ylabel(r"$k_y$ (\AA$^{-1}$)")
ax[1, 0].set_ylabel(r"$k_y$ (\AA$^{-1}$)")
ax[2, 0].set_ylabel(r"$k_y$ (\AA$^{-1}$)")
ax[3, 0].set_ylabel(r"$k_y$ (\AA$^{-1}$)")
ax[3, 0].set_xlabel(r"$k_x$ (\AA$^{-1}$)")
ax[3, 1].set_xlabel(r"$k_x$ (\AA$^{-1}$)")
ax[1, 1].set_xticks([-xlim, 0, xlim])

ax[1,1].tick_params(labelbottom=True)

plt.tight_layout()
plt.savefig("kwf.png", dpi=600)
plt.savefig("kwf.pdf", dpi=600)
# plt.show()