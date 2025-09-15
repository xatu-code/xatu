#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt

def configure_matplotlib(fontsize=15, markersize=110, use_latex=True):
    """Set global matplotlib rcParams."""
    if use_latex:
        plt.rcParams.update({
            "text.usetex": True,
            "font.family": "serif"
        })
    plt.rcParams["axes.labelsize"]   = fontsize
    plt.rcParams["axes.titlesize"]   = fontsize
    plt.rcParams["xtick.labelsize"]  = fontsize
    plt.rcParams["ytick.labelsize"]  = fontsize
    plt.rcParams["legend.fontsize"]  = fontsize
    plt.rcParams["font.size"]        = fontsize


def read_exciton_file(path):
    """
    Read exciton file: energy and BSE conductivities only.
    Returns:
      energy, sigma_xx, sigma_yy, sigma_zz
    """
    energy = []
    sigma_xx = []
    sigma_yy = []
    sigma_zz = []

    with open(path, 'r') as f:
        for raw in f:
            parts = raw.split()
            if not parts:
                break
            vals = [float(v) for v in parts]
            energy.append(vals[0])
            sigma_xx.append(vals[1])
            sigma_yy.append(vals[5])
            sigma_zz.append(vals[9])

    return np.array(energy), np.array(sigma_xx), np.array(sigma_yy), np.array(sigma_zz)


def read_sp_file(path):
    """
    Read single-particle conductivity file.
    Returns:
      sigma_sp_xx, sigma_sp_yy, sigma_sp_zz
    """
    sigma_sp_xx = []
    sigma_sp_yy = []
    sigma_sp_zz = []
    with open(path, 'r') as f:
        for raw in f:
            vals = [float(v) for v in raw.split()]
            sigma_sp_xx.append(vals[1])
            sigma_sp_yy.append(vals[5])
            sigma_sp_zz.append(vals[9])
    return np.array(sigma_sp_xx), np.array(sigma_sp_yy), np.array(sigma_sp_zz)


def read_oscillator_file(path):
    """
    Read oscillator strengths from third file.
    Returns:
      E_exc, Vx, Vy, Vz
    """
    E_exc = []
    Vx = []
    Vy = []
    Vz = []
    with open(path, 'r') as f:
        for raw in f:
            parts = raw.split()
            if not parts:
                continue
            vals = [float(v) for v in parts]
            E_exc.append(vals[0])
            Vx.append(abs(vals[1] + 1j*vals[2])**2)
            Vy.append(abs(vals[3] + 1j*vals[4])**2)
            Vz.append(abs(vals[5] + 1j*vals[6])**2)
    return np.array(E_exc), np.array(Vx), np.array(Vy), np.array(Vz)


def print_header(plot_type, output_file):
    """Prints a standardized header with current plot info."""
    print("===================================")
    print("       PLOTTING CONDUCTIVITY       ")
    print("===================================")
    print(f"Selected plot type: {plot_type}")
    print(f"Output file: {output_file}")
    print()


def print_usage():
    print("Usage:")
    print("  script.py <exciton_file> <sp_file> [oscillator_file]")
    print()
    print("Running with two files does plot type 1.")
    print("Add a third file to get plot type 2.")


def plot_type1(energy, sigma_xx, sigma_yy, sigma_zz,
               sigma_sp_xx, sigma_sp_yy, sigma_sp_zz,
               output="conductivity.png"):
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    ax.plot(energy, sigma_xx,    "r-",  label="BSE_xx", alpha=0.6)
    ax.plot(energy, sigma_yy,    "orange", label="BSE_yy", alpha=0.6)
    ax.plot(energy, sigma_zz,    "g-",  label="BSE_zz", alpha=0.6)
    ax.plot(energy, sigma_sp_xx, "b--", label="IPA_xx")
    ax.plot(energy, sigma_sp_yy, "c--", label="IPA_yy")
    ax.plot(energy, sigma_sp_zz, ls="--", c="navy", label="IPA_zz")

    ax.set_xlabel("E (eV)", fontsize=18)
    ax.set_ylabel(r"$\sigma$ ($e^2/\hbar$)", fontsize=18)
    ax.legend()
    ax.set_xlim(energy[0], energy[-1])
    plt.tight_layout()
    plt.savefig(output, dpi=600)


def plot_type2(energy, sigma_xx, sigma_yy, sigma_zz,
               sigma_sp_xx, sigma_sp_yy, sigma_sp_zz,
               E_exc, Vx, Vy, Vz,
               output="conductivity_oscillators.png"):
    fig, axes = plt.subplots(2, 1, figsize=(6.4, 4.8))
    ax1 = axes[0]
    ax1.plot(energy, sigma_xx,    "r-",     label="BSE_xx")
    ax1.plot(energy, sigma_yy,    "orange", label="BSE_yy")
    ax1.plot(energy, sigma_zz,    "g-",     label="BSE_zz")
    ax1.plot(energy, sigma_sp_xx, "b--",    label="IPA_xx")
    ax1.plot(energy, sigma_sp_yy, "c--",    label="IPA_yy")
    ax1.plot(energy, sigma_sp_zz, ls="--",  c="navy", label="IPA_zz")
    ax1.set_ylabel(r"$\sigma$ ($e^2/\hbar$)", fontsize=18)

    ax2 = axes[1]

    # mask weak oscilaltors (plotting everything is costly and messy)
    mask_vx = Vx > 0.1
    mask_vy = Vy > 0.1
    mask_vz = Vz > 0.5

    if mask_vx.any():
        ax2.bar(E_exc[mask_vx], E_exc[mask_vx]*Vx[mask_vx], width=0.03, alpha=0.3, label="xx")
    if mask_vy.any():
        ax2.bar(E_exc[mask_vy], E_exc[mask_vy]*Vy[mask_vy], width=0.03, alpha=0.3, label="yy")
    if mask_vz.any():
        ax2.bar(E_exc[mask_vz], E_exc[mask_vz]*Vz[mask_vz], width=0.03, alpha=0.3, label="zz")

    ax2.set_ylabel(r"$|V^\alpha|^2$", fontsize=18)
    ax2.set_xlabel("E (eV)", fontsize=18)

    for ax in axes:
        ax.set_xlim(energy[0], energy[-1])
        ax.tick_params(labelsize=12)

    ax1.legend()
    ax2.legend()
    plt.subplots_adjust(hspace=0.1)
    plt.tight_layout()
    plt.savefig(output, dpi=600)


def main():
    if len(sys.argv) < 3:
        print_usage()
        sys.exit(1)

    exciton_path = sys.argv[1]
    sp_path      = sys.argv[2]
    has_osc      = len(sys.argv) == 4
    plot_type    = 2 if has_osc else 1
    output_file  = "conductivity_oscillators.png" if has_osc else "conductivity.png"

    print_header(plot_type, output_file)
    configure_matplotlib()

    energy, sigma_xx, sigma_yy, sigma_zz = read_exciton_file(exciton_path)
    sigma_sp_xx, sigma_sp_yy, sigma_sp_zz = read_sp_file(sp_path)

    if plot_type == 1:
        plot_type1(energy, sigma_xx, sigma_yy, sigma_zz,
                   sigma_sp_xx, sigma_sp_yy, sigma_sp_zz,
                   output_file)
    else:
        E_exc, Vx, Vy, Vz = read_oscillator_file(sys.argv[3])
        plot_type2(energy, sigma_xx, sigma_yy, sigma_zz,
                   sigma_sp_xx, sigma_sp_yy, sigma_sp_zz,
                   E_exc, Vx, Vy, Vz,
                   output_file)

    # plt.show()

if __name__ == "__main__":
    main()
