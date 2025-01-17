import numpy as np
import matplotlib.pyplot as plt

# Read the data from the file, skipping rows with inconsistent columns
data = np.genfromtxt('hBN_ex.dat', unpack=True, invalid_raise=False)

# Ensure the data has at least 10 columns to avoid index errors
if data.shape[0] != 10:
    raise ValueError(f"Expected 10 columns, but got {data.shape[0]}.")

# Extract the data
omega = data[0]
sigma_columns = data[1:]

# Labels for the sigma components
sigma_labels = [r'σ_xx', r'σ_xy', r'σ_xz', r'σ_yx', r'σ_yy', r'σ_yz', r'σ_zx', r'σ_zy', r'σ_zz']

# Plot all sigma components
plt.figure(figsize=(10, 6))
for i in range(len(sigma_columns)):
    plt.plot(omega, sigma_columns[i], label=sigma_labels[i])

# Plot customization
plt.xlabel('Omega (ω)')
plt.ylabel('Sigma Components')
plt.title('Optical Conductivity Components')
plt.legend()
plt.grid(True)

# Save the figure
plt.savefig('hBN_tb.png')
plt.show()
