import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Generate spherical coordinates
theta = np.linspace(0, np.pi, 180)      # Elevation angle
phi = np.linspace(0, 2 * np.pi, 360)    # Azimuth angle
theta_2d = np.linspace(0, 2 * np.pi, 360)  # For 2D pattern
theta_grid, phi_grid = np.meshgrid(theta, phi)

# Placeholder gain pattern
r_3d = 10 * np.abs(np.cos(theta_grid))**3
r_2d = 10 * np.abs(np.cos(theta_2d))**3

# Convert spherical to Cartesian coordinates for 3D plot
x = r_3d * np.sin(theta_grid) * np.cos(phi_grid)
y = r_3d * np.sin(theta_grid) * np.sin(phi_grid)
z = r_3d * np.cos(theta_grid)

# Create 2D and 3D plots
fig = plt.figure(figsize=(14, 6))

# 2D Polar Plot
ax1 = fig.add_subplot(1, 2, 1, polar=True)
ax1.plot(theta_2d, r_2d, linewidth=2)
ax1.set_title("2D Radiation Pattern (Polar)")

# 3D Radiation Pattern
ax2 = fig.add_subplot(1, 2, 2, projection='3d')
ax2.plot_surface(x, y, z, cmap='viridis', alpha=0.9, edgecolor='none')
ax2.set_title("3D Radiation Pattern (Placeholder)")
ax2.set_xlabel("X")
ax2.set_ylabel("Y")
ax2.set_zlabel("Z")

plt.tight_layout()
plt.show()
