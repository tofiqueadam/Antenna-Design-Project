import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

# Try to set an interactive backend for displaying plots
try:
    plt.switch_backend('TkAgg')  #  TkAgg for interactive rendering
except ImportError:
    print("Warning: TkAgg backend unavailable. Install Tkinter (e.g., pip install tk). Falling back to default backend.")

# Constants
lambda_val = 1.0  # Normalized wavelength (m)
k = 2 * np.pi / lambda_val  # Wave number
spacing_start = 0.1 * lambda_val
spacing_end = 2.0 * lambda_val
num_spacing_points = 100
spacings = np.linspace(spacing_start, spacing_end, num_spacing_points)
num_elements_list = [2, 4, 8, 16]  # Array sizes to analyze

# Function to calculate array factor
def array_factor(theta, n_elements, d, k):
    """
    Computes the radiation pattern of a uniform linear dipole array.
    Args:
        theta: Elevation angle in radians
        n_elements: Number of dipole elements
        d: Inter-element spacing (m)
        k: Wavenumber (2π/λ)
    Returns:
        float: Power pattern of the array (|AF|²)
    """
    if np.sin(theta) == 0 and d == 0:  # Avoid division by zero
        return n_elements**2
    kd = k * d
    af = np.sin(n_elements * kd * np.cos(theta) / 2) / np.sin(kd * np.cos(theta) / 2)
    return af**2

# Function to calculate directivity
def array_directivity(n_elements, d_lambda):
    """
    Computes the maximum directivity of a dipole array.
    Args:
        n_elements: Number of radiating elements
        d_lambda: Element spacing in wavelengths
    Returns:
        float: Peak directivity (unitless)
    """
    d = d_lambda * lambda_val
    # Numerical integration of radiation pattern
    integrand = lambda theta: array_factor(theta, n_elements, d, k) * np.sin(theta)
    integral, _ = quad(integrand, 0, np.pi, epsabs=1e-8)
    if integral == 0:  # Prevent division by zero
        return np.nan
    # Directivity formula for broadside linear arrays
    directivity = 2 * n_elements**2 / integral
    return directivity

# Calculate directivity for each array size and spacing
directivity_data = {}
peak_directivities = {}
for n_elements in num_elements_list:
    directivity_values = []
    for d in spacings:
        d_lambda = d / lambda_val
        directivity = array_directivity(n_elements, d_lambda)
        directivity_values.append(directivity)
    directivity_data[n_elements] = directivity_values
    # Find peak directivity for this array size
    valid_directivities = [d for d in directivity_values if not np.isnan(d)]
    peak_directivities[n_elements] = max(valid_directivities) if valid_directivities else np.nan
    
# Plotting the directivity curves
plt.figure(figsize=(12, 7))
for n_elements, directivity_values in directivity_data.items():
    plt.plot(spacings / lambda_val, directivity_values, 
             label=f'{n_elements} Elements', linewidth=2)
plt.axvline(x=1.0, color='gray', linestyle='--', 
            alpha=0.5, label='Grating Lobe Threshold (d = λ)')
plt.xlabel('Element Spacing (wavelengths)')
plt.ylabel('Peak Directivity (dB)')
plt.title('Dipole Array Performance: Directivity vs. Element Spacing')
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()
plt.tight_layout()

# Save plot
plt.savefig('dipole_array_directivity_analysis.png')

# Display plot
plt.show(block=True)

# This part of code generates the report of the above code.

