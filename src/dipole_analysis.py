import numpy as np
import matplotlib.pyplot as plt
from scipy.special import sici
from scipy.constants import pi

# --------------------------
# Parameter Selection
# --------------------------
# Lengths of dipole in terms of wavelength (L/λ) from 0.1 to 2.5 in steps of 0.01
lengths = np.linspace(0.1, 2.5, 241)

# --------------------------
# Analytical Functions
# --------------------------
def dipole_impedance(L):
    """
    Calculate input impedance (resistance and reactance) for a center-fed dipole.
    Uses approximated analytical expressions based on Balanis' formulation.

    Parameters:
        L : float
            Length of dipole in wavelengths.

    Returns:
        R : float
            Input resistance in ohms.
        X : float
            Input reactance in ohms.
    """
    kl = 2 * pi * L  # Electrical length (k*l)
    Cin, _ = sici(2 * kl)     # Cosine integral at 2kl
    Si_kl, Ci_kl = sici(kl)   # Sine and cosine integral at kl

    # Real part of impedance (resistance)
    R = 30 * (
        0.5772 + np.log(kl) - Cin
        + 0.5 * np.sin(kl) * (Si_kl - 0.5 * Si_kl)
        - 0.5 * np.cos(kl) * (0.5772 + np.log(kl / 2) + Cin - 2 * Si_kl)
    )

    # Imaginary part of impedance (reactance)
    X = 30 * (
        2 * Si_kl
        + np.cos(kl) * (2 * Si_kl - Si_kl)
        - np.sin(kl) * (2 * Ci_kl + 0.5772 + np.log(kl / 2) - Ci_kl)
    )
    return R, X

def dipole_directivity(L):
    """
    Compute approximate maximum directivity of a dipole.

    Parameters:
        L : float
            Dipole length in wavelengths.

    Returns:
        D : float
            Directivity in dB.
    """
    if L < 0.637:
        return 10 * np.log10(1.5)  # Short dipole constant directivity ~1.76 dB
    else:
        return 10 * np.log10(2.41 * L / (0.5 * (1 - np.exp(-1.386 * L))))

# --------------------------
# Main Computation
# --------------------------
R = np.zeros_like(lengths)  # Resistance array
X = np.zeros_like(lengths)  # Reactance array
D = np.zeros_like(lengths)  # Directivity array

# Loop over each length and compute impedance and directivity
for i, L in enumerate(lengths):
    R[i], X[i] = dipole_impedance(L)
    D[i] = dipole_directivity(L)

# --------------------------
# Plotting Results
# --------------------------
# Figure 1: Input Impedance vs Length
plt.figure(figsize=(8, 5))
plt.plot(lengths, R, label='R_in (Ω)', color='blue')
plt.plot(lengths, X, '--', label='X_in (Ω)', color='orange')
plt.axhline(0, color='black', linestyle=':')

# Mark resonance points (where reactance crosses zero)
zeros = np.where(np.diff(np.sign(X)))[0]
for z in zeros:
    if abs(X[z]) < 10:  # Filter near-zero crossings
        plt.plot(lengths[z], 0, 'ro')
        plt.annotate(f"{lengths[z]:.2f}λ", (lengths[z], 0), textcoords='offset points', xytext=(5, 5))

plt.xlabel('L / λ')
plt.ylabel('Impedance (Ω)')
plt.title('Figure 1: Dipole Input Impedance vs Length')
plt.grid(True)
plt.legend()
plt.tight_layout()

# Figure 2: Directivity vs Length
plt.figure(figsize=(8, 5))
plt.plot(lengths, D, label='Directivity', color='green')

# Mark peaks in directivity
peaks = np.r_[True, D[1:] > D[:-1]] & np.r_[D[:-1] > D[1:], True]
for i, p in enumerate(peaks):
    if p and lengths[i] > 0.2:
        plt.plot(lengths[i], D[i], 'ro')
        plt.annotate(f"{D[i]:.1f} dB\n({lengths[i]:.2f}λ)", (lengths[i], D[i]), textcoords='offset points', xytext=(5, 5))

plt.xlabel('L / λ')
plt.ylabel('Directivity (dB)')
plt.title('Figure 2: Dipole Directivity vs Length')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()