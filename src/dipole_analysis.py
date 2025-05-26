import numpy as np
import matplotlib.pyplot as plt
from scipy.special import sici
from scipy.constants import pi

#--------------------------
# Parameter Selection
#--------------------------
# Define a range of dipole lengths (L) in wavelengths (λ)
# Sweep from 0.1λ to 2.5λ in 241 steps (increment of 0.01λ)
lengths = np.linspace(0.1, 2.5, 241)  # L/λ sweep

#--------------------------
# Analytical Functions
#--------------------------
def dipole_impedance(L):
    """
    Compute the input resistance (R) and reactance (X)
    for a center-fed dipole of length L (in wavelengths).
    Uses analytical approximations involving sine and cosine integrals.
    
    Parameters:
      L : float
          Dipole length in wavelengths (λ).
    
    Returns:
      R : float
          Input resistance in ohms (Ω).
      X : float
          Input reactance in ohms (Ω).
    """
    # Electrical length: k * L, where k = 2π/λ and L is in wavelengths
    kl = 2 * pi * L

    # Compute cosine integral at 2*kl (Ci(2kl)) and sine integral at kl (Si(kl))
    Cin, _   = sici(2 * kl)   # Ci(2kl), ignore second return (Si)
    Si2kl, _ = sici(kl)       # Si(kl), ignore second return (Ci)

    # Calculate input resistance R using:
    # R = 30 [ γ + ln(kl) - Ci(2kl)
    #          + 0.5·sin(kl)·(Si(kl) - 0.5·Si(kl))
    #          - 0.5·cos(kl)·(γ + ln(kl/2) + Ci(2kl) - 2·Si(kl)) ]
    # where γ ≈ 0.5772 is Euler's constant
    R = 30 * (
        0.5772
        + np.log(kl)             # Natural log of kl
        - Cin                    # Subtract Ci(2kl)
        + 0.5 * np.sin(kl) * (Si2kl - 0.5 * Si2kl)  # Sin term
        - 0.5 * np.cos(kl) *     # Cos term
          (0.5772 + np.log(kl/2) + Cin - 2 * Si2kl)
    )

    # Compute cosine integral at kl for reactance term
    _, Ci_kl = sici(kl)        # Ci(kl), ignore first return (Si)

    # Calculate input reactance X using:
    # X = 30 [ 2·Si(kl)
    #          + cos(kl)·(2·Si(kl) - Si(kl))
    #          - sin(kl)·(2·Ci(kl) + γ + ln(kl/2) - Ci(kl)) ]
    X = 30 * (
        2 * Si2kl                              # 2·Si(kl)
        + np.cos(kl) * (2 * Si2kl - Si2kl)     # Cosine component
        - np.sin(kl) *                         # Sine component
          (2 * Ci_kl + 0.5772 + np.log(kl/2) - Ci_kl)
    )

    return R, X

def dipole_directivity(L):
    """
    Compute the maximum directivity D (in dB) for a dipole
    of length L (in wavelengths).

    Uses:
      - Short dipole approximation (L < 0.637λ): D ≈ 1.5 linear → 1.76 dB.
      - Empirical formula for longer dipoles.
    
    Parameters:
      L : float
          Dipole length in wavelengths (λ).
    
    Returns:
      D : float
          Maximum directivity in decibels (dB).
    """
    if L < 0.637:
        # Short dipole region
        return 10 * np.log10(1.5)  # Convert linear directivity to dB
    else:
        # Empirical directivity formula for longer dipoles
        return 10 * np.log10(
            2.41 * L / (0.5 * (1 - np.exp(-1.386 * L)))
        )

#--------------------------
# Main Computation
#--------------------------
# Preallocate arrays to store R, X, and D values for each length
R = np.zeros_like(lengths)
X = np.zeros_like(lengths)
D = np.zeros_like(lengths)

# Loop over each dipole length, compute impedance and directivity
for i, L in enumerate(lengths):
    R[i], X[i] = dipole_impedance(L)
    D[i]       = dipole_directivity(L)

#--------------------------
# Plotting the Results
#--------------------------

# Figure 1: Input Impedance vs. Length
plt.figure(figsize=(8, 5))
plt.plot(lengths, R, label='R_in (Resistance)')
plt.plot(lengths, X, '--', label='X_in (Reactance)')
plt.axhline(0, color='k', linestyle=':', linewidth=0.8)  # Zero reactance line
plt.xlabel('Length (L/λ)')
plt.ylabel('Impedance (Ω)')
plt.title('Figure 1: Dipole Input Impedance')
plt.grid(True)
plt.legend()

# Mark resonant points where reactance crosses zero
zero_crossings = np.where(np.diff(np.sign(X)))[0]
for z in zero_crossings:
    if abs(X[z]) < 10:  # Only mark significant zero crossings
        L_res = lengths[z]
        plt.plot(L_res, 0, 'ro')  # Highlight resonance marker
        plt.annotate(f'{L_res:.2f}λ',
                     (L_res, 0),
                     textcoords='offset points',
                     xytext=(5, 5))

# Figure 2: Directivity vs. Length
plt.figure(figsize=(8, 5))
plt.plot(lengths, D, label='D_max')
plt.xlabel('Length (L/λ)')
plt.ylabel('Directivity (dB)')
plt.title('Figure 2: Dipole Maximum Directivity')
plt.grid(True)

# Mark directivity peaks for clarity
peaks = np.r_[True, D[1:] > D[:-1]] & np.r_[D[:-1] > D[1:], True]
for i, is_peak in enumerate(peaks):
    if is_peak and lengths[i] > 0.2:
        plt.plot(lengths[i], D[i], 'bo')
        plt.annotate(f'{D[i]:.1f} dB\n({lengths[i]:.2f}λ)',
                     (lengths[i], D[i]),
                     textcoords='offset points',
                     xytext=(5, 5))

plt.show()
