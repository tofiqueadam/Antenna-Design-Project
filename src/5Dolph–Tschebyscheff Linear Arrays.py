import numpy as np
import matplotlib.pyplot as plt
from scipy.special import chebyt

# Parameters
N = 10       # Number of elements
d = 0.5      # Element spacing (wavelengths)
theta = np.linspace(-np.pi/2, np.pi/2, 1000)  # Angle range
sll_levels = [20, 30, 40]  # Side-lobe levels (dB)

def dolph_chebyshev_weights(N, sll_db):
    """Calculate Dolph-Chebyshev weights"""
    R = 10**(sll_db/20)
    x0 = np.cosh(np.arccosh(R)/(N-1))
    
    # Array indices
    n = np.arange(N)
    weights = np.zeros(N, dtype=complex)
    
    # Sum over Chebyshev polynomial samples
    for m in range(N):
        theta_m = (2*m + 1)*np.pi/(2*N)
        x_m = x0 * np.cos(theta_m)
        weights += chebyt(N-1)(x_m) * np.exp(1j*2*np.pi*m*n/N)
    
    return np.real(weights)/np.max(np.abs(weights))

def array_factor(weights, theta, d):
    """Compute array factor"""
    N = len(weights)
    n = np.arange(N) - (N-1)/2  # Symmetric indices
    psi = 2*np.pi*d*np.sin(theta)
    AF = np.zeros(len(theta), dtype=complex)
    
    for i, angle in enumerate(psi):
        AF[i] = np.sum(weights * np.exp(1j*n*angle))
    
    return np.abs(AF)

# Compute patterns
plt.figure(figsize=(12, 7))
for sll in sll_levels:
    weights = dolph_chebyshev_weights(N, sll)
    AF = array_factor(weights, theta, d)
    plt.plot(np.degrees(theta), 20*np.log10(AF/np.max(AF)+1e-10), 
             label=f'D-T {sll}dB')

# Uniform array comparison
uniform_weights = np.ones(N)
AF_uniform = array_factor(uniform_weights, theta, d)
plt.plot(np.degrees(theta), 20*np.log10(AF_uniform/np.max(AF_uniform)+1e-10), 
         'k--', label='Uniform')

# Plot formatting
plt.title('Dolph-Tschebyscheff Array Patterns (N=10, d=λ/2)')
plt.xlabel('Angle (degrees)')
plt.ylabel('Normalized Array Factor (dB)')
plt.ylim(-50, 5)
plt.grid(True)
plt.legend()
plt.show()