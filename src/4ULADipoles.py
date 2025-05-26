import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

# Try to set an interactive backend for displaying plots
try:
    plt.switch_backend('TkAgg')  # Use TkAgg for interactive rendering
except ImportError:
    print("Warning: TkAgg backend unavailable. Install Tkinter (e.g., pip install tk). Falling back to default backend.")


    - Theoretical maximum directivity approaches 2N (e.g., {2*2:.1f}, {2*4:.1f}, {2*8:.1f}, {2*16:.1f} for N = 2, 4, 8, 16) at optimal spacing.
- **Element Spacing (d)**:
    - At small spacings (d < 0.5 lambda), mutual coupling reduces directivity by limiting the effective array aperture.
    - Optimal directivity occurs around d approximately 0.5–0.8 lambda, balancing beamwidth and gain without grating lobes.
    - For d > lambda, grating lobes emerge, distributing power to secondary lobes and reducing main lobe directivity.
    - The grating lobe onset at d = lambda is marked on the plot, showing a drop in directivity beyond this point.
- **Trends**:
    - Larger arrays (e.g., N = 16) exhibit sharper directivity peaks, indicating higher sensitivity to spacing variations.
    - Smaller arrays (e.g., N = 2) show flatter directivity curves, less affected by spacing changes.
- **Practical Implications**:
    - For high directivity, use d approximately 0.5–0.8 lambda with larger N, avoiding grating lobes.
    - Spacings d > lambda are unsuitable for applications needing a single main lobe (e.g., radar, communications).
- **Plot Details**:
    - Directivity curves are plotted for N = 2, 4, 8, 16 over spacings from 0.1 to 2.0 lambda.
    - Plot saved as 'dipole_array_directivity_analysis.png' and displayed.
"""
print(summary)
