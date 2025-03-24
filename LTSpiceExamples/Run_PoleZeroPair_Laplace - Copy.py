# simulate_and_plot.py
from PyLTSpice import SimCommander  # For running LTSpice
import ltspice                      # For reading .raw files
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os
import sys                          # Added for sys.exit()
import mplcursors                   # For interactive cursor

# Switch to an interactive backend
matplotlib.use('Qt5Agg')  # Separate window

# Configuration
LTSPICE_PATH = "C:/Program Files/ADI/LTspice/LTspice.exe"
ASC_FILE = "PoleZeroPair_Laplace.asc"
RAW_FILE = "PoleZeroPair_Laplace_1.raw"    # Updated to match LTspice output

# Check if the .asc file exists
if not os.path.exists(ASC_FILE):
    print(f"Error: {ASC_FILE} not found in the current directory!")
    sys.exit(1)                    # Use sys.exit() instead of exit()

# Run the simulation
print("Starting LTSpice simulation...")
lt = SimCommander(ASC_FILE)
lt.run()
lt.wait_completion()
print("Simulation complete. Output saved as PoleZeroPair_1.raw")

# Load the raw file
if not os.path.exists(RAW_FILE):
    print(f"Error: {RAW_FILE} not found! Simulation may have failed.")
    sys.exit(1)                    # Use sys.exit() instead of exit()

l = ltspice.Ltspice(RAW_FILE)
l.parse()

# Extract frequency, magnitude, and phase data
freq = l.get_frequency()
vout = l.get_data('V(out)')  # Matches FLAG 288 96 OUT
mag = 20 * np.log10(np.abs(vout))  # Magnitude in dB
phase = np.angle(vout, deg=True)   # Phase in degrees

# Create subplots: magnitude and phase
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

# Plot Magnitude
mag_plot, = ax1.semilogx(freq, mag, label='Magnitude (dB)', color='blue')
ax1.grid(True, which="both", ls="--")
ax1.set_ylabel('Magnitude (dB)')
ax1.set_title('Frequency Response: Gain and Phase')

# Markers for poles and zeros on magnitude plot
ax1.axvline(x=10, color='red', linestyle='--', label='Pole at 0 Hz')
ax1.plot(10, mag[0], 'rx', markersize=10, label='_')
zero_freq = 10e3
zero_idx = np.argmin(np.abs(freq - zero_freq))  # Fixed typo: freq - zero_freq
ax1.axvline(x=zero_freq, color='green', linestyle='--', label='Zero at 10kHz')
ax1.plot(zero_freq, mag[zero_idx], 'go', markersize=10, label='_')
pole_freq = 100e3
pole_idx = np.argmin(np.abs(freq - pole_freq))
ax1.axvline(x=pole_freq, color='red', linestyle='--', label='Pole at 100kHz')
ax1.plot(pole_freq, mag[pole_idx], 'rx', markersize=10, label='_')
ax1.legend()

# Plot Phase
phase_plot, = ax2.semilogx(freq, phase, label='Phase (deg)', color='purple')
ax2.grid(True, which="both", ls="--")
ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Phase (deg)')
ax2.set_ylim(-180, 180)  # Limit phase to ±180°

# Markers for poles and zeros on phase plot
ax2.axvline(x=10, color='red', linestyle='--', label='_')
ax2.plot(10, phase[0], 'rx', markersize=10, label='_')
ax2.axvline(x=zero_freq, color='green', linestyle='--', label='_')
ax2.plot(zero_freq, phase[zero_idx], 'go', markersize=10, label='_')
ax2.axvline(x=pole_freq, color='red', linestyle='--', label='_')
ax2.plot(pole_freq, phase[pole_idx], 'rx', markersize=10, label='_')

# Add interactive cursor
cursor = mplcursors.cursor([mag_plot, phase_plot], hover=True)
@cursor.connect("add")
def on_add(sel):
    x, y = sel.target
    sel.annotation.set_text(f'Freq: {x:.1e} Hz\nValue: {y:.2f}')

# Adjust layout and show
plt.tight_layout()
plt.show()

# Run the cleanup script
os.system("python LTSpice_Cleanup.py")