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
ASC_FILE = "PoleZeroPair_Circuit.asc"
RAW_FILE = "PoleZeroPair_Circuit_1.raw"    # Updated to match LTspice output

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
vout = l.get_data('V(out)')  # Matches FLAG 448 96 out
mag = 20 * np.log10(np.abs(vout))  # Magnitude in dB
phase = np.angle(vout, deg=True)   # Phase in degrees

# Calculate pole and zero frequencies from component values
R1 = 500e3    # 500kΩ
R2 = 10e3     # 10kΩ
C1 = 1000e-12 # 1000pF
fp = 1/(2*np.pi*(R1 + R2)*C1)  # Pole frequency ~312Hz
fz = 1/(2*np.pi*R2*C1)         # Zero frequency ~15.9kHz

# Calculate frequency and value of maximum phase deflection
f_max_phase = np.sqrt(fp * fz)  # Geometric mean ~2.2kHz
max_phase = np.arctan(np.sqrt(fp/fz)) - np.arctan(np.sqrt(fz/fp))  # ~74°
max_phase_deg = np.degrees(max_phase)  # Convert to degrees

# Find indices for pole, zero, and max phase frequencies
fp_idx = np.argmin(np.abs(freq - fp))
fz_idx = np.argmin(np.abs(freq - fz))
f_max_idx = np.argmin(np.abs(freq - f_max_phase))

# Create figure with three subplots: magnitude, phase, and schematic
fig = plt.figure(figsize=(12, 12))  # Increased height to accommodate third subplot
# Define subplot layout: 3 rows, 1 column
# Adjust heights to give more space to magnitude and phase plots
gs = fig.add_gridspec(3, 1, height_ratios=[2, 2, 1])  # 2:2:1 ratio for heights

# Magnitude plot (first subplot)
ax1 = fig.add_subplot(gs[0, 0])
mag_plot, = ax1.semilogx(freq, mag, label='Magnitude (dB)', color='blue')
ax1.grid(True, which="both", ls="--")
ax1.set_ylabel('Magnitude (dB)')
ax1.set_title('Frequency Response: Gain and Phase')

# Markers for pole and zero on magnitude plot
ax1.axvline(x=fp, color='red', linestyle='--', label=f'Pole at {fp:.0f}Hz')
ax1.plot(fp, mag[fp_idx], 'rx', markersize=10)
ax1.axvline(x=fz, color='green', linestyle='--', label=f'Zero at {fz/1000:.1f}kHz')
ax1.plot(fz, mag[fz_idx], 'go', markersize=10)
ax1.legend()

# Phase plot (second subplot)
ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
phase_plot, = ax2.semilogx(freq, phase, label='Phase (deg)', color='purple')
ax2.grid(True, which="both", ls="--")
ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Phase (deg)')
ax2.set_ylim(-180, 180)  # Limit phase to ±180°

# Markers for pole, zero, and max phase on phase plot
ax2.axvline(x=fp, color='red', linestyle='--')
ax2.plot(fp, phase[fp_idx], 'rx', markersize=10)
ax2.axvline(x=fz, color='green', linestyle='--')
ax2.plot(fz, phase[fz_idx], 'go', markersize=10)
ax2.axvline(x=f_max_phase, color='orange', linestyle='--', 
            label=f'Max Phase at {f_max_phase/1000:.1f}kHz ({max_phase_deg:.0f}°)')
ax2.plot(f_max_phase, phase[f_max_idx], 'yo', markersize=10)
ax2.legend()

# Schematic plot (third subplot)
ax3 = fig.add_subplot(gs[2, 0])
schematic_img = plt.imread("PoleZeroPair_Circuit.jpg")
ax3.imshow(schematic_img)
ax3.axis('off')  # Hide axes for the schematic
ax3.set_title('Schematic')

# Add interactive cursor to magnitude and phase plots
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