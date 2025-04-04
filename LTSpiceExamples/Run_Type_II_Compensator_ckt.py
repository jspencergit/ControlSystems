# Run_Type_II_Compensator_ckt.py
from PyLTSpice import SimCommander  # For running LTSpice
import ltspice                      # For reading .raw files
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os
import sys                          # For sys.exit()

# Switch to an interactive backend
matplotlib.use('TkAgg')

# Configuration
LTSPICE_PATH = "C:/Program Files/ADI/LTspice/LTspice.exe"
ASC_FILE = "Type_II_Compensator_ckt.asc"
RAW_FILE = "Type_II_Compensator_ckt_1.raw"  # Matches LTspice output
SCHEMATIC_IMG = "Type_II_Compensator_ckt.jpg"  # Assumed schematic image

# Check if the .asc file exists
if not os.path.exists(ASC_FILE):
    print(f"Error: {ASC_FILE} not found in the current directory!")
    sys.exit(1)

# Modify the .asc file to start AC analysis at 100 Hz
with open(ASC_FILE, 'r') as file:
    asc_content = file.read()
asc_content = asc_content.replace('.ac oct 1000 10 1Meg', '.ac oct 1000 100 1Meg')
with open(ASC_FILE, 'w') as file:
    file.write(asc_content)

# Run the simulation
print("Starting LTSpice simulation...")
lt = SimCommander(ASC_FILE)
lt.run()
if not lt.wait_completion(timeout=30):  # Wait up to 30 seconds
    print("Simulation timed out!")
    sys.exit(1)
print("Simulation complete. Output saved as", RAW_FILE)

# Load the raw file
if not os.path.exists(RAW_FILE):
    print(f"Error: {RAW_FILE} not found! Simulation may have failed.")
    sys.exit(1)

l = ltspice.Ltspice(RAW_FILE)
l.parse()

# Extract frequency, magnitude, and phase data
freq = l.get_frequency()
vout = l.get_data('V(out)')  # Matches FLAG 576 96 out
mag = 20 * np.log10(np.abs(vout))  # Magnitude in dB
phase = np.angle(vout, deg=True)   # Phase in degrees

# Calculate pole and zero frequencies from component values
R1 = 500e3    # 500 kΩ
R2 = 10e3     # 10 kΩ
C1 = 1000e-12 # 1000 pF
C2 = 100e-12  # 100 pF
fp1 = 1 / (2 * np.pi * (R1 + R2) * C1)  # Pole 1 frequency ~312 Hz
fz = 1 / (2 * np.pi * R2 * C1)          # Zero frequency ~15.9 kHz
fp2 = 1 / (2 * np.pi * R2 * C2)         # Pole 2 frequency ~159 kHz

# Calculate frequency of maximum phase (between fz and fp2)
f_max_phase = np.sqrt(fz * fp2)  # Geometric mean ~50.3 kHz

# Find indices for pole, zero, and max phase frequencies
fp1_idx = np.argmin(np.abs(freq - fp1))
fz_idx = np.argmin(np.abs(freq - fz))
fp2_idx = np.argmin(np.abs(freq - fp2))
f_max_idx = np.argmin(np.abs(freq - f_max_phase))

# Calculate phase boost as the difference from -90°
phase_at_max = phase[f_max_idx]  # Phase at 50.3 kHz
phase_boost_deg = phase_at_max - (-90)  # e.g., -33 - (-90) = 57°

# Create figure with three subplots: magnitude, phase, and schematic
fig = plt.figure(figsize=(12, 12))
gs = fig.add_gridspec(3, 1, height_ratios=[2, 2, 1])  # 2:2:1 ratio for heights

# Magnitude plot (first subplot)
ax1 = fig.add_subplot(gs[0, 0])
mag_plot, = ax1.semilogx(freq, mag, label='Magnitude (dB)', color='blue')
ax1.grid(True, which="both", ls="--")
ax1.set_ylabel('Magnitude (dB)')
ax1.set_title('Type II Compensator Frequency Response: Gain and Phase')

# Markers for poles and zero on magnitude plot
ax1.axvline(x=fp1, color='red', linestyle='--', label=f'Pole 1 at {fp1:.0f} Hz')
ax1.plot(fp1, mag[fp1_idx], 'rx', markersize=10)
ax1.axvline(x=fz, color='green', linestyle='--', label=f'Zero at {fz/1000:.1f} kHz')
ax1.plot(fz, mag[fz_idx], 'go', markersize=10)
ax1.axvline(x=fp2, color='purple', linestyle='--', label=f'Pole 2 at {fp2/1000:.0f} kHz')
ax1.plot(fp2, mag[fp2_idx], 'mo', markersize=10)
ax1.legend(loc='upper right')  # Explicitly set legend position

# Phase plot (second subplot)
ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
phase_plot, = ax2.semilogx(freq, phase, label='Phase (deg)', color='purple')
ax2.grid(True, which="both", ls="--")
ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Phase (deg)')
ax2.set_ylim(-180, 90)  # Set max y-axis to 90°

# Markers for poles, zero, and phase boost on phase plot
ax2.axvline(x=fp1, color='red', linestyle='--')
ax2.plot(fp1, phase[fp1_idx], 'rx', markersize=10)
ax2.axvline(x=fz, color='green', linestyle='--')
ax2.plot(fz, phase[fz_idx], 'go', markersize=10)
ax2.axvline(x=fp2, color='purple', linestyle='--')
ax2.plot(fp2, phase[fp2_idx], 'mo', markersize=10)
ax2.axvline(x=f_max_phase, color='orange', linestyle='--', 
            label=f'Phase Boost at {f_max_phase/1000:.1f} kHz ({phase_boost_deg:.0f}°)')
ax2.plot(f_max_phase, phase[f_max_idx], 'yo', markersize=10)
ax2.legend(loc='upper right')  # Explicitly set legend position

# Schematic plot (third subplot)
ax3 = fig.add_subplot(gs[2, 0])
if not os.path.exists(SCHEMATIC_IMG):
    ax3.text(0.5, 0.5, "Schematic image not found", ha="center", va="center")
    ax3.axis('off')
else:
    ax3.imshow(plt.imread(SCHEMATIC_IMG))
    ax3.axis('off')
ax3.set_title('Type II Compensator Schematic')

# Custom interactive cursor for magnitude plot
mag_dot, = ax1.plot([], [], 'ko', markersize=5)  # Black dot for magnitude plot
mag_text = ax1.text(0.5, 0.85, '', transform=ax1.transAxes, ha='center', va='top', fontsize=10)

def on_move_mag(event):
    if event.inaxes != ax1:
        mag_dot.set_data([], [])
        mag_text.set_text('')
        fig.canvas.draw_idle()
        return
    x = event.xdata
    if x is None:
        return
    idx = np.argmin(np.abs(freq - x))
    mag_dot.set_data([freq[idx]], [mag[idx]])
    mag_text.set_text(f'Freq: {freq[idx]:.1e} Hz\nValue: {mag[idx]:.2f} dB')
    fig.canvas.draw_idle()

# Custom interactive cursor for phase plot
phase_dot, = ax2.plot([], [], 'ko', markersize=5)  # Black dot for phase plot
phase_text = ax2.text(0.5, 0.85, '', transform=ax2.transAxes, ha='center', va='top', fontsize=10)

def on_move_phase(event):
    if event.inaxes != ax2:
        phase_dot.set_data([], [])
        phase_text.set_text('')
        fig.canvas.draw_idle()
        return
    x = event.xdata
    if x is None:
        return
    idx = np.argmin(np.abs(freq - x))  # Fixed: Use x instead of fz
    phase_dot.set_data([freq[idx]], [phase[idx]])
    phase_text.set_text(f'Freq: {freq[idx]:.1e} Hz\nValue: {phase[idx]:.2f}°')
    fig.canvas.draw_idle()

# Connect the event handlers
fig.canvas.mpl_connect('motion_notify_event', on_move_mag)
fig.canvas.mpl_connect('motion_notify_event', on_move_phase)

# Adjust layout and show
plt.tight_layout()
plt.show()

# Run the cleanup script
if os.path.exists("LTSpice_Cleanup.py"):
    os.system("python LTSpice_Cleanup.py")
else:
    print("Warning: Cleanup script not found.")