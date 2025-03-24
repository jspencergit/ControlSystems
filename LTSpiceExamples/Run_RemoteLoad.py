# Run_RemoteLoad.py
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
ASC_FILE = "Classical_Feedback_Example_FrequencyResponse_RemoteLoad.asc"
RAW_FILE_BASE = "Classical_Feedback_Example_FrequencyResponse_RemoteLoad"

# List of R_ISO values to simulate (in ohms)
R_ISO_VALUES = [1e-3, 1e6]  # 1mΩ, 1MΩ
COLORS = ['blue', 'red']  # Colors for each R_ISO value

# Check if the .asc file exists
if not os.path.exists(ASC_FILE):
    print(f"Error: {ASC_FILE} not found in the current directory!")
    sys.exit(1)

# Initialize SimCommander
lt = SimCommander(ASC_FILE)

# Create subplots: magnitude and phase
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

# Lists to store bandwidth and phase margin for display
bandwidths = []
phase_margins = []

# Loop over R_ISO values to run simulations and plot
for idx, r_iso in enumerate(R_ISO_VALUES):
    # Set the parameter R_ISO
    lt.set_parameter('R_ISO', r_iso)
    
    # Run the simulation
    print(f"Starting LTSpice simulation for R_ISO = {r_iso}...")
    lt.run()
    lt.wait_completion()
    print(f"Simulation complete for R_ISO = {r_iso}.")
    
    # Load the raw file (PyLTSpice may append a number to the raw file)
    raw_file = f"{RAW_FILE_BASE}_{idx+1}.raw"
    if not os.path.exists(raw_file):
        print(f"Error: {raw_file} not found! Simulation may have failed.")
        sys.exit(1)
    
    # Parse the raw file
    l = ltspice.Ltspice(raw_file)
    l.parse()
    
    # Extract frequency and node voltages
    freq = l.get_frequency()
    v_remote_sense = l.get_data('V(remote_sense)')
    v_fb = l.get_data('V(fb)')
    
    # Compute the transfer function V(Remote_Sense)/V(FB)
    transfer = v_remote_sense / v_fb
    mag = 20 * np.log10(np.abs(transfer))  # Magnitude in dB
    phase = np.angle(transfer, deg=True)   # Phase in degrees
    
    # Calculate Bandwidth (frequency where magnitude crosses 0dB)
    gain_crossover_idx = np.where(mag <= 0)[0]
    if len(gain_crossover_idx) > 0:
        bandwidth_freq = freq[gain_crossover_idx[0]]
        phase_at_crossover = phase[gain_crossover_idx[0]]
    else:
        bandwidth_freq = freq[-1]  # Fallback if no crossing
        phase_at_crossover = phase[-1]
    bandwidths.append(bandwidth_freq)
    
    # Calculate Phase Margin (phase at gain crossover - 0°)
    phase_margin = phase_at_crossover - 0
    phase_margins.append(phase_margin)
    
    # Plot Magnitude
    mag_plot, = ax1.semilogx(freq, mag, label=f'R_ISO = {r_iso if r_iso >= 1 else r_iso*1000:.0f}{"MΩ" if r_iso >= 1 else "mΩ"}', color=COLORS[idx])
    ax1.grid(True, which="both", ls="--")
    ax1.set_ylabel('Magnitude (dB)')
    ax1.set_title('Frequency Response: V(Remote_Sense)/V(FB)')
    
    # Add bandwidth marker on magnitude plot
    ax1.axvline(x=bandwidth_freq, color=COLORS[idx], linestyle='--', alpha=0.5)
    
    # Plot Phase
    phase_plot, = ax2.semilogx(freq, phase, label=f'R_ISO = {r_iso if r_iso >= 1 else r_iso*1000:.0f}{"MΩ" if r_iso >= 1 else "mΩ"}', color=COLORS[idx])
    ax2.grid(True, which="both", ls="--")
    ax2.set_xlabel('Frequency (Hz)')
    ax2.set_ylabel('Phase (deg)')
    ax2.set_ylim(-180, 180)  # Limit phase to ±180°
    
    # Add gain crossover marker on phase plot
    ax2.axvline(x=bandwidth_freq, color=COLORS[idx], linestyle='--', alpha=0.5)

# Add legends to the plots
# Move magnitude legend to upper-right corner
ax1.legend(loc='upper right')
# Phase legend remains in lower-left corner
ax2.legend(loc='lower left')

# Add text box with bandwidth and phase margin information in lower-left corner
textstr = '\n'.join([
    f'R_ISO = 1mΩ:',
    f'  Bandwidth: {bandwidths[0]/1000:.2f} kHz',
    f'  Phase Margin: {phase_margins[0]:.1f}°',
    '',  # Add an extra newline for spacing
    f'R_ISO = 1MΩ:',
    f'  Bandwidth: {bandwidths[1]/1000:.2f} kHz',
    f'  Phase Margin: {phase_margins[1]:.1f}°'
])
props = dict(boxstyle='round', facecolor='white', alpha=0.8)
ax1.text(0.05, 0.05, textstr, transform=ax1.transAxes, fontsize=10,
         verticalalignment='bottom', horizontalalignment='left', bbox=props)

# Add interactive cursor
cursor = mplcursors.cursor([mag_plot, phase_plot], hover=True)
@cursor.connect("add")
def on_add(sel):
    x, y = sel.target
    sel.annotation.set_text(f'Freq: {x:.1e} Hz\nValue: {y:.2f}')

# Adjust layout and show
plt.tight_layout()
plt.show()

# Run the cleanup script (if you have one)
os.system("python LTSpice_Cleanup.py")