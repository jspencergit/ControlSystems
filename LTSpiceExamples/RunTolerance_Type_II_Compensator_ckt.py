# RunTolerance_Type_II_Compensator_ckt.py
from PyLTSpice import SimCommander  # For running LTSpice
import ltspice                      # For reading .raw files
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os
import sys                          # For sys.exit()
import re                           # For regular expressions

# Switch to an interactive backend
matplotlib.use('TkAgg')

# Configuration
LTSPICE_PATH = "C:/Program Files/ADI/LTspice/LTspice.exe"
ASC_FILE = "Type_II_Compensator_ckt.asc"
NET_FILE = "Type_II_Compensator_ckt.net"  # Netlist file
RAW_FILE_BASE = "Type_II_Compensator_ckt"  # Base name for .raw files
SCHEMATIC_IMG = "Type_II_Compensator_ckt.jpg"  # Assumed schematic image

# Parameter tolerance
R1_NOMINAL = 500e3  # 500 kΩ
R1_TOLERANCE = 0.10  # 10% tolerance
R1_MIN = R1_NOMINAL * (1 - R1_TOLERANCE)  # 450 kΩ
R1_MAX = R1_NOMINAL * (1 + R1_TOLERANCE)  # 550 kΩ
R1_STEPS = 11  # Number of steps for R1 variation (including min and max)

# Check if the .asc file exists
if not os.path.exists(ASC_FILE):
    print(f"Error: {ASC_FILE} not found in the current directory!")
    sys.exit(1)

# Clean up old .raw files to avoid confusion
for i in range(1, R1_STEPS + 1):
    raw_file = f"{RAW_FILE_BASE}_{i}.raw"
    if os.path.exists(raw_file):
        os.remove(raw_file)
        print(f"Removed old raw file: {raw_file}")

# Step 1: Generate the netlist from the .asc file (only once)
print("Generating netlist from .asc file...")
lt = SimCommander(ASC_FILE)
lt.run()  # This generates the .net file
if not lt.wait_completion(timeout=30):
    print("Failed to generate netlist!")
    sys.exit(1)

# Verify the netlist was generated
if not os.path.exists(NET_FILE):
    print(f"Error: {NET_FILE} not found! Netlist generation may have failed.")
    sys.exit(1)

# Function to modify the .net file by updating the Rout parameter
def update_net_file(r1_value):
    with open(NET_FILE, 'r') as file:
        net_content = file.read()
    
    # Debug: Print the netlist content before modification
    print(f"Netlist content before modification for R1 = {r1_value/1000:.1f} kΩ:")
    print(net_content)
    
    # Format R1 value in kΩ as an integer (e.g., 450k, 500k, etc.)
    r1_kohm = r1_value / 1000
    r1_str = f"{int(r1_kohm)}k"
    
    # Replace the .param Rout value in the netlist
    pattern = r'\.param Rout=\S+'
    replacement = f'.param Rout={r1_str}'
    new_net_content = re.sub(pattern, replacement, net_content)
    
    # Verify the replacement was successful
    if new_net_content == net_content:
        print(f"Error: Failed to update Rout parameter in .net file. Expected to replace '{pattern}' with '{replacement}'.")
        print("Current .net file content:")
        print(net_content)
        sys.exit(1)
    
    with open(NET_FILE, 'w') as file:
        file.write(new_net_content)
    
    # Debug: Print the netlist content after modification
    print(f"Netlist content after modification for R1 = {r1_value/1000:.1f} kΩ:")
    print(new_net_content)

# Run simulations for a range of R1 values
r1_values = np.linspace(R1_MIN, R1_MAX, R1_STEPS)  # e.g., [450k, 460k, ..., 550k]
freq = None
mag_data = []
phase_data = []

# Create a single SimCommander instance for the netlist file
lt = SimCommander(NET_FILE)  # Use the netlist file directly

for i, r1 in enumerate(r1_values, 1):  # Start counting from 1
    print(f"Simulating for R1 = {r1/1000:.1f} kΩ...")
    # Update the .net file by changing the Rout parameter
    update_net_file(r1)
    
    # Reset the netlist to ensure LTspice uses the updated file
    lt.reset_netlist()
    
    # Run the simulation using the modified netlist
    lt.run()  # This uses the modified .net file directly
    if not lt.wait_completion(timeout=30):  # Wait up to 30 seconds
        print(f"Simulation failed for R1 = {r1/1000:.1f} kΩ!")
        sys.exit(1)
    
    # Construct the correct .raw file name
    raw_file = f"{RAW_FILE_BASE}_{i}.raw"
    
    # Load the raw file
    if not os.path.exists(raw_file):
        print(f"Error: {raw_file} not found for R1 = {r1/1000:.1f} kΩ! Simulation may have failed.")
        sys.exit(1)
    
    try:
        l = ltspice.Ltspice(raw_file)
        l.parse()
    except Exception as e:
        print(f"Failed to parse raw file {raw_file} for R1 = {r1/1000:.1f} kΩ: {e}")
        sys.exit(1)
    
    # Extract frequency, magnitude, and phase data
    if freq is None:
        freq = l.get_frequency()  # Same for all simulations
    vout = l.get_data('V(out)')  # Matches FLAG 576 96 out
    if vout is None:
        print(f"Error: Could not extract V(out) data for R1 = {r1/1000:.1f} kΩ!")
        sys.exit(1)
    mag = 20 * np.log10(np.abs(vout))  # Magnitude in dB
    phase = np.angle(vout, deg=True)   # Phase in degrees
    
    mag_data.append(mag)
    phase_data.append(phase)

# Convert to numpy arrays for easier manipulation
mag_data = np.array(mag_data)  # Shape: (R1_STEPS, freq_points)
phase_data = np.array(phase_data)

# Debug: Print the magnitude and phase data to check for variation
print("Magnitude data (dB) for each R1 value:")
for i, r1 in enumerate(r1_values):
    print(f"R1 = {r1/1000:.1f} kΩ: {mag_data[i][:5]} ...")  # Print first 5 values for brevity
print("Phase data (deg) for each R1 value:")
for i, r1 in enumerate(r1_values):
    print(f"R1 = {r1/1000:.1f} kΩ: {phase_data[i][:5]} ...")  # Print first 5 values for brevity

# Calculate nominal, min, and max for magnitude and phase
nominal_idx = (R1_STEPS - 1) // 2  # Middle index for nominal (500 kΩ)
mag_nominal = mag_data[nominal_idx]
phase_nominal = phase_data[nominal_idx]
mag_min = np.min(mag_data, axis=0)  # Min across all R1 values
mag_max = np.max(mag_data, axis=0)  # Max across all R1 values
phase_min = np.min(phase_data, axis=0)
phase_max = np.max(phase_data, axis=0)

# Debug: Check the difference between min and max to ensure tolerance bands are visible
print(f"Magnitude tolerance range (dB): min = {mag_min[0]:.2f}, max = {mag_max[0]:.2f}, diff = {(mag_max[0] - mag_min[0]):.2f}")
print(f"Phase tolerance range (deg): min = {phase_min[0]:.2f}, max = {phase_max[0]:.2f}, diff = {(phase_max[0] - phase_min[0]):.2f}")

# Calculate pole and zero frequencies for nominal R1 (for markers)
R1 = R1_NOMINAL
R2 = 10e3     # 10 kΩ
C1 = 1000e-12 # 1000 pF
C2 = 100e-12  # 100 pF
fp1 = 1 / (2 * np.pi * (R1 + R2) * C1)  # Pole 1 frequency ~312 Hz
fz = 1 / (2 * np.pi * R2 * C1)          # Zero frequency ~15.9 kHz
fp2 = 1 / (2 * np.pi * R2 * C2)         # Pole 2 frequency ~159 kHz

# Calculate frequency of maximum phase (between fz and fp2) for nominal R1
f_max_phase = np.sqrt(fz * fp2)  # Geometric mean ~50.3 kHz

# Find indices for pole, zero, and max phase frequencies
fp1_idx = np.argmin(np.abs(freq - fp1))
fz_idx = np.argmin(np.abs(freq - fz))
fp2_idx = np.argmin(np.abs(freq - fp2))
f_max_idx = np.argmin(np.abs(freq - f_max_phase))

# Calculate phase boost as the difference from -90° for nominal R1
phase_at_max = phase_nominal[f_max_idx]  # Phase at 50.3 kHz
phase_boost_deg = phase_at_max - (-90)  # e.g., -33 - (-90) = 57°

# Create figure with three subplots: magnitude, phase, and schematic
fig = plt.figure(figsize=(12, 12))
gs = fig.add_gridspec(3, 1, height_ratios=[2, 2, 1])  # 2:2:1 ratio for heights

# Magnitude plot (first subplot)
ax1 = fig.add_subplot(gs[0, 0])
# Plot nominal curve
mag_plot, = ax1.semilogx(freq, mag_nominal, label='Magnitude (dB)', color='blue', linewidth=2)
# Add shaded region for tolerance range
ax1.fill_between(freq, mag_min, mag_max, color='blue', alpha=0.2, label='R1 Tolerance Range')
ax1.grid(True, which="both", ls="--")
ax1.set_ylabel('Magnitude (dB)')
ax1.set_title('Type II Compensator Frequency Response: Gain and Phase\nParameter Varied: R1, Tolerance: ±10%')

# Markers for poles and zero on magnitude plot
ax1.axvline(x=fp1, color='red', linestyle='--', label=f'Pole 1 at {fp1:.0f} Hz')
ax1.plot(fp1, mag_nominal[fp1_idx], 'rx', markersize=10)
ax1.axvline(x=fz, color='green', linestyle='--', label=f'Zero at {fz/1000:.1f} kHz')
ax1.plot(fz, mag_nominal[fz_idx], 'go', markersize=10)
ax1.axvline(x=fp2, color='purple', linestyle='--', label=f'Pole 2 at {fp2/1000:.0f} kHz')
ax1.plot(fp2, mag_nominal[fp2_idx], 'mo', markersize=10)
ax1.legend(loc='upper right')

# Phase plot (second subplot)
ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
# Plot nominal curve
phase_plot, = ax2.semilogx(freq, phase_nominal, label='Phase (deg)', color='purple', linewidth=2)
# Add shaded region for tolerance range
ax2.fill_between(freq, phase_min, phase_max, color='purple', alpha=0.2, label='R1 Tolerance Range')
ax2.grid(True, which="both", ls="--")
ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Phase (deg)')
ax2.set_ylim(-180, 90)  # Set max y-axis to 90°

# Markers for poles, zero, and phase boost on phase plot
ax2.axvline(x=fp1, color='red', linestyle='--')
ax2.plot(fp1, phase_nominal[fp1_idx], 'rx', markersize=10)
ax2.axvline(x=fz, color='green', linestyle='--')
ax2.plot(fz, phase_nominal[fz_idx], 'go', markersize=10)
ax2.axvline(x=fp2, color='purple', linestyle='--')
ax2.plot(fp2, phase_nominal[fp2_idx], 'mo', markersize=10)
ax2.axvline(x=f_max_phase, color='orange', linestyle='--', 
            label=f'Phase Boost at {f_max_phase/1000:.1f} kHz ({phase_boost_deg:.0f}°)')
ax2.plot(f_max_phase, phase_nominal[f_max_idx], 'yo', markersize=10)
ax2.legend(loc='upper right')

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
mag_text = ax1.text(0.05, 0.05, '', transform=ax1.transAxes, ha='left', va='bottom', fontsize=10)

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
    mag_dot.set_data([freq[idx]], [mag_nominal[idx]])
    # Calculate nominal, max, and min at this frequency
    nominal = mag_nominal[idx]
    max_diff = mag_max[idx] - nominal
    min_diff = nominal - mag_min[idx]
    mag_text.set_text(f'Magnitude\nFreq: {freq[idx]:.1e} Hz\nValue: {nominal:.2f} dB +{max_diff:.2f} dB -{min_diff:.2f} dB')
    fig.canvas.draw_idle()

# Custom interactive cursor for phase plot
phase_dot, = ax2.plot([], [], 'ko', markersize=5)  # Black dot for phase plot
phase_text = ax2.text(0.05, 0.05, '', transform=ax2.transAxes, ha='left', va='bottom', fontsize=10)

def on_move_phase(event):
    if event.inaxes != ax2:
        phase_dot.set_data([], [])
        phase_text.set_text('')
        fig.canvas.draw_idle()
        return
    x = event.xdata
    if x is None:
        return
    idx = np.argmin(np.abs(freq - x))
    phase_dot.set_data([freq[idx]], [phase_nominal[idx]])
    # Calculate nominal, max, and min at this frequency
    nominal = phase_nominal[idx]
    max_diff = phase_max[idx] - nominal
    min_diff = nominal - phase_min[idx]
    phase_text.set_text(f'Phase\nFreq: {freq[idx]:.1e} Hz\nValue: {nominal:.2f}° +{max_diff:.2f}° -{min_diff:.2f}°')
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