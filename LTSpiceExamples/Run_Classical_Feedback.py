# Run_Classical_Feedback.py
from PyLTSpice import SimCommander  # For running LTSpice
import ltspice                      # For reading .raw files
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os
import sys                          # Added for sys.exit()
import mplcursors                   # For interactive cursor
import time                         # For adding a delay

# Switch to an interactive backend
#matplotlib.use('Qt5Agg')  # Separate window.  Useful for Spyder.  Requires "pip install PyQt5"
matplotlib.use('TkAgg')

# Configuration
LTSPICE_PATH = "C:/Program Files/ADI/LTspice/LTspice.exe"  # Update this path if needed
ASC_FILE = "Classical_Feedback_Example_Bode.asc"
RAW_FILE_BASE = "Classical_Feedback_Example_Bode"

# Print current working directory for debugging
print(f"Current working directory: {os.getcwd()}")

# Check if the .asc file exists
if not os.path.exists(ASC_FILE):
    print(f"Error: {ASC_FILE} not found in the current directory!")
    sys.exit(1)

# Initialize SimCommander
lt = SimCommander(ASC_FILE)

# Run the simulation
print("Starting LTSpice simulation...")
lt.run()
lt.wait_completion()
print("Simulation complete.")

# Add a small delay to ensure the .raw file is written
time.sleep(1)  # Wait 1 second

# Look for the .raw file (PyLTSpice may append a suffix like _1, _2, etc.)
raw_file = None
for i in range(1, 10):  # Check for files like _1.raw, _2.raw, etc.
    potential_raw_file = f"{RAW_FILE_BASE}_{i}.raw"
    print(f"Checking for {potential_raw_file}...")
    if os.path.exists(potential_raw_file):
        raw_file = potential_raw_file
        break

# If no suffixed file is found, try the base name
if raw_file is None:
    raw_file = f"{RAW_FILE_BASE}.raw"
    print(f"Checking for {raw_file}...")
    if not os.path.exists(raw_file):
        print(f"Error: No suitable .raw file found! Simulation may have failed.")
        print("Listing files in directory:")
        print(os.listdir(os.getcwd()))
        sys.exit(1)

print(f"Using raw file: {raw_file}")

# Parse the raw file
l = ltspice.Ltspice(raw_file)
l.parse()

# Extract frequency and node voltages
freq = l.get_frequency()
v_output = l.get_data('V(output)')
v_ith = l.get_data('V(ith)')
v_error = l.get_data('V(error)')
v_fb = l.get_data('V(fb)')

# Check if node voltages were extracted successfully
if any(v is None for v in [freq, v_output, v_ith, v_error, v_fb]):
    print("Error: Failed to extract one or more node voltages from the raw file.")
    print(f"freq: {freq is not None}")
    print(f"v_output: {v_output is not None}")
    print(f"v_ith: {v_ith is not None}")
    print(f"v_error: {v_error is not None}")
    print(f"v_fb: {v_fb is not None}")
    sys.exit(1)

# Debug: Print node voltages at selected frequencies
print("Debug: Node voltages at selected frequencies:")
for i in range(0, len(freq), len(freq)//5):  # Sample 5 points
    print(f"Freq: {freq[i]:.1e} Hz, V(Error): {v_error[i]:.2e}, V(FB): {v_fb[i]:.2e}")

# Compute transfer functions
transfer_functions = {
    "Loop Gain": v_output / v_fb,  # Top-left (index 0)
    "Plant": v_output / v_ith,     # Top-right (index 1)
    "Compensator": v_ith / v_error, # Bottom-left (index 2)
    "Feedback": v_error / v_fb     # Bottom-right (index 3) - Plotting V(Error)/V(FB)
}

# Debug: Compute and print Feedback transfer function at selected frequencies
print("Debug: Feedback transfer function V(Error)/V(FB) at selected frequencies:")
for i in range(0, len(freq), len(freq)//5):  # Sample 5 points
    tf = transfer_functions["Feedback"][i]
    mag = 20 * np.log10(np.abs(tf))
    phase = np.angle(tf, deg=True)
    adjusted_phase = phase - 180
    # Wrap the adjusted phase to [-180, 180]
    adjusted_phase = (adjusted_phase + 180) % 360 - 180
    print(f"Freq: {freq[i]:.1e} Hz, Magnitude: {mag:.2f} dB, Phase: {phase:.2f}°, Adjusted Phase (phase - 180): {adjusted_phase:.2f}°")

# Create figure with 2x2 subplots
fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
axes = axes.flatten()  # Flatten to easily index as [0, 1, 2, 3]
titles = ["Loop Gain: V(Output)/V(FB)", "Plant: V(Output)/V(Ith)", 
          "Compensator: V(Ith)/V(Error)", "Feedback: V(Error)/V(FB)"]
colors = ['purple', 'blue', 'green', 'red']

# Lists to store secondary axes and crossover frequency
secondary_axes = []  # To store ax2 for each subplot
loop_gain_crossover_freq = None  # To store the 0 dB crossover frequency of Loop Gain

# Plot each transfer function
for idx, (key, tf) in enumerate(transfer_functions.items()):
    ax = axes[idx]
    
    # Calculate magnitude and phase
    mag = 20 * np.log10(np.abs(tf))  # Magnitude in dB
    phase = np.angle(tf, deg=True)   # Phase in degrees
    
    # Adjust phase for Feedback plot: phase - 180
    if key == "Feedback":
        phase = phase - 180
        # Wrap the phase to [-180, 180] for plotting
        phase = (phase + 180) % 360 - 180
    
    # Calculate Bandwidth and Gain Margin for Loop Gain only
    if key == "Loop Gain":
        # Bandwidth (frequency where magnitude crosses 0 dB)
        gain_crossover_idx = np.where(mag <= 0)[0]
        if len(gain_crossover_idx) > 0:
            loop_gain_crossover_freq = freq[gain_crossover_idx[0]]
            phase_at_crossover = phase[gain_crossover_idx[0]]
        else:
            loop_gain_crossover_freq = freq[-1]  # Fallback if no crossing
            phase_at_crossover = phase[-1]
        
        # Debug: Print the phase at crossover
        print(f"Phase at crossover frequency ({loop_gain_crossover_freq/1000:.2f} kHz): {phase_at_crossover:.1f}°")
        
        # Calculate Phase Margin (phase at gain crossover - 0)
        phase_margin = phase_at_crossover - 0
        print(f"Calculated Phase Margin: {phase_margin:.1f}°")
        
        # Calculate Gain Margin (gain at the frequency where phase crosses -180°)
        phase_crossover_idx = np.where(phase <= -180)[0]
        if len(phase_crossover_idx) > 0:
            gain_margin_freq = freq[phase_crossover_idx[0]]
            gain_margin = -mag[phase_crossover_idx[0]]  # Negative of the gain (since gain is in dB)
        else:
            gain_margin = float('inf')  # If phase never crosses -180°, gain margin is infinite
            gain_margin_freq = freq[-1]
    
    # For other plots, get gain and phase at the Loop Gain crossover frequency
    else:
        crossover_idx = np.abs(freq - loop_gain_crossover_freq).argmin()
        gain_at_crossover = mag[crossover_idx]
        phase_at_crossover = phase[crossover_idx]
    
    # Plot Magnitude (primary y-axis)
    mag_plot, = ax.semilogx(freq, mag, color=colors[idx], label='Gain')
    ax.grid(True, which="both", ls="--")
    ax.set_ylabel('Magnitude (dB)', color=colors[idx])
    ax.tick_params(axis='y', labelcolor=colors[idx])
    ax.set_title(titles[idx])
    
    # Set custom y-axis limits for Feedback plot (magnitude)
    if key == "Feedback":
        ax.set_ylim(-20, 10)  # Min: -20 dB, Max: 10 dB
    
    # Add vertical marker at the Loop Gain crossover frequency
    ax.axvline(x=loop_gain_crossover_freq, color='black', linestyle='--', alpha=0.5)
    
    # Create secondary y-axis for Phase
    ax2 = ax.twinx()
    phase_plot, = ax2.semilogx(freq, phase, color=colors[idx], linestyle='-.', label='Phase')
    ax2.set_ylabel('Phase (deg)', color=colors[idx])
    
    # Set custom y-axis limits for Feedback plot (phase)
    if key == "Feedback":
        ax2.set_ylim(-90, 90)  # Adjusted to show phase increasing from 0° to 90°
    else:
        ax2.set_ylim(-180, 180)  # Default phase limits for other plots: ±180°
    
    ax2.tick_params(axis='y', labelcolor=colors[idx])
    
    # Store the secondary axis
    secondary_axes.append(ax2)
    
    # Add text box
    props = dict(boxstyle='round', facecolor='white', alpha=0.8)
    if key == "Loop Gain":
        # For Loop Gain, show BW, PM, and GM
        if gain_margin == float('inf'):
            textstr = f'BW: {loop_gain_crossover_freq/1000:.2f} kHz\nPM: {phase_margin:.1f}°\nGM: ∞ dB'
        else:
            textstr = f'BW: {loop_gain_crossover_freq/1000:.2f} kHz\nPM: {phase_margin:.1f}°\nGM: {gain_margin:.1f} dB'
    else:
        # For other plots, show Gain and Phase at crossover frequency (use "BW" instead of "CF")
        textstr = f'Gain at BW: {gain_at_crossover:.2f} dB\nPhase at BW: {phase_at_crossover:.1f}°'
    ax.text(0.95, 0.95, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', horizontalalignment='right', bbox=props)

# Set x-label on the bottom subplots
axes[2].set_xlabel('Frequency (Hz)')
axes[3].set_xlabel('Frequency (Hz)')

# Add interactive cursor using stored secondary axes
cursor = mplcursors.cursor([ax.lines[0] for ax in axes] + [ax2.lines[0] for ax2 in secondary_axes], hover=True)
@cursor.connect("add")
def on_add(sel):
    x, y = sel.target
    sel.annotation.set_text(f'Freq: {x:.1e} Hz\nValue: {y:.2f}')

# Adjust layout and show
plt.tight_layout()
plt.show()

# Run the cleanup script
os.system("python LTSpice_Cleanup.py")