import numpy as np
import cmath

# Define zabc matrix manually (ohms/mile) taken from 402.py
zabc = np.array([
    [1.3238 + 1.3569j, 0.2066 + 0.4591j, 0.2101 + 0.5779j],
    [0.2066 + 0.4591j, 1.3294 + 1.3471j, 0.2130 + 0.5015j],
    [0.2101 + 0.5779j, 0.2130 + 0.5015j, 1.3368 + 1.3343j]], dtype=complex)

# Constants
j = 1j
length = 6000  # feet
Zabc = zabc * length / 5280  # Convert to ohms

# Load Impedance
Zload = 25+j*15
ZL = Zload*np.eye(3, dtype=complex)

Es = np.array([
    12470 / np.sqrt(3),
    12470 / np.sqrt(3) * cmath.exp(-j * 2 * np.pi / 3),
    12470 / np.sqrt(3) * cmath.exp(j * 2 * np.pi / 3)
]).reshape((3, 1))

Iabc = np.linalg.solve((Zabc+ZL), Es)
Vdrop = np.matmul(Zabc, Iabc)
Vabc = Es - Vdrop

# Convert to magnitude and phase
VLabc_mag = np.abs(Vabc)
VLabc_phase = np.angle(Vabc, deg=True)
Iabc_mag = np.abs(Iabc)
Iabc_phase = np.angle(Iabc, deg=True)


# Finding complex power for each load
Sabc = np.zeros((3, 1), dtype=complex)
for i in range(3):
    Sabc[i, 0] = VLabc_mag[i, 0]*cmath.exp(j*np.deg2rad(VLabc_phase[i, 0]))*np.conjugate(Iabc_mag[i, 0]*cmath.exp(j*np.deg2rad(Iabc_phase[i, 0])))

# Print results
phases = ['A', 'B', 'C']
for i, ph in enumerate(phases):
    print(f"Phase {ph}:")
    print(f"  VLN = {VLabc_mag[i, 0]:.2f}∠{VLabc_phase[i, 0]:.2f}° V")
    print(f"  ILoad = {Iabc_mag[i, 0]:.2f}∠{Iabc_phase[i, 0]:.2f}° A")
    print(f"  S = {Sabc[i, 0]/10e3} kVA")
    print(f"  Voltage Drop (%): {(np.abs(Vdrop[i, 0]/Es[i, 0])) * 100:.2f}\n")
