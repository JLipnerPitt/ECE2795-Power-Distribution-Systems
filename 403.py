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

Tol = 0.0001
start = np.zeros((3, 1), dtype=complex)

Iabc = start.copy()
Vold = start.copy()
error = start.copy()
Vdrop = np.zeros((3, 1))

# Load power (complex power per phase) and source voltages
SL = np.array([5000000 * cmath.exp(j * cmath.acos(0.85))] * 3).reshape((3, 1))
Es = np.array([
    12470 / np.sqrt(3),
    12470 / np.sqrt(3) * cmath.exp(-j * 2 * np.pi / 3),
    12470 / np.sqrt(3) * cmath.exp(j * 2 * np.pi / 3)
]).reshape((3, 1))

# Iteration
for n in range(1, 31):
    VLabc = Es - Zabc @ Iabc

    for i in range(3):
        error[i] = abs(abs(VLabc[i]) - abs(Vold[i]))

    Err = np.max(error)
    if Err < Tol:
        break

    Vold = VLabc.copy()

    for i in range(3):
        Iabc[i] = SL[i] / VLabc[i]

    for i in range(3):
        Vdrop[i] = abs(abs(Es[i]) - abs(VLabc[i])) / abs(Es[i])

# Convert to magnitude and phase
VLabc_mag = np.abs(VLabc)
VLabc_phase = np.angle(VLabc, deg=True)
Iabc_mag = np.abs(Iabc)
Iabc_phase = np.angle(Iabc, deg=True)


# Finding complex power for each load
Sabc = np.zeros((3, 1), dtype=complex)
for i in range(3):
    Sabc[i, 0] = VLabc_mag[i, 0]*cmath.exp(j*np.deg2rad(VLabc_phase[i, 0]))*np.conjugate(Iabc_mag[i, 0]*cmath.exp(j*np.deg2rad(Iabc_phase[i, 0])))

# Print results
phases = ['A', 'B', 'C']
print(f"Converged in {n} iterations\n")
for i, ph in enumerate(phases):
    print(f"Phase {ph}:")
    print(f"  VLN = {VLabc_mag[i, 0]:.2f}∠{VLabc_phase[i, 0]:.2f}° V")
    print(f"  ILoad = {Iabc_mag[i, 0]:.2f}∠{Iabc_phase[i, 0]:.2f}° A")
    print(f"  S = {Sabc[i, 0]/10e3} kVA")
    print(f"  Voltage Drop (%): {Vdrop[i, 0] * 100:.2f}\n")
