# Distribution System Modelling and Analysis, Example 3.5
# Written by William Kersting and Robert Kerestes

import numpy as np

# Line impedances (ohms)
Sc = -200e3j
capacitor_node = 3

Zline = np.array([
    0.1739 + 0.3564j,
    0.1449 + 0.3970j,
    0.1159 + 0.2376j,
    0.2028 + 0.4158j
])

# Load powers in kVA (complex)
SL = np.array([
    0 + 0j,
    450 + 218j,
    1020 + 632j,
    855 + 281j,
])

count = len(SL)  # number of loads

# Initialization
Tol = 0.0001
Vs = 7200  # volts
maxit = 20

V_old = np.zeros(count, dtype=complex)
IL = np.zeros(count, dtype=complex)
I = np.zeros(count, dtype=complex)
V = np.zeros(count, dtype=complex)
Error = np.ones(count) * 1000
Error[0] = 0

# Set source voltage
V[0] = Vs

# Power flow iterations
for n in range(1, maxit + 1):
    # Forward sweep: update voltages
    for k in range(1, count):
        V[k] = V[k - 1] - Zline[k - 1] * I[k - 1]
        Error[k] = abs(abs(V[k]) - abs(V_old[k])) / Vs
        V_old[k] = V[k]

    Emax = np.max(Error)
    if Emax < Tol:
        break

    # Backward sweep: update currents
    for i in range(count-1, 0, -1):
        #  node capacitor bank is installed at
        if i == capacitor_node:
            IL[i] = np.conj(SL[i] * 1000 / V[i])
            Ic = Sc/V[i]
            I[i - 1] = IL[i] + I[i] - Ic
            continue

        IL[i] = np.conj(SL[i] * 1000 / V[i])
        I[i - 1] = IL[i] + I[i]

# Print results
print(f"Total iterations: {n}")

print("\nNodal Voltages:")
for i in range(count):
    mag = np.abs(V[i])
    ang = np.angle(V[i], deg=True)
    print(f"V{i+1}: {mag:.3f} ∠ {ang:.2f}° V")

print("\nLine Currents:")
for i in range(count):
    mag = np.abs(I[i])
    ang = np.angle(I[i], deg=True)
    print(f"I{i+1}: {mag:.3f} ∠ {ang:.2f}° A")

print("\nLoad Currents:")
for i in range(count):
    mag = np.abs(IL[i])
    ang = np.angle(IL[i], deg=True)
    print(f"IL{i+1}: {mag:.3f} ∠ {ang:.2f}° A")

# Percent voltage drop
Vrop = (np.abs(V[0]) - np.abs(V[count-1])) / np.abs(V[0]) * 100
print(f"\nPercent voltage drop from node 1 to node {count}: {Vrop:.2f}%")
