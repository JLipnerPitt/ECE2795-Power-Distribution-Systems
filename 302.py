# Distribution System Modelling and Analysis, Example 3.4 Revised Routine
# Written by William Kersting and Robert Kerestes

import numpy as np

# Define complex line impedances (ohms)
Zline = np.array([
    (0.350000 + 0.617613j)*0.5,
    (0.350000 + 0.617613j)*0.65,
    (0.350000 + 0.617613j)*0.9,
])

# Define complex load impedances (ohms)
ZL = np.array([
    99999,  # very large (practically open)
    93.3 + 45.2j,
    36.7 + 22.8j,
    54.7 + 18.0j,
])

n = len(Zline)  # number of line impedances

# Source voltage
Vs = 12.47e3  # V

# Initial voltages at each bus (V)
V = np.array([
    12.47e3,
    12.47e3,
    12.47e3,
    12.47e3,
], dtype=complex)

m = len(V)  # number of buses

# Initialize current arrays
IL = np.zeros(m, dtype=complex)
I = np.zeros(m, dtype=complex)

# Backward sweep
for i in range(n, 0, -1):  # from n down to 0 (Python index n to 0)
    IL[i] = V[i] / ZL[i]
    I[i - 1] = I[i] + IL[i]
    V[i - 1] = V[i] + Zline[i - 1] * I[i - 1]

# Normalize voltages and currents to match source voltage
Ratio = Vs / V[0]
V *= Ratio
I *= Ratio
IL *= Ratio

# Display results
print("Nodal Voltages:")
for i in range(m):
    mag = np.abs(V[i])
    ang = np.angle(V[i], deg=True)
    print(f"V{i+1}: {mag:.3f} ∠ {ang:.2f}° V")

print("\nLine Currents:")
for i in range(m):
    mag = np.abs(I[i])
    ang = np.angle(I[i], deg=True)
    print(f"I{i+1}: {mag:.3f} ∠ {ang:.2f}° A")

print("\nLoad Currents:")
for i in range(m):
    mag = np.abs(IL[i])
    ang = np.angle(IL[i], deg=True)
    print(f"IL{i+1}: {mag:.3f} ∠ {ang:.2f}° A")
