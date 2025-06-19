import numpy as np

# Define constants
j = 1j
f = 60.0  # Hz

# Phase conductor data
phase = {
    'GMR': 0.0244,
    'resistance': 0.306,
    'diameter': 0.721,
    'ncond': 3
}

# Neutral conductor data
neutral = {
    'GMR': 0.00814,
    'resistance': 0.592,
    'diameter': 0.563,
    'ncond': 1
}

ncond = phase['ncond'] + neutral['ncond']

# Pre‐allocate matrices
Dshunt = np.zeros((ncond, ncond), dtype=float)
S = np.zeros((ncond, ncond), dtype=float)
Pprim = np.zeros((ncond, ncond), dtype=float)

# Conductor coordinates (complex numbers in feet)
d = np.array([0 + j*29, 2.5 + j*29, 7 + j*29, 4 + j*25])

# Build Dshunt (self‐distances = diameter/24, mutual = geometric distance)
for i in range(ncond):
    for k in range(ncond):
        if i == k:
            if i < ncond - 1:
                Dshunt[i, k] = phase['diameter'] / 24.0
            else:
                Dshunt[i, k] = neutral['diameter'] / 24.0
        else:
            Dshunt[i, k] = abs(d[i] - d[k])

# Build image distances S
for i in range(ncond):
    for k in range(ncond):
        S[i, k] = abs(d[i] - np.conj(d[k]))

# Primitive potential coefficient matrix
# 11.17689 = 1/(2*pi*ε0) in appropriate units
for i in range(ncond):
    for k in range(ncond):
        Pprim[i, k] = 11.17689 * np.log(S[i, k] / Dshunt[i, k])

# Partition for Kron reduction
Pij = Pprim[:phase['ncond'], :phase['ncond']]
Pin = Pprim[:phase['ncond'], phase['ncond']:]
Pnj = Pprim[phase['ncond']:, :phase['ncond']]
Pnn = Pprim[phase['ncond'], phase['ncond']]

# Kron reduction: Pabc = Pij – Pin * (1/Pnn) * Pnj
Pabc = Pij - (Pin @ Pnj) / Pnn

# Phase capacitance matrix (µF/mile)
Cabc = np.linalg.inv(Pabc)

# Shunt admittance matrix (µS/mile)
yabc = j * 2 * np.pi * f * Cabc

# Display results
print("Image distance matrix S (feet):")
print(S)
print("\nPrimitive potential coefficient matrix Pprim (mile/µF):")
print(Pprim)
print("\nKron‐reduced phase potential coefficient matrix Pabc (mile/µF):")
print(Pabc)
print("\nPhase capacitance matrix Cabc (µF/mile):")
print(Cabc)
print("\nShunt admittance matrix yabc (µS/mile):")
print(yabc)