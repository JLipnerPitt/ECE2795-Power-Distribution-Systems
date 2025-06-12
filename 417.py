import numpy as np

# Constants
j = 1j

# Concentric neutral data for line 1
cn1 = {
    'GMRc': 0.0171,
    'rc': 0.41,
    'GMRs': 0.00208,
    'rs': 14.8722,
    'dod': 1.29,
    'ds': 0.0641,
    'k': 13,
    'ncond': 6
}

# Concentric neutral data for line 2
cn2 = {
    'GMRc': 0.0171,
    'rc': 0.41,
    'GMRs': 0.00208,
    'rs': 14.8722,
    'dod': 1.29,
    'ds': 0.0641,
    'k': 13,
    'ncond': 6
}

# Additional neutral data
neutral1 = {'GMR': 0.00814, 'r': 0.592, 'ncond': 1}  # neutral for line1 4/0 6/1 ACSR
neutral2 = {'GMR': 0.0051, 'r': 0.895, 'ncond': 1}  # neutral for line2 2/0 ACSR
ncond = cn1['ncond'] + cn2['ncond'] + neutral1['ncond'] + neutral2['ncond']

# Geometry and derived data
R1 = (cn1['dod'] - cn1['ds']) / 24
R2 = (cn2['dod'] - cn2['ds']) / 24

GMRcn1 = (cn1['GMRs'] * cn1['k'] * R1**(cn1['k'] - 1))**(1 / cn1['k'])
GMRcn2 = (cn2['GMRs'] * cn2['k'] * R2**(cn2['k'] - 1))**(1 / cn2['k'])

rcn1 = cn1['rs'] / cn1['k']
rcn2 = cn2['rs'] / cn2['k']

print(f"R1 = {R1}, R2 = {R2}")
print(f"GMRcn1 = {GMRcn1}, GMRcn2 = {GMRcn2}")
print(f"rcn1 = {rcn1}, rcn2 = {rcn2}")

# Distance vector
d = np.array([
    0 + j*0, 6/12 + j*0, 12/12 + j*0, 
    6/12 + j*(-24/12), 12/12 + j*(-10/12), -24/12 + 0j,
    0 + j*R1, 6/12 + j*R1, 12/12 + j*R1,
    6/12 + j*(R2 - 24/12), 12/12 + j*(R2 - 24/12), 0 + j*(R2 - 24/12),
    16/12 +j*0, 16/12 + j*(-24/12)
])

# Resistance vector r
r = np.zeros(ncond)
for i in range(ncond):
    if i < cn1['ncond'] // 2:
        r[i] = cn1['rc']
    elif i < cn1['ncond'] // 2 + cn2['ncond'] // 2:
        r[i] = cn2['rc']
    elif i < cn1['ncond'] + cn2['ncond'] // 2:
        r[i] = rcn1
    elif i < cn1['ncond'] + cn2['ncond']:
        r[i] = rcn2
    elif i == ncond-1:
        r[i] = neutral1['r']
    else:
        r[i] = neutral2['r']

# Distance matrix D
D = np.zeros((ncond, ncond), dtype=complex)
for i in range(ncond):
    for k in range(ncond):
        if i == k:
            if i < cn1['ncond'] // 2:
                D[i, k] = cn1['GMRc']
            elif i < cn1['ncond'] // 2 + cn2['ncond'] // 2:
                D[i, k] = cn2['GMRc']
            elif i < cn1['ncond'] + cn2['ncond'] // 2:
                D[i, k] = GMRcn1
            elif i < cn1['ncond'] + cn2['ncond']:
                D[i, k] = GMRcn2
            elif i == ncond-1:
                D[i, k] = neutral1['GMR']
            else:
                D[i, k] = neutral2['GMR']
        else:
            D[i, k] = abs(d[i] - d[k])

# Primitive impedance matrix
zprim = np.zeros((ncond, ncond), dtype=complex)
for i in range(ncond):
    for k in range(ncond):
        log_term = np.log(1 / D[i, k]) + 7.93402
        if i == k:
            zprim[i, k] = r[i] + 0.0953 + j * 0.12134 * log_term
        else:
            zprim[i, k] = 0.0953 + j * 0.12134 * log_term

# Partitioning
p_total = cn1['ncond'] // 2 + cn2['ncond'] // 2
zij = zprim[:p_total, :p_total]
zin = zprim[:p_total, p_total:]
znj = zprim[p_total:, :p_total]
znn = zprim[p_total:, p_total:]

# Kron reduction
zabc = zij - zin @ np.linalg.inv(znn) @ znj

# Neutral transformation matrix
tn = -np.linalg.inv(znn) @ znj

# Partitioned phase matrices
z11abc = zabc[0:3, 0:3]
z12abc = zabc[0:3, 3:6]
z21abc = zabc[3:6, 0:3]
z22abc = zabc[3:6, 3:6]

# Output
np.set_printoptions(precision=4, suppress=True)

print("\nDistance matrix D =\n", D)
print("\nPartitioned phase impedance matrices:\n")
print("[z11abc] =\n", z11abc)
print("\n[z12abc] =\n", z12abc)
print("\n[z21abc] =\n", z21abc)
print("\n[z22abc] =\n", z22abc, "\n")

print('The "Kron" reduced phase impedance matrix in ohms\mile is\n')
print("[zabc] = \n", zabc, "\n")