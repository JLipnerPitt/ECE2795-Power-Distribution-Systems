import numpy as np

# Constants
j = 1j

# Conductor data for line 1 and 2
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
cn2 = cn1.copy()

# Neutral conductor data
neutral = {'GMR': 0.01579, 'r': 0.303, 'ncond': 1}

ncond = cn1['ncond'] + cn2['ncond'] + neutral['ncond']

# Geometry and derived data
R1 = (cn1['dod'] - cn1['ds']) / 24
R2 = R1

GMRcn1 = (cn1['GMRs'] * cn1['k'] * R1**(cn1['k'] - 1))**(1 / cn1['k'])
GMRcn2 = GMRcn1

rcn1 = cn1['rs'] / cn1['k']
rcn2 = rcn1

print("R1 =\n", R1)
print("\nGMRcn1 =\n", GMRcn1)
print("\nrcn =\n", rcn1)

# Distance vector
d = np.array([
    0 + j*0, 4/12 + j*0, 8/12 + j*0,
    4/12 + j*(-10/12), 0 + j*(-10/12), 8/12 + j*(-10/12),
    0 + j*R1, 4/12 + j*R1, 8/12 + j*R1,
    4/12 + j*(R2 - 10/12), 0 + j*(R2 - 10/12), 8/12 + j*(R2 - 10/12),
    10/12 + j*(-5/12)
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
    else:
        r[i] = neutral['r']

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
            else:
                D[i, k] = neutral['GMR']
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
print("\n[z22abc] =\n", z22abc)
