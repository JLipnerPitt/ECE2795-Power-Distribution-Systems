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
    'dod': 1.10,
    'ds': 0.0641,
    'k': 7,
    'ncond': 6
}

# Additional neutral data
neutral1 = {'GMR': 0.00814, 'r': 0.592, 'ncond': 1}  # neutral for line1
neutral2 = {'GMR': 0.0051, 'r': 0.895, 'ncond': 1}  # neutral for line2
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
    0 + j*0, 4/12 + j*0, 8/12 + j*0,
    4/12 + j*(-10/12), 0 + j*(-10/12), 8/12 + j*(-10/12),
    0 + j*R1, 4/12 + j*R1, 8/12 + j*R1,
    4/12 + j*(R2 - 10/12), 0 + j*(R2 - 10/12), 8/12 + j*(R2 - 10/12),
    10/12 + j*(-5/12)
])
print(d)

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


