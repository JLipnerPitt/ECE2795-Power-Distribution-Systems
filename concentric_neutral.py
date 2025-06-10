import numpy as np

# Constants
j = 1j

# Concentric neutral data
cn = {
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
neutral = {'ncond': 0}  # No extra neutral

# Radius between concentric strands
R = (cn['dod'] - cn['ds']) / 24
print("R =\n", R)

# GMR of concentric neutral
GMRcn = (cn['GMRs'] * cn['k'] * R**(cn['k'] - 1))**(1 / cn['k'])
print("\nGMRcn =\n", GMRcn)

# Resistance per strand
rcn = cn['rs'] / cn['k']
print("\nrcn =\n", rcn)

ncond = cn['ncond'] + neutral['ncond']

# Resistance vector
r = np.zeros(ncond)
for i in range(ncond):
    if i < cn['ncond'] // 2:
        r[i] = cn['rc']
    else:
        r[i] = rcn

# Distance vector
d = np.array([
    0 + j * 0,
    0.5 + j * 0,
    1 + j * 0,
    0 + j * R,
    0.5 + j * R,
    1 + j * R
])

# Distance matrix D
D = np.zeros((ncond, ncond), dtype=complex)
for i in range(ncond):
    for k in range(ncond):
        if i == k:
            D[i, k] = cn['GMRc'] if i < cn['ncond'] // 2 else GMRcn
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
p = cn['ncond'] // 2
zij = zprim[:p, :p]
zin = zprim[:p, p:]
znj = zprim[p:, :p]
znn = zprim[p:, p:]

# Kron reduction
zabc = zij - zin @ np.linalg.inv(znn) @ znj

# Neutral transformation matrix
tn = -np.linalg.inv(znn) @ znj

# Sequence transformation
a = np.exp(j * 2 * np.pi / 3)
As = np.array([
    [1, 1, 1],
    [1, a**2, a],
    [1, a, a**2]
], dtype=complex)
As_inv = np.linalg.inv(As)
z012 = As_inv @ zabc @ As

# Output
np.set_printoptions(precision=4, suppress=True)

print("\nThe distance matrix D is:\n")
print(D)

print("\nThe primitive impedance matrix zprim in ohms/mile is:\n")
print(zprim)

print("\nSubmatrices in partitioned form:")
print("[zij] =\n", zij)
print("[zin] =\n", zin)
print("[znj] = [zin] =\n", znj)
print("[znn] =\n", znn)

print("\nThe 'Kron' reduced phase impedance matrix zabc in ohms/mile is:\n")
print(zabc)

print("\nThe sequence impedance matrix z012 in ohms/mile is:\n")
print(z012)
