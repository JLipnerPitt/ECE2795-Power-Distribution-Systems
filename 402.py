import numpy as np

# Constants
j = 1j

# Phase conductor data
phase_GMR = 0.00446
phase_resistance = 1.12
phase_ncond = 3

# Neutral conductor data
neutral_GMR = 0.00446
neutral_resistance = 1.12
neutral_ncond = 1

ncond = phase_ncond + neutral_ncond

# Distance vector d
d = np.array([
    0 + j * 29,
    7 + j * 29,
    2.5+ j * 29,   
    4 + j * 25
])

# Resistance vector r
r = np.array([phase_resistance]*phase_ncond + [neutral_resistance])

# Distance matrix D
D = np.zeros((ncond, ncond), dtype=complex)
for i in range(ncond):
    for k in range(ncond):
        if i == k:
            D[i, k] = phase_GMR if i < phase_ncond else neutral_GMR
        else:
            D[i, k] = abs(d[i] - d[k])

# Primitive impedance matrix zprim
zprim = np.zeros((ncond, ncond), dtype=complex)
for i in range(ncond):
    for k in range(ncond):
        log_term = np.log(1 / D[i, k]) + 7.93402
        if i == k:
            zprim[i, k] = r[i] + 0.0953 + j * 0.12134 * log_term
        else:
            zprim[i, k] = 0.0953 + j * 0.12134 * log_term

# Partitioning zprim
zij = zprim[:phase_ncond, :phase_ncond]
zin = zprim[:phase_ncond, phase_ncond:]
znj = zprim[phase_ncond:, :phase_ncond]
znn = zprim[phase_ncond:, phase_ncond:][0, 0]  # Scalar

# Kron reduction
zabc = zij - (zin @ znj) / znn

# Neutral transformation matrix
tn = -znj / znn


# Output
np.set_printoptions(precision=4, suppress=True)

print("The distance matrix is\n")
print("D = \n", D, "\n")

print("The primitive impedance matrix in ohms/mile is\n")
print("[z] = \n", zprim, "\n")

print("In partitioned form, the submatrices in ohms/mile are\n")
print("[zij] = \n", zij, "\n")
print("[zin] = \n", zin, "\n")
print("[znj] = \n", znj, "\n")
print("[znn] = \n", znn, "\n")

print('The "Kron" reduced phase impedance matrix in ohms/mile is\n')
print("[zabc] = \n", zabc, "\n")

print("The neutral transformation matrix is\n")
print("[tn] = \n", tn)
