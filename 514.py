import numpy as np

# Complex unit
j = 1j

# Concentric neutral data for line 1
cn1 = {
    'GMRc': 0.0171,
    'dc': 0.567,
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
    'dc': 0.567,
    'rc': 0.41,
    'GMRs': 0.00208,
    'rs': 14.8722,
    'dod': 1.29,
    'ds': 0.0641,
    'k': 13,
    'ncond': 6
}

# Additional neutral (unused here)
neutral = {'GMR': 0.01579, 'dn': 0.522, 'r': 0.303, 'ncond': 1}

# Radii of conductor and strand (inches)
RDc1 = cn1['dc'] / 2.0
RDs1 = cn1['ds'] / 2.0
RDc2 = cn2['dc'] / 2.0
RDs2 = cn2['ds'] / 2.0

# Radii to concentric neutrals (feet → inches)
R1 = (cn1['dod'] - cn1['ds']) / 24.0 * 12.0
R2 = (cn2['dod'] - cn2['ds']) / 24.0 * 12.0

# Number of phase conductors per line
n1 = cn1['ncond'] // 2
n2 = cn2['ncond'] // 2
n_total = n1 + n2

# Pre-allocate the phase admittance matrix
yabc = np.zeros((n_total, n_total), dtype=complex)

# Helper for the admittance formula
def yc(R, RDc, RDs, k):
    num = j * 77.3619
    den = np.log(R / RDc) - (1.0 / k) * np.log(k * RDs / R)
    return num / den

# Fill block for line 1 (rows 0..n1-1, cols 0..n1-1)
for i in range(n1):
    for k in range(n1):
        if i == k:
            yabc[i, k] = yc(R1, RDc1, RDs1, cn1['k'])
        # off-diagonals remain zero

# Fill block for line 2 (rows n1..n1+n2-1, cols n1..n1+n2-1)
for i in range(n1, n1 + n2):
    for k in range(n1, n1 + n2):
        if i == k:
            yabc[i, k] = yc(R2, RDc2, RDs2, cn2['k'])
        # off-diagonals remain zero

# Display results
print(f"R1 (inches) = {R1}")
print(f"RDc1 (inches) = {RDc1}")
print(f"RDs1 (inches) = {RDs1}")
print(f"yc (line 1 self-admittance) = {yabc[0,0]}\n")
print("The phase admittance matrix [yabc] is:")
print(yabc)
