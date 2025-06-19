import numpy as np

# Define constants
j = 1j

# Concentric neutral data for the first line
cn = {
    'GMRc': 0.0171,
    'rc': 0.41,
    'diameter': 0.567,
    'GMRs': 0.00208,
    'rs': 14.8722,
    'dod': 1.29,
    'ds': 0.0641,
    'k': 13,
    'ncond': 6
}

# No additional neutral in this example
neutral = {'ncond': 0}

# Compute R = (dod − ds)/24 in feet, then convert to inches
R = (cn['dod'] - cn['ds']) / 24.0
R *= 12.0  # inches

# Radii of conductor and strand (in inches)
RDc = cn['diameter'] / 2.0
RDs = cn['ds'] / 2.0

print("R (inches) =")
print(R)

# Number of phase conductors (half of concentric neutral pairs)
n_phase = cn['ncond'] // 2

# Pre‐allocate phase admittance matrix
yabc = np.zeros((n_phase, n_phase), dtype=complex)

# Compute primitive admittance matrix
# 77.3619 = 2πfε0 in appropriate units
for i in range(n_phase):
    for k in range(n_phase):
        if i == k:
            numerator = j * 77.3619
            denom = (np.log(R / RDc)
                     - (1.0 / cn['k']) * np.log(cn['k'] * RDs / R))
            yabc[i, k] = numerator / denom
        else:
            yabc[i, k] = 0.0

print("\nThe phase admittance matrix [yabc] for the three-phase line is:")
print(yabc)