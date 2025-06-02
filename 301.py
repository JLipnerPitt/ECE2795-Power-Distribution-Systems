# Distribution System Modelling and Analysis, Example 3.1
# Written by William Kersting and Robert Kerestes

import numpy as np

# Defining the phase conductor data
phase_GMR = 0.0217  # feet
phase_resistance = 0.350  # ohms/mile

# Defining the spacings between conductors (in feet)
Dab = 2.5
Dbc = 7.0
Dca = 2.5

# Computing the equivalent spacing
Deq = (Dab * Dbc * Dca) ** (1 / 3)

# Using equation 3.4 to calculate line impedance
zline = phase_resistance + 1j * 0.12134 * np.log(Deq / phase_GMR)

# Display results
print(f"The equivalent spacing in feet is: {Deq:.3f}")
print(f"The line impedance for this line is: {zline.real:.3f} + j{zline.imag:.3f} ohms/mile")
