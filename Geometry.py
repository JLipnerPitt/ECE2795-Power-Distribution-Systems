"""
Module to implement geometry functionality

Filename: Geometry.py
Author: Justin Lipner
Date: 2025-01-23
"""

import numpy as np
from Conductor import Conductor

class Geometry:
    """
    Subclass geometry for distribution lines
    """
    def __init__(self, name: str, d: list[complex], nphases: int, phase_conductor: Conductor, neutral_conductor: Conductor):
        # pos = [[], []]
        self.name = name
        self.nphases = nphases
        self.phase_conductor = phase_conductor
        self.neutral_conductor = neutral_conductor
        self.ncond = len(d)
        self.nphases = nphases
        self.d = d
        self.D = self.calc_D()  # in meters
        self.Dshunt = self.calc_Dshunt()


    def calc_D(self):
        D = np.zeros((self.ncond, self.ncond), dtype=complex)
        for i in range(self.ncond):
            for k in range(self.ncond):
                if i == k:
                    D[i, k] = self.phase_conductor.GMR if i < self.nphases else self.neutral_conductor.GMR
                else:
                    D[i, k] = abs(self.d[i] - self.d[k])
        return D
    

    def calc_Dshunt(self):
        Dshunt = np.zeros((self.ncond, self.ncond), dtype=float)
        for i in range(self.ncond):
            for k in range(self.ncond):
                if i == k:
                    if i < self.ncond - 1:
                        Dshunt[i, k] = self.phase_conductor.radius
                    else:
                        Dshunt[i, k] = self.neutral_conductor.radius
                else:
                    Dshunt[i, k] = abs(self.d[i] - self.d[k])
        return Dshunt


# validation tests
if __name__ == '__main__':
    j = 1j
    from Geometry import Geometry
    from Conductor import Conductor
    phase_conductor = Conductor("1/0_ACSR", 0.398, 0.0446, 1.12, 230)
    neutral_conductor = Conductor("1/0_ACSR", 0.398, 0.0446, 1.12, 230)
    geometry1 = Geometry("Geometry 1", [0+j*29, 7+j*29, 2.5+j*29, 4+j*25], 3, phase_conductor, neutral_conductor)
    print("D =", geometry1.D, "m")
