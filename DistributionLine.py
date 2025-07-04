import numpy as np
from Bus import Bus
from Geometry import Geometry
from Settings import settings
j = 1j

class DistributionLine:
    """
    DistributionLine class to hold distribution line information
    """

    def __init__(self, name, bus1: Bus, bus2: Bus, geometry: Geometry, length: float):
        self.name = name
        self.bus1 = bus1
        self.bus2 = bus2
        self.geometry = geometry
        self.length = length
        self.freq = settings.freq
        self.powerbase = settings.powerbase
        #self.Zbase = self.bus1.base_kv**2/self.powerbase
        self.Zprim = self.calc_Zprim()
        self.Zabc, self.tn = self.partition_Zprim()
        self.Pprim = self.calc_Pprim()
        self.Yabc = self.partition_Pprim()
    

    def calc_Zprim(self):
        """
        Calculate line series resistance from bundle information and length
        :return: Series resistance in pu (float)
        """
        # Resistance vector r
        r = np.array([self.geometry.phase_conductor.resistance]*self.geometry.nphases + [self.geometry.neutral_conductor.resistance])

        # Primitive impedance matrix
        zprim = np.zeros((self.geometry.ncond, self.geometry.ncond), dtype=complex)

        for i in range(self.geometry.ncond):
            for k in range(self.geometry.ncond):
                log_term = np.log(1 / self.geometry.D[i, k]) + 7.93402
                if i == k:
                    zprim[i, k] = r[i] + 0.0953 + j * 0.12134 * log_term
                else:
                    zprim[i, k] = 0.0953 + j * 0.12134 * log_term

        return zprim
    

    def partition_Zprim(self):
        # Partitioning zprim
        zij = self.Zprim[:self.geometry.nphases, :self.geometry.nphases]
        zin = self.Zprim[:self.geometry.nphases, self.geometry.nphases:]
        znj = self.Zprim[self.geometry.nphases:, :self.geometry.nphases]
        znn = self.Zprim[self.geometry.nphases:, self.geometry.nphases:][0, 0]  # Scalar

        # Kron reduction
        zabc = zij - (zin @ znj) / znn
        Zabc = self.length*zabc

        # Neutral transformation matrix
        tn = -znj / znn

        mask = np.array(self.geometry.phases, dtype=bool)
        # zero out unwanted rows
        Zabc[~mask, :] = 0
        # zero out unwanted columns
        Zabc[:, ~mask] = 0
        return Zabc, tn


    def calc_Pprim(self):
        # Pre‐allocate matrices
        S = np.zeros((self.geometry.ncond, self.geometry.ncond), dtype=complex)
        Pprim = np.zeros((self.geometry.ncond, self.geometry.ncond), dtype=complex)

        # Build image distances S
        for i in range(self.geometry.ncond):
            for k in range(self.geometry.ncond):
                S[i, k] = abs(self.geometry.d[i] - np.conj(self.geometry.d[k]))

       
        # Primitive potential coefficient matrix
        # 11.17689 = 1/(2*pi*ε0) in appropriate units
        for i in range(self.geometry.ncond):
            for k in range(self.geometry.ncond):
                Pprim[i, k] = 11.17689 * np.log(S[i, k] / self.geometry.Dshunt[i, k])

        return Pprim
    

    def partition_Pprim(self):
        # Partition for Kron reduction
        Pij = self.Pprim[:self.geometry.nphases, :self.geometry.nphases]
        Pin = self.Pprim[:self.geometry.nphases, self.geometry.nphases:]
        Pnj = self.Pprim[self.geometry.nphases:, :self.geometry.nphases]
        Pnn = self.Pprim[self.geometry.nphases, self.geometry.nphases]

        # Kron reduction: Pabc = Pij – Pin * (1/Pnn) * Pnj
        Pabc = Pij - (Pin @ Pnj) / Pnn
        # Phase capacitance matrix (µF/mile)
        Cabc = np.linalg.inv(Pabc)

        # Shunt admittance matrix (µS/mile)
        yabc = j * 2 * np.pi * self.freq * Cabc
        Yabc = self.length*yabc*10e-6

        mask = np.array(self.geometry.phases, dtype=bool)
        # zero out unwanted rows
        Yabc[~mask, :] = 0
        # zero out unwanted columns
        Yabc[:, ~mask] = 0
        return Yabc


# validation tests
if __name__ == '__main__':
    from Geometry import Geometry
    from Conductor import Conductor
    from DistributionLine import DistributionLine
    
    phase_conductor = Conductor("1/0_ACSR", 0.398, 0.00446, 1.12, 230)
    neutral_conductor = Conductor("1/0_ACSR", 0.398, 0.00446, 1.12, 230)
    geometry1 = Geometry("Geometry 1", [0+j*29, 7+j*29, 2.5+j*29, 4+j*25], phase_conductor, neutral_conductor, [0, 1, 0])
    line1 = DistributionLine("OH1", "bus1", "bus2", geometry1, 1.893939)

    print("The primitive impedance matrix in ohms/mile is\n")
    print("[z] = \n", line1.Zprim, "\n")

    print('The "Kron" reduced phase impedance matrix in ohms/mile is\n')
    print("[zabc] = \n", line1.Zabc, "\n")

    print("\nShunt admittance matrix yabc (S/mile):")
    print(line1.Yabc)
    print()

    
    phase_conductor = Conductor("336400_26/7_ACSR", 0.721, 0.0244, 0.306, 530)
    neutral_conductor = Conductor("4/0_6/1_ACSR", 0.563, 0.00814, 0.592, 340)
    geometry2 = Geometry("Geometry 2", [0 + j*29, 7+j*29, 2.5 + j*29, 4+j*25], phase_conductor, neutral_conductor)
    line2 = DistributionLine("OH1", "bus1", "bus2", geometry2, 1.893939)

    print("The primitive impedance matrix in ohms/mile is\n")
    print("[z] = \n", line2.Zprim, "\n")

    print('The "Kron" reduced phase impedance matrix in ohms/mile is\n')
    print("[zabc] = \n", line2.Zabc, "\n")

    print("\nShunt admittance matrix yabc (S/mile):")
    print(line2.Yabc)
    
