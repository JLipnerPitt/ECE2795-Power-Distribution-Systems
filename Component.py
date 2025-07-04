#  This class contains various components used in electrical circuits. 
#  Component is a parent class for all the child "component" classes.
from Settings import settings
from math import sin, cos
import numpy as np
j = 1j

class Load:
    """
    Unbalanced, multiphase Load.
    """
    def __init__(self,
                 name: str,
                 bus: str,
                 kva: list[float],
                 pf: list[float],
                 type: str = 'PQ',
                 connection: str = 'Y',
                 phases=None):
        """
        :param name:           Load name
        :param bus:            Bus this load is connected to
        :param real_power:     list of real powers for phase a, b, c
        :param pf:             list of each phases power factor
        :param connection:     'Y' or 'Δ' defaults to 'Y'
        :param phases:         list of phase labels (defaults to keys of real_power)
        """
        self.name       = name
        self.bus        = bus
        self.connection = connection
        self.phases     = ['A','B','C'] if phases is None else phases
        self.kva = np.array(kva)*1e3
        self.pf = np.array(pf)
        self.type = type
        self.S = self.calc_S()
    

    def calc_S(self):
        S = np.array([self.kva[0]*(cos(self.pf[0])+j*sin(self.pf[0])), self.kva[1]*(cos(self.pf[1])+j*sin(self.pf[1])), self.kva[2]*(cos(self.pf[2])+j*sin(self.pf[2]))])
        return S




class Generator:
    """
    Class to represent generator objects
    """
    def __init__(self, name: str, bus: str, voltage: float, real_power: float, sub_transient_reactance = 0.0, neg_impedance = 0.0, zero_impedance = 0.0, gnd_impedance = None, var_limit = float('inf')):
        """
        Constructor for Generator class
        :param name: Name of generator
        :param bus: Bus connection
        :param voltage: Generator voltage
        :param real_power: Real power generation
        :param sub_transient_reactance: Sub transient reactance for faults
        :param neg_impedance: Negative impedance
        :param zero_impedance: Zero impedance
        :param gnd_impedance: Ground impedance
        :param var_limit: VAR Limit
        """
        self.name = name
        self.bus = bus
        self.voltage = voltage
        self.real_power = real_power*1e6
        self.reactive_power = 0.
        self.X0 = self.calc_X0(zero_impedance)
        self.X1 = self.calc_X1(sub_transient_reactance)
        self.X2 = self.calc_X2(neg_impedance)
        self.Zn = gnd_impedance
        self.Y0prim = self.calc_Y0prim() if self.X0 != 0.0 else 0.0
        self.var_limit = var_limit
    

    def calc_X0(self, X0):
        """
        Return imaginary reactance X0
        :param X0: Reactance
        :return:
        """
        X0 = 1j*X0*settings.powerbase/self.real_power
        return X0


    def calc_X1(self, X1):
        """
        Return imaginary reactance X1
        :param X1: Reactance
        :return:
        """
        X1 = 1j*X1*settings.powerbase/self.real_power
        return X1


    def calc_X2(self, X2):
        """
        Return imaginary reactance X2
        :param X2: Reactance
        :return:
        """
        X2 = 1j*X2*settings.powerbase/self.real_power
        return X2


    def set_power(self, real: float, reactive: float):
        """
        Set function for power
        :param real: Real power
        :param reactive: Reactive power
        :return:
        """
        self.real_power = real*1e6
        self.reactive_power = reactive*1e6

        
    def calc_Y0prim(self):
        """
        Generator primitive admittance matrix
        :return:
        """
        if self.Zn == None:
            Y0prim = 0
        elif self.Zn >= 0:
            Y0prim = 1/(3*self.Zn+self.X0)
        
        return Y0prim


# validation tests
if __name__ == '__main__':
    from Component import Load
    load1 = Load("Load1", "Bus3", [1275, 1800, 2375], [0.85, 0.9, 0.95])
    print(load1.kva)
    print(load1.pf)
    print(load1.S)