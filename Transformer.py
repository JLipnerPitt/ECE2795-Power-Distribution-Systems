"""
Module to implement transformer functionality

Filename: Transformer.py
Author: Justin Lipner
Date: 2025-02-03
"""

import numpy as np
from Bus import Bus
from math import atan, sin, cos
from Settings import settings


class Transformer:
    """
    Transformer class to hold transformer information
    """

    def __init__(self, name: str, connection: str, bus1: Bus, bus2: Bus, Vprim: float, Vsec: float, power_rating: float,
                 impedance_percent: float, x_over_r_ratio: float, gnd_impedance=None):
        """
        Constructor for Transformer objects
        :param name: Name of transformer
        :param bus1: First bus connection
        :param bus2: Second bus connection
        :param power_rating: Power rating
        :param impedance_percent: Impedance percent
        :param x_over_r_ratio: X/R Ratio
        :param gnd_impedance: Impedance that grounds the Wye side
        """
        self.name = name
        self.type = connection
        self.bus1 = bus1
        self.bus2 = bus2
        self.Vprim = Vprim*1e3
        self.Vsec = Vsec*1e3
        self.n = self.Vprim/self.Vsec
        self.power_rating = power_rating*1e6
        self.impedance_percent = impedance_percent
        self.Ztabc = self.calc_Zt()
  

    def calc_Zt(self):
        #X = X*settings.powerbase/self.power_rating  # updating pu to system power base
        Z = self.Vprim**2/self.power_rating 
        Ztabc = np.array([Z, Z, Z], dtype=complex)
        return Ztabc




# validation tests 
if __name__ == '__main__':
    from Transformer import Transformer
    from Bus import Bus
    from Settings import settings

    settings.set_powerbase(100e6)
    bus1 = Bus("bus1", 15e3, 1)
    bus2 = Bus("bus2", 30e3, 2)
    power_rating = 125e6
    impedance_percent = 8.5
    x_over_r_ratio = 10
    transformer1 = Transformer("T1", bus1, bus2, power_rating, impedance_percent, x_over_r_ratio)

    print(f"Name: {transformer1.name}, from {transformer1.bus1.name} to {transformer1.bus2.name},", 
          f"Rating = {transformer1.power_rating/1e6} MVA")
    print(f"Z = {transformer1.Zpu}, Y = {transformer1.Ypu}")
    print(f"Yprim = {transformer1.yprim}")
