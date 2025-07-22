"""
Module to implement transformer functionality

Filename: Transformer.py
Author: Justin Lipner
Date: 2025-02-03
"""

import numpy as np
from math import sqrt
j = 1j


class Transformer:
    """
    Transformer class to hold transformer information
    """

    def __init__(self, name: str, bus1: str, bus2: str, VprimLL: float, VsecLL: float, power_rating: list[float],
                 resistance_percent: list[float], reactance_percent: list[float]):
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
        self.bus1 = bus1
        self.bus2 = bus2
        self.VprimLN = VprimLL*1e3/sqrt(3)
        self.VsecLL = VsecLL*1e3
        self.power_rating = np.array(power_rating)*1e3
        self.resistance_percent = resistance_percent
        self.reactance_percent = reactance_percent

        self.n = self.VprimLN/self.VsecLL
        self.Ztabc = self.calc_Zt()
        self.at = self.n*np.array([[1, -1, 0], [0, 1, -1], [-1, 0, 1]])
        self.bt = self.calc_bt()
        self.ct = np.zeros((3, 3))
        self.dt = 1/(self.n*3)*np.array([[1, -1, 0], [1, 2, 0], [-2, -1, 0]])
        self.At = 1/(self.n*3)*np.array([[2, 1, 0], [0, 2, 1], [1, 0, 2]])
        self.Bt = self.calc_Bt()


    def calc_Zt(self):
        #X = X*settings.powerbase/self.power_rating  # updating pu to system power base
        Zbase = [self.VsecLL**2/self.power_rating[0], self.VsecLL**2/self.power_rating[1], self.VsecLL**2/self.power_rating[2]]
        Zab = Zbase[0]*(self.resistance_percent[0] + j*self.reactance_percent[0])
        Zbc = Zbase[1]*(self.resistance_percent[1] + j*self.reactance_percent[1])
        Zca = Zbase[2]*(self.resistance_percent[2] + j*self.reactance_percent[2])
        Ztabc = np.diag(np.array([Zab, Zbc, Zca], dtype=complex))
        return Ztabc
    

    def calc_bt(self):
        bt1 = [self.Ztabc[0][0], -self.Ztabc[0][0], 0]
        bt2 = [self.Ztabc[1][1], 2*self.Ztabc[1][1], 0]
        bt3 = [-2*self.Ztabc[2][2], -self.Ztabc[2][2], 0]
        bt = (self.n/3)*np.array([bt1, bt2, bt3])
        return bt


    def calc_Bt(self):
        Bt1 = [2*self.Ztabc[0][0]+self.Ztabc[1][1], 2*(self.Ztabc[1][1]-self.Ztabc[0][0]), 0]
        Bt2 = [2*(self.Ztabc[1][1]-self.Ztabc[2][2]), 4*self.Ztabc[1][1]-self.Ztabc[2][2], 0]
        Bt3 = [self.Ztabc[0][0]-4*self.Ztabc[2][2], -self.Ztabc[0][0]-2*self.Ztabc[2][2], 0]
        Bt = (1/9)*np.array([Bt1, Bt2, Bt3])
        return Bt



# validation tests 
if __name__ == '__main__':
    from Transformer import Transformer

    transformer1 = Transformer("T1", "bus2", "bus3", 12.47, 0.24, [100, 50, 50], [0.01, 0.015, 0.015], [0.04, 0.035, 0.035])
    print(f"nt = {transformer1.n}")
    print()
    print(f"Ztabc = {transformer1.Ztabc}")
    print()
    print(f"at = {transformer1.at}")
    print()
    print(f"bt = {transformer1.bt}")
    print()
    print(f"dt = {transformer1.dt}")
    print()
    print(f"At = {transformer1.At}")
    print()
    print(f"Bt = {transformer1.Bt}")