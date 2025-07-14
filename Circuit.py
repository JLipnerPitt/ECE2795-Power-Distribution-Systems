"""
Module to implement circuit/system functionality
Disclaimer: ChatGPT used for assistance

Filename: Circuit.py
Author: Justin Lipner, Bailey Stout
Date: 2025-02-03
"""

from Component import Load, Generator
import numpy as np
from Bus import Bus
from DistributionLine import DistributionLine
from Geometry import Geometry
from Conductor import Conductor
from Transformer import Transformer
from Settings import settings
from math import sin, cos
import pandas as pd

#  This class "creates" circuits.
class Circuit:
    """
    Circuit class to hold information about system
    """
    def __init__(self, name: str):
        """
        Constructor for the circuit class
        :param name: Name of circuit
        """
        self.name = name
        self.powerbase = settings.powerbase

        self.buses = {}
        self.conductors = {}
        self.geometries = {}
        self.distribution_lines = {}
        self.transformers = {}
        self.loads = {}
        self.generators = {}
        self.reactors = {}
        self.capacitors = {}

        self.count = 0
        self.slack_bus = str
        self.slack_index = int
        self.pq_indexes = []
        self.pv_indexes = []
        self.pq_and_pv_indexes = []
        self.bus_order = []

        self.voltages = None
        self.currents = None
        
        self.changed = False


    def change_power_base(self, p: float):
        """
        Changes the systems base power.
        :param p: New base power.
        :return:
        """
        settings.set_powerbase(p)
        self.powerbase = p*1e6


    def change_frequency(self, f: float):
        """
        Changes the systems base frequency.
        :param f: New system frequency.
        :return:
        """
        settings.set_freq(f)

    
    def add_bus(self, name: str, voltage: float):
        """
        Adds a bus to the system.
        :param name: Name of base
        :param voltage: Rated voltage of bus
        :return:
        """
        if name in self.buses:
            print(f"{name} already exists. No changes to circuit")

        else:
            self.count += 1
            bus = Bus(name, voltage, self.count)
            self.buses.update({name: bus})
            self.pq_indexes.append(self.count)
            self.bus_order.append(self.count)


    def add_load(self, name: str, bus: str, real: list[float], pf: list[float], type: str = 'PQ', connection: str = 'Y', phases=None):
        """
        Adds a load to system.
        :param name: Name of load
        :param bus: Bus connection
        :param real: Load's real power usage
        :param reactive: Load's reactive power usage
        :return:
        """
        if name in self.loads:
            print(f"{name} already exists. No changes to circuit.")
            return
        
        if bus not in self.buses:
            print(f"{bus} does not exist. No changes to circuit.")
            return

        load = Load(name, bus, real, pf, type, connection, phases)
        self.loads.update({name: load})
        #self.buses[bus].set_power(-real*1e6, -reactive*1e6)


    '''
    def add_dline_from_parameters(self, name: str, bus1: str, bus2: str, R: float, X: float, B: float):
        """
        Adds a transmission line to system.
        :param name: Name of transmission line
        :param bus1: First bus connection
        :param bus2: Second bus connection
        :param R: per unit resistance
        :param X: per unit reactance
        :param B: per unit shunt admittance
        :return:
        """
        
        if name in self.distribution_lines:
            print(f"{name} already exists. No changes to circuit")
            return
        
        tline = DistributionLine.from_parameters(name, self.get_bus(bus1), self.get_bus(bus2), R, X, B)
        self.transmission_lines.update({name: tline})
        self.changed = True
    '''
    

    def add_generator(self, name: str, bus: str, voltage: float, real_power: float, pos_imp = 0.0, neg_imp = 0.0, zero_imp = 0.0, gnd_imp = 0.0, var_limit = float('inf')):
        """
        Adds a generator to system.
        :param name: Name of transformer
        :param bus: Bus connection
        :param voltage: Operating voltage in pu
        :param real_power: Power rating
        :param pos_imp: Positive sequence impedance
        :param neg_imp: Negative sequence impedance
        :param zero_imp: Zero sequence impedance
        :param gnd_imp: Ground sequence impedance
        :param var_limit: Maximum VARs the generator can safely output
        :return:
        """
        if name in self.generators:
            print(f"{name} already exists. No changes to circuit")
            return
        
        if bus not in self.buses:
            print(f"{bus} does not exist. No changes to circuit")
            return
    
        if len(self.generators) == 0:
            gen = Generator(name, bus, voltage, real_power, pos_imp, neg_imp, zero_imp, gnd_imp, var_limit)
            self.generators.update({name: gen})
            self.buses[bus].type = "Slack"
            self.slack_bus = bus
            self.slack_index = self.buses[bus].index
            self.pq_indexes.remove(self.buses[bus].index)
            self.buses[bus].set_power(real_power*1e6, 0)
        
        else:
            gen = Generator(name, bus, voltage, real_power, pos_imp, neg_imp, zero_imp, gnd_imp, var_limit)
            self.generators.update({name: gen})
            self.buses[bus].type = "PV"
            self.pq_indexes.remove(self.buses[bus].index)
            self.pv_indexes.append(self.buses[bus].index)
            self.buses[bus].set_power(real_power*1e6, 0)


    def add_conductor(self, name: str, diam: float, GMR: float, resistance: float, ampacity: float):
        """
        Adds conductor to circuit object for repeated use
        :param name: Name of conductor
        :param diam: Diameter of conductor in inches
        :param GMR: GMR of conductor in feet
        :param resistance: Resistance of conductor at 50°C, 60 Hz
        :param ampacity: Rated ampacity of conductor
        :return:
        """
        if name in self.conductors:
            print(f"{name} already exists. No changes to circuit")

        else:
            conductor = Conductor(name, diam, GMR, resistance, ampacity)
            self.conductors.update({name: conductor})


    def add_geometry(self, name: str, d: list[complex], phase_conductor: Conductor, neutral_conductor: Conductor, phases=[1, 1, 1]):
        """
        Adds geometry to circuit object for repeated use
        :param name: Name of geometry
        :param d: Coordinates of each conductor
        :param nphases: Number of phases
        :return:
        """
        if name in self.geometries:
            print("Name already exists. No changes to circuit")
    
        else:
            geometry = Geometry(name, d, phase_conductor, neutral_conductor, phases)
            self.geometries.update({name: geometry})


    def add_dline_from_geometry(self, name: str, bus1: str, bus2: str, geometry: str, length: float):
        """
        Adds a transmission line to system.
        :param name: Name of transmission line
        :param bus1: First bus connection
        :param bus2: Second bus connection
        :param bundle: Bundle information passed via subclass
        :param geometry: Geometry information passed via subclass
        :param length: Length of transmission line in miles
        :return:
        """
        
        if name in self.distribution_lines:
            print(f"{name} already exists. No changes to circuit")
            return
        
        dline = DistributionLine(name, self.get_bus(bus1), self.get_bus(bus2), self.get_geometry(geometry), length)
        self.distribution_lines.update({name: dline})
        self.changed = True


    def add_transformer(self, name: str, bus1: str, bus2: str, Vprim: float, Vsec: float, power_rating: list[float],
                 resistance_percent: list[float], reactance_percent: list[float]):
        
        if name in self.transformers:
            print(f"{name} already exists. No changes to circuit")
            return

        transformer = Transformer(name, bus1, bus2, Vprim, Vsec, power_rating, resistance_percent, reactance_percent)
        self.transformers.update({name: transformer})
        self.changed = True


    def get_conductor(self, name: str):
        """
        Retrieves the name of the specified conductor.
        :param name: Conductor name
        :return:
        """
        return self.conductors[name]
    

    def get_bus(self, name: str):
        """
        Retrieves the name of the specified bus.
        :param name: Bus name
        :return:
        """
        return self.buses[name]
    

    def get_geometry(self, name: str):
        """
        Retrieves the name of the specified geometry.
        :param name: Geometry name
        :return:
        """
        return self.geometries[name]


    def get_line_impedance(self, line: str):
        return self.distribution_lines[line].Zabc
    

    def get_line_shunt_admittance(self, line: str):
        return self.distribution_lines[line].Yabc


    def change_slack(self, old: str, new: str):
        """
        Changes the system's slack bus.
        :param old: Old slack bus.
        :param new: New slack bus.
        :return:
        """
        if self.buses[new].type != "PV":
            print(f"Cannot make '{self.buses[new].name}' a slack bus because it has no generator connection. No changes made to circuit.")
            return
        
        self.buses[old].set_type("PV")
        self.buses[new].set_type("Slack")
        self.slack_bus = new
        self.slack_index = self.buses[new].index
        self.pv_indexes.remove(self.buses[new].index)
        self.pv_indexes.append(self.buses[old].index)
    
    
    def do_fbsweep(self):
        from Solution import LIT
        solution = LIT(self)
        self.voltages, self.currents = solution.lit()
        self.print_data()
    

    def update_voltages_and_angles(self):
        """
        Updates the voltages and angles at each bus with the values calculated in the power flow results.
        :return:
        """
        d = self.x[self.x.index.str.startswith("d")]
        V = self.x[self.x.index.str.startswith("V")]

        for bus in self.buses:
            index = self.buses[bus].index-1
            self.buses[bus].set_bus_v(V.iloc[index, 0])
            self.buses[bus].set_angle(d.iloc[index, 0])

            
    def print_data(self):
        """
        Prints necessary information from system.
        :return:
        """
        for i in range(len(self.buses)):
            print(f"[VLGabc]{i+1} =", np.abs(self.voltages[f"V{i+1}"]))
        
        print()
        
        for i in range(len(self.buses)):
            print(f"[Iabc]{i+1} =", np.abs(self.currents[f"I{i+1}"]))


    

# validation tests
if __name__ == '__main__':
    
    import Validations
    #Validations.CreateProject1()
    #Validations.CreateProject2()
    Validations.CreateProject4()