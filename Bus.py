"""
Module to implement bus functionality

Filename: Bus.py
Author: Justin Lipner
Date: 2025-01-23
"""


class Bus:
    """
    Bus class supporting unbalanced, multiphase operation.
    """
    def __init__(self, name: str, base_kv: float, index: int, phases=None):
        """
        :param name:   Name of bus
        :param base_kv: Base voltage of bus (line-to-neutral) in kV
        :param index:  Index of bus
        :param phases: Iterable of phase labels, e.g. ['A','B','C']
        """
        self.name = name
        self.base_kv = base_kv * 1e3  # convert to V
        self.index = index

        # default to three-phase
        self.phases = ['A','B','C'] if phases is None else phases

        # initialize per-phase quantities
        self.Vpu   = {phase: 1.0     for phase in self.phases}
        self.V     = {phase: self.base_kv for phase in self.phases}
        self.angle = {phase: 0.0     for phase in self.phases}
        self.P     = {phase: 0.0     for phase in self.phases}  # real power injection (W)
        self.Q     = {phase: 0.0     for phase in self.phases}  # reactive power injection (var)
        self.type  = "PQ"  # you might also want per-phase types


# validation tests
if __name__ == '__main__':
    from Bus import Bus
    bus1 = Bus("Bus 1", 20, 1)
    print(bus1.name, bus1.base_kv, bus1.index, bus1.V)