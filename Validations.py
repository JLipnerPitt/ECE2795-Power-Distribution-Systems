import numpy as np

from Circuit import Circuit
from Settings import settings
j = 1j

def Create4NodeSystem():
  circ = Circuit("4nodebus")
  print(circ.name)

  circ.add_bus("bus1", 12.47)
  circ.add_bus("bus2", 12.47)
  #circ.add_bus("bus3", 12.47)

  circ.add_conductor("1/0_ACSR", 0.398, 0.00446, 1.12, 230)
  circ.add_conductor("1/0_ACSR_neutral", 0.398, 0.00446, 1.12, 230)
  circ.add_geometry("Geometry 1", [0+j*29, 7+j*29, 2.5+j*29, 4+j*25], 3, circ.conductors["1/0_ACSR"], circ.conductors["1/0_ACSR_neutral"])
  line1 = circ.add_dline_from_geometry("line1", "bus1", "bus2", "Geometry 1", 0.3787879)
  #line2 = circ.add_dline_from_geometry("line2", "bus2", "bus3", "Geometry 1", 0.4734848)
  circ.add_load("load1", "bus3", [1800, 1800, 1800], [0.9, 0.9, 0.9])
  
  circ.do_fbsweep()


def CreateProject1():
  circ = Circuit("Project1")
  print(circ.name)

  circ.add_bus("bus1", 12.47)
  circ.add_bus("bus2", 12.47)
  circ.add_conductor("1/0_ACSR", 0.398, 0.00446, 1.12, 230)
  circ.add_conductor("1/0_ACSR_neutral", 0.398, 0.00446, 1.12, 230)
  circ.add_geometry("Geometry 1", [0+j*29, 7+j*29, 2.5+j*29, 4+j*25], 3, circ.conductors["1/0_ACSR"], circ.conductors["1/0_ACSR_neutral"])
  line1 = circ.add_dline_from_geometry("line1", "bus1", "bus2", "Geometry 1", 0.3787879)
  circ.add_load("load1", "bus2", [1000, 800, 1200], [0.9, 0.85, 0.95])
  circ.do_fbsweep()