import numpy as np

from Circuit import Circuit
from Settings import settings
j = 1j

def CreateProject4():
  circ = Circuit("Project4")
  print(circ.name)

  circ.add_bus("bus1", 12.47)
  circ.add_bus("bus2", 12.47)
  circ.add_bus("bus3", 4.16)
  circ.add_bus("bus4", 4.16)
  circ.add_generator("Gen1", "bus1", 12.47, 20)

  circ.add_conductor("336400_26/7_ACSR", 0.721, 0.0244, 0.306, 230)
  circ.add_conductor("4/0_6/1_ACSR", 0.563, 0.00814, 0.592, 230)
  circ.add_geometry("Geometry 1", [0+j*29, 2.5+j*29, 7+j*29, 4+j*25], circ.conductors["336400_26/7_ACSR"], circ.conductors["4/0_6/1_ACSR"], [0, 1, 0])
  circ.add_dline_from_geometry("line1", "bus1", "bus2", "Geometry 1", 0.3787879)
  circ.add_dline_from_geometry("line2", "bus3", "bus4", "Geometry 1", 0.4734848)
  circ.add_transformer("T1", "bus2", "bus3", 7.20, 2.40, [6000, 6000, 6000], [0.01, 0.01, 0.01], [0.06, 0.06, 0.06])
  circ.add_load("load1", "bus4", [1800, 1800, 1800], [0.9, 0.9, 0.9])

  circ.do_fbsweep()
  