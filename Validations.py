import numpy as np

from Circuit import Circuit
from Settings import settings
j = 1j

def Create4NodeSystem():
  circ = Circuit("4nodebus")
  print(circ.name)

  circ.add_bus("bus1", 1)
  circ.add_bus("bus2", 1)
  circ.add_bus("bus3", 1)
  circ.add_bus("bus4", 1)

  phase_conductor = circ.add_conductor("1/0_ACSR", 0.398, 0.00446, 1.12, 230)
  print(circ.conductors["1/0_ACSR"].GMR)
  neutral_conductor = circ.add_conductor("1/0_ACSR_neutral", 0.398, 0.00446, 1.12, 230)
  #geometry1 = circ.add_geometry("Geometry 1", [0+j*29, 7+j*29, 2.5+j*29, 4+j*25], 3, phase_conductor, neutral_conductor)
  #line1 = circ.add_dline_from_geometry("OH1", "bus1", "bus2", geometry1, 1.893939)
