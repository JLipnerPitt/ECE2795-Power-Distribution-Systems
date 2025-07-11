import numpy as np

from Circuit import Circuit
from Settings import settings
j = 1j


def CreateProject1():
  circ = Circuit("Project1")
  print(circ.name)

  circ.add_bus("bus1", 12.47)
  circ.add_generator("Gen1", "bus1", 12.47, 20)
  circ.add_bus("bus2", 12.47)
  circ.add_conductor("1/0_ACSR", 0.398, 0.00446, 1.12, 230)
  circ.add_conductor("1/0_ACSR_neutral", 0.398, 0.00446, 1.12, 230)
  circ.add_geometry("Geometry 1", [0+j*29, 7+j*29, 2.5+j*29, 4+j*25], circ.conductors["1/0_ACSR"], circ.conductors["1/0_ACSR_neutral"])
  circ.add_dline_from_geometry("line1", "bus1", "bus2", "Geometry 1", 1.893939)
  circ.add_load("load1", "bus2", [1000, 800, 1200], [0.9, 0.85, 0.95])
  circ.do_fbsweep1()
  print()


def CreateProject2():
  circ = Circuit("Project2")
  print(circ.name)

  circ.add_bus("bus1", 12.47)
  circ.add_bus("bus2", 12.47)
  circ.add_bus("bus3", 12.47)
  circ.add_generator("Gen1", "bus1", 12.47, 20)

  circ.add_conductor("1/0_ACSR", 0.398, 0.00446, 1.12, 230)
  circ.add_conductor("1/0_ACSR_neutral", 0.398, 0.00446, 1.12, 230)
  circ.add_geometry("Geometry 1", [0+j*29, 7+j*29, 2.5+j*29, 4+j*25], circ.conductors["1/0_ACSR"], circ.conductors["1/0_ACSR_neutral"], [1, 1, 1])
  circ.add_dline_from_geometry("line1", "bus1", "bus2", "Geometry 1", 1.893939)

  circ.add_conductor("336400_26/7_ACSR", 0.721, 0.0244, 0.306, 230)
  circ.add_conductor("4/0_6/1_ACSR", 0.563, 0.00814, 0.592, 230)
  circ.add_geometry("Geometry 2", [0+j*29, 7+j*29, 2.5+j*29, 4+j*25], circ.conductors["336400_26/7_ACSR"], circ.conductors["4/0_6/1_ACSR"], [0, 1, 0])
  circ.add_dline_from_geometry("line2", "bus2", "bus3", "Geometry 2", 1.893939)

  circ.add_load("load1", "bus2", [1000, 800, 1200], [0.9, 0.85, 0.95])
  circ.add_load("load2", "bus3", [0, 200, 0], [0, 0.9, 0], type='Z')
  circ.do_fbsweep2()
  print()


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
  circ.add_geometry("Geometry 1", [0+j*29, 7+j*29, 2.5+j*29, 4+j*25], circ.conductors["336400_26/7_ACSR"], circ.conductors["4/0_6/1_ACSR"], [0, 1, 0])
  circ.add_dline_from_geometry("line1", "bus1", "bus2", "Geometry 1", 0.3787879)
  circ.add_dline_from_geometry("line2", "bus3", "bus4", "Geometry 1", 0.4734848)
  circ.add_transformer("T1", "bus2", "bus3", 7.20, 2.40, [6000, 6000, 6000], [0.01, 0.01, 0.01], [0.06, 0.06, 0.06])
  circ.add_load("load1", "bus4", [1800, 1800, 1800], [0.9, 0.9, 0.9])
  