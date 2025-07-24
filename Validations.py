from Circuit import Circuit
import numpy as np
j = 1j

def CreateProject4_Balanced_Loads():
  circ = Circuit("Project4_balanced")
  print(circ.name)

  circ.add_bus("bus1", 12.47)
  circ.add_bus("bus2", 12.47)
  circ.add_bus("bus3", 4.16)
  circ.add_bus("bus4", 4.16)
  circ.add_generator("Gen1", "bus1", 12.47, 20)

  circ.add_conductor("336400_26/7_ACSR", 0.721, 0.0244, 0.306, 230)
  circ.add_conductor("4/0_6/1_ACSR", 0.563, 0.00814, 0.592, 230)
  circ.add_geometry("Geometry 1", [0+j*29, 2.5+j*29, 7+j*29, 4+j*24], circ.conductors["336400_26/7_ACSR"], circ.conductors["4/0_6/1_ACSR"], [1, 1, 1])
  circ.add_dline_from_geometry("line1", "bus1", "bus2", "Geometry 1", 0.3787879)
  circ.add_dline_from_geometry("line2", "bus3", "bus4", "Geometry 1", 0.4734848, True)

  '''
  circ.distribution_lines["line1"].Zabc = 0.3787879*np.array([
    [0.4576 + 1.0780j, 0.1559 + 0.5017j, 0.1535 + 0.3849j],
    [0.1559 + 0.5017j, 0.4666 + 1.0482j, 0.1580 + 0.4236j],
    [0.1535 + 0.3849j, 0.1580 + 0.4236j, 0.4615 + 1.0651j]], dtype=complex)
  
  circ.distribution_lines["line2"].Zabc = 0.4734848*np.array([
    [0.4013 + 1.4133j, 0.0953 + 0.8515j, 0.0953 + 0.7266j],
    [0.0953 + 0.8515j, 0.4013 + 1.4133j, 0.0953 + 0.7802j],
    [0.0953 + 0.7266j, 0.0953 + 0.7802j, 0.4013 + 1.4133j]], dtype=complex)
  '''

  circ.add_transformer("T1", "bus2", "bus3", 12.47, 4.16, [6000, 6000, 6000], [0.01, 0.01, 0.01], [0.06, 0.06, 0.06])
  circ.add_load("load1", "bus4", [2000, 2000, 2000], [0.9, 0.9, 0.9])

  circ.do_fbsweep()
  print()


def CreateProject4_Unbalanced_Loads():
  circ = Circuit("Project4_unbalanced")
  print(circ.name)

  circ.add_bus("bus1", 12.47)
  circ.add_bus("bus2", 12.47)
  circ.add_bus("bus3", 4.16)
  circ.add_bus("bus4", 4.16)
  circ.add_generator("Gen1", "bus1", 12.47, 20)

  circ.add_conductor("336400_26/7_ACSR", 0.721, 0.0244, 0.306, 230)
  circ.add_conductor("4/0_6/1_ACSR", 0.563, 0.00814, 0.592, 230)
  circ.add_geometry("Geometry 1", [0+j*29, 2.5+j*29, 7+j*29, 4+j*24], circ.conductors["336400_26/7_ACSR"], circ.conductors["4/0_6/1_ACSR"], [1, 1, 1])
  circ.add_dline_from_geometry("line1", "bus1", "bus2", "Geometry 1", 0.3787879)
  circ.add_dline_from_geometry("line2", "bus3", "bus4", "Geometry 1", 0.4734848, True)
  circ.add_transformer("T1", "bus2", "bus3", 12.47, 4.16, [6000, 6000, 6000], [0.01, 0.01, 0.01], [0.06, 0.06, 0.06])
  circ.add_load("load1", "bus4", [1500, 2000, 2500], [0.85, 0.9, 0.95])

  #print(circ.distribution_lines["line1"].Zabc)

  circ.do_fbsweep()
  print()


def CreateExample8_3():
  from math import sqrt
  circ = Circuit("Exmaple 8.3")
  print(circ.name)

  circ.add_bus("bus1", 12.47)
  circ.add_bus("bus2", 240)
  circ.add_transformer("T1", "bus1", "bus2", 12.47, 0.240, [100, 50, 50], [0.01, 0.015, 0.015], [0.04, 0.035, 0.035])
  circ.add_load("load1", "bus2", [100, 50, 50], [0.9, 0.8, 0.8])

  Iabc = {}
  Voltages = {}
  Di = np.array([[1, 0, -1], [-1, 1, 0], [0, -1, 1]])
  Dv = np.array([[1, -1, 0], [0, 1, -1], [-1, 0, 1]])
  W = (1/3)*np.array([[2, 1, 0], [0, 2, 1], [1, 0, 2]])

  Iabc.update({"I12": np.array([0, 0, 0], dtype=complex)})
  Voltages.update({"V1": (1e3/sqrt(3))*np.array([12.47, 12.47*np.exp(-j*2*np.pi/3), 12.47*np.exp(j*2*np.pi/3)], dtype=complex)})
  Voltages.update({"V2": np.array([240, 240*np.exp(-j*2*np.pi/3), 240*np.exp(j*2*np.pi/3)], dtype=complex)})
  print("VLLabc =", np.abs(Voltages["V2"]))

  iters = 1
  VLLabc2 = Voltages["V2"]

  # calculating delta load currents
  load1 = circ.loads["load1"]
  S = load1.S
  print("SDabc =", S)                                                 
  IDabc = np.conjugate(S/Voltages["V2"])
  print("IDabc =", np.abs(IDabc))

  # calculating secondary line currents
  Iabc = Di @ IDabc
  print("Iabc =", np.abs(Iabc))

  # calculating the equivalent secondary line to neutral voltages
  VLNabc2 = W @ VLLabc2
  print("VLNabc2 =", np.abs(VLNabc2))

  VLNABC2 = circ.transformers["T1"].at @ VLNabc2 + circ.transformers["T1"].bt @ Iabc
  print("VLNABC2 =", np.abs(VLNABC2))

  VLLABC2 = Dv @ VLNABC2
  print("VLLABC2 =", np.abs(VLLABC2))

  # calculating IABC, current into the primary side of the transformer
  IABC = np.matmul(circ.transformers["T1"].dt, Iabc)
  print("IABC =", np.abs(IABC))
  print()


  