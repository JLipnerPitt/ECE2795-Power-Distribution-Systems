from Circuit import Circuit
from math import sqrt
import numpy as np
j = 1j

class LIT:

  def __init__(self, circ: Circuit):
    self.circ = circ
    self.numbuses = len(self.circ.buses)
    self.a = {}
    self.b = {}
    self.c = {}
    self.d = {}
    self.A = {}
    self.B = {}
    self.Di = np.array([[1, 0, -1], [-1, 1, 0], [0, -1, 1]])
    self.W = (1/3)*np.array([[2, 1, 0], [0, 2, 1], [1, 0, 2]])
    self.Voltages = {}
    self.Iabc = {}

  
  def setup(self):
    for line in self.circ.distribution_lines.values():
      self.a.update({line.name: np.identity(3) + 0.5*np.matmul(line.Zabc, line.Yabc)})
      self.A.update({line.name: np.linalg.inv(self.a[line.name])})
      self.b.update({line.name: line.Zabc})
      self.B.update({line.name: np.matmul(np.linalg.inv(self.a[line.name]), self.b[line.name])})
      self.c.update({line.name: line.Yabc + 0.25*np.matmul(line.Yabc, np.matmul(line.Zabc, line.Yabc))})
      self.d.update({line.name: np.identity(3) + 0.5*np.matmul(line.Yabc, line.Zabc)})
    
    self.Iabc.update({"I12": np.array([0, 0, 0], dtype=complex)})
    self.Iabc.update({"I23": np.array([0, 0, 0], dtype=complex)})
    self.Iabc.update({"I34": np.array([0, 0, 0], dtype=complex)})

    self.Voltages.update({"V1": (1e3/sqrt(3))*np.array([12.47, 12.47*np.exp(-j*2*np.pi/3), 12.47*np.exp(j*2*np.pi/3)], dtype=complex)})
    self.Voltages.update({"V2": (1e3/sqrt(3))*np.array([12.47, 12.47*np.exp(-j*2*np.pi/3), 12.47*np.exp(j*2*np.pi/3)], dtype=complex)})
    self.Voltages.update({"V3": (1e3/sqrt(3))*np.array([4.16, 4.16*np.exp(-j*2*np.pi/3), 4.16*np.exp(j*2*np.pi/3)], dtype=complex)})
    self.Voltages.update({"V4": (1e3/sqrt(3))*np.array([4.16, 4.16*np.exp(-j*2*np.pi/3), 4.16*np.exp(j*2*np.pi/3)], dtype=complex)})


  def lit(self):
    self.setup()
    iters = 15
    for i in range(iters):
      self.forward_sweep()
      self.backward_sweep()
    
    return self.Voltages, self.Iabc


  def forward_sweep(self):
    # 1) line 1:  V2 = A1·V1 – B1·I12
    VLGabc2 = self.A["line1"] @ self.Voltages["V1"] - self.B["line1"] @ self.Iabc["I12"]

    # 2) transformer: V3 = At·V2 – Bt·I23
    VLLabc3 = self.circ.transformers["T1"].At @ self.Voltages["V2"] - self.circ.transformers["T1"].Bt @ self.Iabc["I23"]

    # 3) line 2:  V4 = A2·V3 – B2·I34
    VLLabc4 = self.A["line2"] @ self.Voltages["V3"] - self.B["line2"] @ self.Iabc["I34"]

    # updating voltages
    self.Voltages["V2"] = VLGabc2
    self.Voltages["V3"] = VLLabc3
    self.Voltages["V4"] = VLLabc4

  
  def backward_sweep(self):
    # calculating delta load currents
    load1 = self.circ.loads["load1"]
    S = load1.S                                  
    VLLabc4 = self.Voltages["V4"]                  
    IDabc = np.conjugate(S/VLLabc4)

    # calculating secondary line currents
    Iabc = self.Di @ IDabc

    # calculating the equivalent secondary line to neutral voltages
    VLGabc4 = self.W @ VLLabc4

    # calculating I34
    I34 = (self.c["line2"] @ VLGabc4
          + self.d["line2"] @ Iabc)
    self.Iabc["I34"] = I34
    
    # calculating IABC, current into the primary side of the transformer
    IABC = np.matmul(self.circ.transformers["T1"].dt, Iabc)
    self.Iabc["I23"] = IABC

    # build I12 from line1 using I23
    I12 = (self.c["line1"] @ self.Voltages["V2"]
          + self.d["line1"] @ IABC)
    self.Iabc["I12"] = I12
    

    '''
    # 1) delta‐load branch currents from line-to-line volts
    Va, Vb, Vc    = self.VLGabc["V4"]
    V_LL          = np.array([Va-Vb, Vb-Vc, Vc-Va], dtype=complex)
    S             = self.circ.loads["load1"].S
    I_branch      = np.conjugate(S / V_LL)           # [Iab, Ibc, Ica]
    Iab, Ibc, Ica = I_branch
    I34           = np.array([Iab-Ica, Ibc-Iab, Ica-Ibc], dtype=complex)

    # 2) transformer primary current
    I23 = self.circ.transformers["T1"].dt @ I34

    # 3) upstream line1 current
    I12 = self.c["line1"] @ self.VLGabc["V2"] \
        + self.d["line1"] @ I23

    # store
    self.Iabc["I34"], self.Iabc["I23"], self.Iabc["I12"] = I34, I23, I12
    '''
