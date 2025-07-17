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
    self.VLGabc = {}
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

    self.VLGabc.update({"V1": (1e3/sqrt(3))*np.array([12.47, 12.47*np.exp(j*-2*np.pi/3), 12.47*np.exp(j*2*np.pi/3)], dtype=complex)})
    self.VLGabc.update({"V2": (1e3/sqrt(3))*np.array([12.47, 12.47*np.exp(j*-2*np.pi/3), 12.47*np.exp(j*2*np.pi/3)], dtype=complex)})
    self.VLGabc.update({"V3": 1e3*np.array([4.16, 4.16*np.exp(j*-2*np.pi/3), 4.16*np.exp(j*2*np.pi/3)], dtype=complex)})
    self.VLGabc.update({"V4": 1e3*np.array([4.16, 4.16*np.exp(j*-2*np.pi/3), 4.16*np.exp(j*2*np.pi/3)], dtype=complex)})


  def lit(self):
    self.setup()
    iters = 500
    for i in range(iters):
      self.forward_sweep()
      self.backward_sweep()
    
    return self.VLGabc, self.Iabc


  def forward_sweep(self):
    VLGabc2 = np.matmul(self.A["line1"], self.VLGabc["V1"]) - np.matmul(self.B["line1"], self.Iabc["I12"])
    VLGabc3 = np.matmul(self.circ.transformers["T1"].At, self.VLGabc["V2"]) - np.matmul(self.circ.transformers["T1"].Bt, self.Iabc["I23"])
    VLGabc4 = np.matmul(self.A["line2"], self.VLGabc["V3"]) - np.matmul(self.B["line2"], self.Iabc["I34"])
    self.VLGabc["V2"] = VLGabc2
    self.VLGabc["V3"] = VLGabc3
    self.VLGabc["V4"] = VLGabc4

  
  def backward_sweep(self):
    # calculating Eload current
    load1 = self.circ.loads["load1"]
    S = load1.S                                  
    VLGabc4 = self.VLGabc["V4"]                  
    Iload = np.conjugate(S/VLGabc4)

    # calculating I34
    I34 = (self.c["line2"] @ self.VLGabc["V4"]
          + self.d["line2"] @ Iload)
    self.Iabc["I34"] = I34
    
    # calculating IABC, current into the primary side of the transformer
    IABC = np.matmul(self.circ.transformers["T1"].dt, I34)
    self.Iabc["I23"] = IABC

    # build I12 from line1 using I23
    I12 = (self.c["line1"] @ self.VLGabc["V2"]
          + self.d["line1"] @ IABC)
    self.Iabc["I12"] = I12

