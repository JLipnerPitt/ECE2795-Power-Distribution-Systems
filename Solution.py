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
    self.Dv = np.array([[1, -1, 0], [0, 1, -1], [-1, 0, 1]])
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
    self.Iabc.update({"I34": np.array([0, 0, 0], dtype=complex)})

    self.Voltages.update({"V1": (1e3/sqrt(3))*np.array([12.47, 12.47*np.exp(j*-2*np.pi/3), 12.47*np.exp(j*2*np.pi/3)], dtype=complex)})
    self.Voltages.update({"V2": (1e3/sqrt(3))*np.array([12.47, 12.47*np.exp(j*-2*np.pi/3), 12.47*np.exp(j*2*np.pi/3)], dtype=complex)})
    self.Voltages.update({"V3": 1e3*np.array([4.16, 4.16*np.exp(j*-2*np.pi/3), 4.16*np.exp(j*2*np.pi/3)], dtype=complex)})
    self.Voltages.update({"V4": 1e3*np.array([4.16, 4.16*np.exp(j*-2*np.pi/3), 4.16*np.exp(j*2*np.pi/3)], dtype=complex)})


  def lit(self):
    self.setup()
    iters = 100
    for i in range(iters):
      self.backward_sweep()
      self.forward_sweep()
    
    return self.Voltages, self.Iabc


  def forward_sweep(self):
    VLNabc2 = self.A["line1"] @ self.Voltages["V1"] - self.B["line1"] @ self.Iabc["I12"]
    
    VLNabc3 = self.circ.transformers["T1"].At @ self.Voltages["V2"] - self.circ.transformers["T1"].Bt @ self.Iabc["I34"]

    VLNabc4 = self.A["line2"] @ self.Voltages["V3"] - self.B["line2"] @ self.Iabc["I34"] 

    self.Voltages["V2"] = VLNabc2
    self.Voltages["V3"] = VLNabc3
    self.Voltages["V4"] = VLNabc4

  
  def backward_sweep(self):
    # calculating bus 4 load current
    load1 = self.circ.loads["load1"]
    S = load1.S                                  
    VLLabc4 = self.Voltages["V4"]                  
    IDabc = np.conjugate(S/VLLabc4)
    Iabc = self.Di @ IDabc
    
    # calculating I12, current into the primary side of the transformer
    I12 = self.circ.transformers["T1"].dt @ Iabc

    # updating currents
    self.Iabc["I34"] = Iabc
    self.Iabc["I12"] = I12