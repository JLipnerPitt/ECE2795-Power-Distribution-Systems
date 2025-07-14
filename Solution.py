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
    
    for i in range(self.numbuses):
      self.Iabc.update({f"I{i+1}": np.array([0, 0, 0], dtype=complex)})

    for bus in self.circ.buses.values():
      V0 = bus.base_kv/sqrt(3)
      self.VLGabc.update({f"V{bus.index}": np.array([V0, V0*np.exp(j*-2*np.pi/3), V0*np.exp(j*2*np.pi/3)], dtype=complex)})


  def lit(self):
    self.setup()
    iters = 200
    for i in range(iters):
      self.forward_sweep_project()
      self.backward_sweep_project()
    
    return self.VLGabc, self.Iabc


  def forward_sweep_project(self):
    VLGabc2 = np.matmul(self.A["line1"], self.VLGabc["V1"]) - np.matmul(self.B["line1"], self.Iabc["I2"])
    VLGabc3 = np.matmul(self.circ.transformers["T1"].At, self.VLGabc["V2"]) - np.matmul(self.circ.transformers["T1"].Bt, self.Iabc["I3"])
    VLGabc4 = np.matmul(self.A["line2"], self.VLGabc["V3"]) - np.matmul(self.B["line2"], self.Iabc["I4"])
    self.VLGabc["V2"] = VLGabc2
    self.VLGabc["V3"] = VLGabc3
    self.VLGabc["V4"] = VLGabc4

  
  def backward_sweep_project(self):
    # calculating bus 4 load current
    load1 = self.circ.loads["load1"]
    S = load1.S                                  
    VLGabc4 = self.VLGabc["V4"]                  
    I4 = np.conjugate(S/VLGabc4)
    self.Iabc["I4"] = I4

    # calculating I3, current in line2
    I3 = (self.c["line2"] @ self.VLGabc["V3"]
          + self.d["line2"] @ I4)
    self.Iabc["I3"] = I3
    
    # calculating I2, current into the primary side of the transformer
    I2 = np.matmul(self.circ.transformers["T1"].dt, I3)
    self.Iabc["I2"] = I2

    # build I1 from line1 using I2
    I1 = (self.c["line1"] @ self.VLGabc["V2"]
          + self.d["line1"] @ I2)
    self.Iabc["I1"] = I1
