from Circuit import Circuit
from math import sqrt, sin, cos
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

    V0 = self.circ.buses[self.circ.slack_bus].base_kv/sqrt(3)
    for i in range(self.numbuses):
      self.VLGabc.update({f"V{i+1}": np.array([V0, V0*np.exp(j*-2*np.pi/3), V0*np.exp(j*2*np.pi/3)], dtype=complex)})
    
  
  def lit1(self):
    self.setup()
    iters = 4
    for i in range(iters):
      self.forward_sweep_project1()
      self.backward_sweep_project1()
    
    return self.VLGabc, self.Iabc


  def lit2(self):
    self.setup()
    iters = 4
    for i in range(iters):
      self.forward_sweep_project2()
      self.backward_sweep_project2()
    
    return self.VLGabc, self.Iabc
  
  
  def forward_sweep_project1(self):
    VLGabc2 = np.matmul(self.A["line1"], self.VLGabc["V1"]) - np.matmul(self.B["line1"], self.Iabc["I2"])
    self.VLGabc["V2"] = VLGabc2


  def forward_sweep_project2(self):
    VLGabc2 = np.matmul(self.A["line1"], self.VLGabc["V1"]) - np.matmul(self.B["line1"], self.Iabc["I2"])
    VLGabc3 = np.matmul(self.A["line2"], self.VLGabc["V2"]) - np.matmul(self.B["line2"], self.Iabc["I3"])
    self.VLGabc["V2"] = VLGabc2
    self.VLGabc["V3"] = VLGabc3


  def backward_sweep_project1(self):
    # bus 2 (PQ load on all 3 phases)
    load2  = self.circ.loads["load1"]
    S2 = load2.S
    V2 = self.VLGabc["V2"]
    Iload2 = np.conjugate(S2 / V2) # vector of 3 currents

    # build I1 from line1 using I2
    I1 = ( self.c["line1"] @ self.VLGabc["V2"]
          + self.d["line1"] @ Iload2 ) + Iload2

    self.Iabc["I1"] = I1 
    self.Iabc["I2"] = Iload2

  
  def backward_sweep_project2(self):
    # bus 3 (constant‐Z, phase B only)
    load3 = self.circ.loads["load2"]
    Sb = load3.S[1]                                  
    Vph = self.circ.buses[self.circ.slack_bus].base_kv/np.sqrt(3)
    Zb = Vph**2 / Sb                             
    I3 = np.array([0, self.VLGabc["V3"][1] / Zb, 0], dtype=complex)

    # bus 2 (PQ load on all 3 phases)
    load2  = self.circ.loads["load1"]
    S2 = load2.S
    V2 = self.VLGabc["V2"]
    Iload2 = np.conjugate(S2 / V2) # vector of 3 currents

    # build I2 by adding line-shunt + load2 + I3 injection
    I2 = ( self.c["line2"] @ self.VLGabc["V3"]
          + self.d["line2"] @ I3
          + Iload2 )

    # build I1 from line1 using I2
    I1 = ( self.c["line1"] @ self.VLGabc["V2"]
          + self.d["line1"] @ I2 )

    self.Iabc["I1"], self.Iabc["I2"], self.Iabc["I3"] = I1, I2, I3
