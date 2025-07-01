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
    self.Iload = {}
    self.Vload = {}
    self.VLGabc = {}
    self.Iabc = {}

  
  def setup(self):
    for line in self.circ.distribution_lines.values():
      self.a.update({line.name: np.identity(3) + 0.5*np.matmul(line.Zabc, line.Yabc)})
      self.b.update({line.name: line.Zabc})
      self.c.update({line.name: line.Yabc + 0.25*np.matmul(line.Yabc, np.matmul(line.Zabc, line.Yabc))})
      self.d.update({line.name: np.identity(3) + 0.5*np.matmul(line.Yabc, line.Zabc)})
    
    self.Iload.update({"Im": np.array([0, 0, 0])})
    self.Vload.update({"Vm": np.array([12.47, 12.47*np.exp(j*-2*np.pi/3), 12.47*np.exp(j*2*np.pi/3)])})

    for i in range(self.numbuses):
      self.VLGabc.update({f"V{i+1}": 1e3*np.array([12.47, 12.47*np.exp(j*-2*np.pi/3), 12.47*np.exp(j*2*np.pi/3)])/sqrt(3)})
    
    for i in range(len(self.circ.distribution_lines)+1):
      self.Iabc.update({f"Iline{i+1}": np.array([0, 0, 0])})
    
    
  def lit(self):
    self.setup()
    tol = 0.0001
    iters = 50
    for i in range(iters):
      self.forward_sweep()
      self.backward_sweep()
    
    print(self.Iabc)
    for v in self.VLGabc.values():
      print(np.abs(v))

      
  def forward_sweep(self):
    for i in range(len(self.circ.distribution_lines)):
      A = np.linalg.inv(self.a[f"line{i+1}"])
      B = np.matmul(np.linalg.inv(self.a[f"line{i+1}"]), self.b[f"line{i+1}"])
      VLGabc = np.matmul(A, self.VLGabc[f"V{i+1}"]) - np.matmul(B, self.Iabc[f"Iline{i+2}"])
      self.VLGabc.update({f"V{i+2}": VLGabc})


  def backward_sweep(self):
    n = len(self.circ.distribution_lines)
    kva = self.circ.loads["load1"].kva
    pf = self.circ.loads["load1"].pf
    Sload = np.array([kva[0]*(cos(pf[0])+j*sin(pf[0])), kva[1]*(cos(pf[1])+j*sin(pf[1])), kva[2]*(cos(pf[2])+j*sin(pf[2]))])
    Iload = np.conjugate(np.divide(Sload, self.VLGabc["V2"]))
    self.Iabc["Iline2"] = Iload

    for i in range(len(self.circ.distribution_lines)-1):
      Iabc = np.matmul(self.c[f"line{i+1}"], self.VLGabc[f"V{i+2}"]) + np.matmul(self.d[f"line{i+1}"], self.Iabc[f"Iline{i+2}"])
      self.Iabc.update({f"Iline{i+1}": Iabc})