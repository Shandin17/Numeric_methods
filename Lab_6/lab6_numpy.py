import numpy as np
from numpy import linalg as la
from numpy.linalg import inv
gamm=5.0/3.0
gm1g = (gamm - 1.0) / gamm
ggm1 =1.0 / gm1g
gp1g = (gamm + 1.0) / gamm
gamm1 = 1.0/gamm

rho_L = 1.0
v_L = 3.0
p_L = 1.0
rho_R = 2.0
v_R = 1.0
p_R = 1.0

def CalcW (rol, vl, pl, ror, vr, pr):
  wl1, wl2, wl3 = rol**0.5, rol**0.5*vl, rol**0.5*(0.5*vl**2.0 + ggm1*pl/rol)
  wr1, wr2, wr3 = ror**0.5, ror**0.5*vr, ror**0.5*(0.5*vr**2.0 + ggm1*pr/ror)
  w1, w2, w3 = 0.5*(wl1 + wr1), 0.5*(wl2 + wr2), 0.5*(wl3 + wr3)
  return w1, w2, w3

def CalcA (w1, w2, w3):
  B = np.array([[2.0*w1, 0.0, 0.0],[w2, w1, 0.0],[w3*gamm1, w2*gm1g, w1*gamm1]])
  C = np.array([[w2, w1, 0.0],[gm1g*w3, gp1g*w2, gm1g*w1],[0.0, w3, w2]])
  A = np.matmul(C,inv(B))
  return A

def FindEigenVects (rol, vl, pl, ror, vr, pr):
  w1,w2,w3 = CalcW (rol, vl, pl, ror, vr, pr)
  A = CalcA(w1,w2,w3)
  return np.linalg.eig(A)


Lambda = FindEigenVects(rho_L,v_L,p_L,rho_R,v_R,p_R)
print(Lambda)
temp = Lambda[0][1]
Lambda[0][1] = Lambda[0][0]
Lambda[0][0] = temp
for j in 0,1,2:
  temp = Lambda[1][j,1]
  Lambda[1][j,1] = Lambda[1][j,0]
  Lambda[1][j,0] = temp

temp = Lambda[0][2]
Lambda[0][2] = Lambda[0][1]
Lambda[0][1] = temp
for j in 0,1,2:
  temp = Lambda[1][j,2]
  Lambda[1][j,2] = Lambda[1][j,1]
  Lambda[1][j,1] = temp

print(Lambda)
ql = [rho_L, rho_L * v_L, rho_L*v_L**2/2 + p_L/(gamm-1)]
qr = [rho_R, rho_R * v_R, rho_R*v_R**2/2 + p_R/(gamm-1)]
Vl = np.linalg.solve(Lambda[1], ql)
Vr = np.linalg.solve(Lambda[1], qr)
ql_dot = Vr[0]*Lambda[1][:,0] + Vl[1]*Lambda[1][:,1] + Vl[2]*Lambda[1][:,2]
#ql_dot = Vr[0]*Lambda[1][0,:] + Vl[1]*Lambda[1][1,:] + Vl[2]*Lambda[1][2,:]
qr_dot = Vr[0]*Lambda[1][:,0] + Vr[1]*Lambda[1][:,1] + Vl[2]*Lambda[1][:,2]
print(ql_dot)
print(qr_dot)
rho = [ql[0], ql_dot[0], qr_dot[0], qr[0]]
v_x = [ql[1]/ql[0], ql_dot[1]/ql_dot[0], qr_dot[1]/qr_dot[0], qr[1]/qr[0]]
p = [(ql[2]-ql[1]**2/(ql[0]*2)) * (gamm-1.0), 
     (ql_dot[2]-ql_dot[1]**2/(ql_dot[0]*2)) * (gamm-1.0),
     (qr_dot[2]-qr_dot[1]**2/(qr_dot[0]*2)) * (gamm-1.0),
     (qr[2]-qr[1]**2/(qr[0]*2)) * (gamm-1.0)]
print('rho = ',rho)
print('v_x = ',v_x)
print('p = ',p)