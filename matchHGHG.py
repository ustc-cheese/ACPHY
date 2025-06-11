# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 21:23:28 2020

@author: HP
"""
import numpy as np
import Matrix as M

global m0,e0,c0

m0 = 9.10938215e-31
e0 = 1.602176487e-19
c0 = 299792458

'''
beam parameters
'''
gamma0 = 800/0.511
emitn = 1e-6
emitx =1e-6
emity = 1e-6
Kf = 10
Kd = 10
'''
Chicane parameters
'''
r56 = 0.00012
imagl = 0.1
ibfield = 0.5
def idril_cal(r56):
    theta = np.arcsin(imagl/(np.sqrt(gamma0**2 - 1) * m0 * c0) * ibfield * e0)**2
    return np.abs(2/3 * imagl-r56/theta/2)
idril =  idril_cal(r56)

'''
modulator parameters
'''
lambda_seed = 266e-9
lambdam = 0.08
periodm = 20
aum = np.sqrt(lambda_seed*2*gamma0**2/lambdam-1)
km = np.pi*2/lambdam
Km = np.sqrt(2)*aum
Kbetam = 1/2*(Km*km/gamma0)**2

'''
radiator parameters
'''
xlamds = 13.5e-9
lambdau = 0.03
periodu = 100
aur = np.sqrt(xlamds*2*gamma0**2/lambdau-1)
ku = np.pi*2/lambdau
Ku = np.sqrt(2)*aur
Kbetau = 1/2*(Ku*ku/gamma0)**2

kf = e0*Kf/m0/c0/np.sqrt(gamma0**2-1)
kd = e0*Kd/m0/c0/np.sqrt(gamma0**2-1)
'''
beamline setup
'''
l1 = 0
lm = periodm*lambdam
l2 = 0
lqf = 4*lambdau
l3 = 16*lambdau
lu = periodu*lambdau
l4 = 17*lambdau
lqd = 4*lambdau
l5 = 16*lambdau
lu = periodu*lambdau
l6 = 17*lambdau

'''
initial Twiss
'''
alphax_ini =-4.259797e-09
alphay_ini = 2.224942e-03 
betax_ini = 7.9509e-05**2/emitn*gamma0
betay_ini = 5.2559e-05**2/emitn*gamma0
gammax_ini = (1+alphax_ini**2)/betax_ini
gammay_ini = (1+alphay_ini**2)/betay_ini
sigmax_ini = np.sqrt(betax_ini*emitn/gamma0)
sigmay_ini = np.sqrt(betax_ini*emitn/gamma0)

'''
Horizonal Matrix
'''
Ux = ((((((((M.TransDrift(l6) @ M.TransUnduH(lu)) @ M.TransDrift(l5)) @ M.TransQuadD(kd, lqd)) 
      @ M.TransDrift(l4)) @ M.TransUnduH(lu)) @ M.TransDrift(l3)) @ M.TransQuadF(kf, lqf)) @ M.TransDrift(l2))
if (Ux[0,0] + Ux[1,1])**2 > 4 or Ux[0,1] < 0:
    print('x-axis no sulution')
    
Uy = ((((((((M.TransDrift(l6) @ M.TransUnduV(Kbetau, lu)) @ M.TransDrift(l5)) @ M.TransQuadF(kf, lqd)) 
      @ M.TransDrift(l4)) @ M.TransUnduV(Kbetau, lu)) @ M.TransDrift(l3)) @ M.TransQuadD(kd, lqf)) @ M.TransDrift(l2))
if (Uy[0,0] + Uy[1,1])**2 > 4 or Uy[0,1] < 0:
    print('y-axis no sulution')

'''
twiss parameter
'''

Ux = np.linalg.inv(Ux)
Uy = np.linalg.inv(Uy)

UMx = np.array([[Ux[0,0]**2,-2*Ux[0,0]*Ux[0,1],Ux[0,1]**2], \
              [-Ux[0,0]*Ux[1,0],1+2*Ux[0,1]*Ux[1,0],-Ux[0,1]*Ux[1,1]], \
              [Ux[1,0]**2,-2*Ux[1,0]*Ux[1,1],Ux[1,1]**2]])
UMy = np.array([[Uy[0,0]**2,-2*Uy[0,0]*Uy[0,1],Uy[0,1]**2], \
              [-Uy[0,0]*Uy[1,0],1+2*Uy[0,1]*Uy[1,0],-Uy[0,1]*Uy[1,1]], \
              [Uy[1,0]**2,-2*Uy[1,0]*Uy[1,1],Uy[1,1]**2]])

Bx = UMx.dot(np.array([betax_ini,alphax_ini,gammax_ini]))    
By = UMy.dot(np.array([betay_ini,alphay_ini,gammay_ini]))
  

alphax = Bx[1]
alphay = By[1]
betax = Bx[0]
betay = By[0]
gammax = Bx[2]
gammay = By[2]
sigmax = np.sqrt(betax*emitn/gamma0)
sigmay = np.sqrt(betay*emitn/gamma0)

print(alphax,alphay,sigmax,sigmay)
'''
entance of radiator
'''
Ax = M.TransQuadF(kf,lqf/2).dot(M.TransDrift(l3))
Ay = M.TransQuadD(kf,lqf/2).dot(M.TransDrift(l3))

Ax = np.linalg.inv(Ax)
Ay = np.linalg.inv(Ay)

Nx = np.array([[Ax[0,0]**2,-2*Ax[0,0]*Ax[0,1],Ax[0,1]**2], \
              [-Ax[0,0]*Ax[1,0],1+2*Ax[0,1]*Ax[1,0],-Ax[0,1]*Ax[1,1]], \
              [Ax[1,0]**2,-2*Ax[1,0]*Ax[1,1],Ax[1,1]**2]])
Ny = np.array([[Ay[0,0]**2,-2*Ay[0,0]*Ay[0,1],Ay[0,1]**2], \
              [-Ay[0,0]*Ay[1,0],1+2*Ay[0,1]*Ay[1,0],-Ay[0,1]*Ay[1,1]], \
              [Ay[1,0]**2,-2*Ay[1,0]*Ay[1,1],Ay[1,1]**2]])

Bx = Nx.dot(np.array([betax,alphax,gammax]))    
By = Ny.dot(np.array([betay,alphay,gammay]))    

alphax0 = Bx[1]
alphay0 = By[1]
betax0 = Bx[0]
betay0 = By[0]
gammax0 = Bx[2]
gammay0 = By[2]
sigmax0 = np.sqrt(betax0*emitn/gamma0)
sigmay0 = np.sqrt(betay0*emitn/gamma0)

print(alphax0,alphay0,sigmax0,sigmay0)

Mx = M.TransChica(imagl,idril,ibfield,gamma0,1).dot(M.TransUnduH(lm))
My = M.TransChica(imagl,idril,ibfield,gamma0,0).dot(M.TransUnduV(Kbetam,lm))

Ax = np.linalg.inv(Mx)
Ay = np.linalg.inv(My)

Nx = np.array([[Ax[0,0]**2,-2*Ax[0,0]*Ax[0,1],Ax[0,1]**2], \
              [-Ax[0,0]*Ax[1,0],1+2*Ax[0,1]*Ax[1,0],-Ax[0,1]*Ax[1,1]], \
              [Ax[1,0]**2,-2*Ax[1,0]*Ax[1,1],Ax[1,1]**2]])
Ny = np.array([[Ay[0,0]**2,-2*Ay[0,0]*Ay[0,1],Ay[0,1]**2], \
              [-Ay[0,0]*Ay[1,0],1+2*Ay[0,1]*Ay[1,0],-Ay[0,1]*Ay[1,1]], \
              [Ay[1,0]**2,-2*Ay[1,0]*Ay[1,1],Ay[1,1]**2]])

Bx = Nx.dot(np.array([betax0,alphax0,gammax0]))    
By = Ny.dot(np.array([betay0,alphay0,gammay0]))    

alphax1 = Bx[1]
alphay1 = By[1]
betax1 = Bx[0]
betay1 = By[0]
gammax1 = Bx[2]
gammay1 = By[2]
sigmax1 = np.sqrt(betax1*emitn/gamma0)
sigmay1 = np.sqrt(betay1*emitn/gamma0)
print(alphax1,alphay1,sigmax1,sigmay1)