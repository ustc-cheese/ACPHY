# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 21:23:28 2020

@author: HP
"""
import numpy as np
import Matrix as M
import os

global m0,e0,c0

m0 = 9.10938215e-31
e0 = 1.602176487e-19
c0 = 299792458

'''
beam parameters
'''
gamma0 = 3053
emitn = 0.6e-6
emitx = 0.6e-6
emity = 0.6e-6
Kf = 16
Kd = 13.5
'''
Chicane parameters
'''
r56_1 =  0.00260627
r56_2 = 0.0000920597
imagl = 0.5
ibfield = 0.5
def idril_cal(r56):
    theta = np.arcsin(imagl/(np.sqrt(gamma0**2 - 1) * m0 * c0) * ibfield * e0)**2
    return np.abs(2/3 * imagl-r56/theta/2)
idril1 =  idril_cal(r56_1)
idril2 =  idril_cal(r56_2)

'''
modulator parameters
'''
lambda_seed = 270e-9
lambdam = 0.08
periodm1 = 20
periodm2 = 10
aum = np.sqrt(lambda_seed*2*gamma0**2/lambdam-1)
km = np.pi*2/lambdam
Km = np.sqrt(2)*aum
Kbetam = 1/2*(Km*km/gamma0)**2

'''
radiator parameters
'''
xlamds = 270e-9/30
lambdau = 0.04
periodu = 74
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
lm1 = periodm1*lambdam
lm2 = periodm2*lambdam
l2 = 0
lqf = 4*lambdau
l3 = 10*lambdau
lu = periodu*lambdau
l4 = 10*lambdau
lqd = 4*lambdau
l5 = 10*lambdau

'''
Horizonal Matrix
'''
Ux = M.TransQuadF(kf,lqf/2).dot(M.TransDrift(l3).dot(M.TransUnduH(lu).dot(M.TransDrift(l4) \
    .dot(M.TransQuadD(kd,lqd).dot(M.TransDrift(l5).dot(M.TransUnduH(lu).dot(M.TransDrift(l4) \
    .dot(M.TransQuadF(kf,lqf/2)))))))))
if (Ux[0,0] + Ux[1,1])**2 > 4 or Ux[0,1] < 0:
    os._exit()
    
Uy = M.TransQuadD(kf,lqf/2).dot(M.TransDrift(l3).dot(M.TransUnduV(Kbetau,lu).dot(M.TransDrift(l4) \
    .dot(M.TransQuadF(kd,lqd).dot(M.TransDrift(l5).dot(M.TransUnduV(Kbetau,lu).dot(M.TransDrift(l4) \
    .dot(M.TransQuadD(kf,lqf/2)))))))))    
if (Uy[0,0] + Uy[1,1])**2 > 4 or Uy[0,1] < 0:
    os._exit()


    
'''
twiss parameter
'''
alphax = (Ux[0,0]-Ux[1,1])/np.sqrt(4-(Ux[0,0]+Ux[1,1])**2)
betax = 2*Ux[0,1]/np.sqrt(4-(Ux[0,0]+Ux[1,1])**2)
gammax = (1 + alphax**2)/betax

alphay = (Uy[0,0]-Uy[1,1])/np.sqrt(4-(Uy[0,0]+Uy[1,1])**2)
betay = 2*Uy[0,1]/np.sqrt(4-(Uy[0,0]+Uy[1,1])**2)
gammay = (1 + alphay**2)/betay

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
'''
entance of modulator
'''
Mx = M.TransChica(imagl,idril2,ibfield,gamma0,1).dot(M.TransUnduH(lm2))
My = M.TransChica(imagl,idril2,ibfield,gamma0,0).dot(M.TransUnduV(Kbetam,lm2))

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

Mx = M.TransChica(imagl,idril1,ibfield,gamma0,1).dot(M.TransUnduH(lm1))
My = M.TransChica(imagl,idril1,ibfield,gamma0,0).dot(M.TransUnduV(Kbetam,lm1))

Ax = np.linalg.inv(Mx)
Ay = np.linalg.inv(My)

Nx = np.array([[Ax[0,0]**2,-2*Ax[0,0]*Ax[0,1],Ax[0,1]**2], \
              [-Ax[0,0]*Ax[1,0],1+2*Ax[0,1]*Ax[1,0],-Ax[0,1]*Ax[1,1]], \
              [Ax[1,0]**2,-2*Ax[1,0]*Ax[1,1],Ax[1,1]**2]])
Ny = np.array([[Ay[0,0]**2,-2*Ay[0,0]*Ay[0,1],Ay[0,1]**2], \
              [-Ay[0,0]*Ay[1,0],1+2*Ay[0,1]*Ay[1,0],-Ay[0,1]*Ay[1,1]], \
              [Ay[1,0]**2,-2*Ay[1,0]*Ay[1,1],Ay[1,1]**2]])

Bx = Nx.dot(np.array([betax1,alphax1,gammax1]))    
By = Ny.dot(np.array([betay1,alphay1,gammay1]))    

alphax2 = Bx[1]
alphay2 = By[1]
betax2 = Bx[0]
betay2 = By[0]
gammax2 = Bx[2]
gammay2 = By[2]
sigmax2 = np.sqrt(betax2*emitn/gamma0)
sigmay2 = np.sqrt(betay2*emitn/gamma0)
print(alphax2,alphay2,sigmax2,sigmay2)