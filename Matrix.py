# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 20:00:58 2020

@author: HP
"""
import numpy as np
import matplotlib.pyplot as mpl

global m0,e0,c0
m0 = 9.10938215e-31
e0 = 1.602176487e-19
c0 = 299792458
    
def TransDrift(s):
    return np.array([[1,s],[0,1]])
    
def TransEdgeX(theta,rho):
    return np.array([[1,0],[np.tan(theta/rho),1]])
    
def TransEdgeY(theta,rho):
    return np.array([[1,0],[-np.tan(theta/rho),1]])
    
def TransQuadF(k,s):
    a = np.cos(np.sqrt(k)*s)
    b = np.sin(np.sqrt(k)*s)/np.sqrt(k)
    c = -np.sqrt(k)*np.sin(np.sqrt(k)*s)
    d = np.cos(np.sqrt(k)*s)
    return np.array([[a,b],[c,d]])
    
def TransQuadD(k,s):
    a = np.cosh(np.sqrt(k)*s)
    b = np.sinh(np.sqrt(k)*s)/np.sqrt(k)
    c = np.sqrt(k)*np.sinh(np.sqrt(k)*s)
    d = np.cosh(np.sqrt(k)*s)
    return np.array([[a,b],[c,d]])
    
def TransSectX(theta,rho):
    return np.array([[np.cos(theta),rho*np.sin(theta)],[-np.sin(theta)/rho,np.cos(theta)]])

def TransSectY(theta, rho):
    return np.array([[1,rho*theta],[0,1]])
    
def TransChica(imagl,idril,ibfield,gamma0,xoy):
    rho = np.sqrt(gamma0**2 - 1)*m0*c0/ibfield/e0
    theta = np.arcsin(imagl/rho)
    if xoy == 1:
        return TransDrift(idril)*TransSectX(theta, rho)*TransEdgeX(theta, rho) \
        *TransDrift(idril)*TransEdgeX(-theta, -rho)*TransSectX(-theta, -rho) \
        *TransDrift(idril)*TransSectX(-theta, -rho)*TransEdgeX(-theta, -rho) \
        *TransDrift(idril)*TransEdgeX(theta,rho)*TransSectX(theta,rho)*TransDrift(idril)
    elif xoy == 0:
        return TransDrift(idril)*TransSectY(theta,rho)*TransEdgeY(theta,rho) \
        *TransDrift(idril)*TransEdgeY(-theta,-rho)*TransSectY(-theta,-rho) \
        *TransDrift(idril)*TransSectY(-theta,-rho)*TransEdgeY(-theta,-rho) \
        *TransDrift(idril)*TransEdgeY(theta,rho)*TransSectY(theta,rho)*TransDrift(idril)

def TransUnduH(s):
    return np.array([[1,s],[0,1]])

def TransUnduV(k,s):
    return TransQuadF(k,s)
         
