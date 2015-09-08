# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 10:41:02 2014

@author: root
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 14:59:21 2013

@author: root
"""

import os
import csv
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import matplotlib
from scipy.optimize import leastsq
import scipy.constants as c
import random
from decimal import *
#import sympy
#from sympy.matrices import Matrix
import scipy.sparse as sparse

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
from matplotlib.colors import colorConverter

import matplotlib
matplotlib.rcParams.update({'font.size': 20, 'font.family': 'serif'})


#import cProfile, pstats, StringIO
#pr = cProfile.Profile()
#pr.enable()


class NElectronMatrix:
    
    def __init__(self,N,bi,Qo,omega0,omegaVector,QVector,Nprec):
        self.params = {'N': N, 'bi': bi, 'Qo':Qo,'omega0':omega0,'omegaVector':omegaVector,'QVector':QVector,'Nprec':Nprec}
        self.N = N
        self.bi =bi
        self.Qo = Qo
        self.omega0 = omega0
        self.omegaVector= omegaVector
        self.Qi= Qi
        self.Nprec = Nprec
       # self.matrix = self.MakeMatrix(N,bi,bo,Qo,omegaVector)
       
    def reset(self,newparams):
        for k in newparams.keys():
            self.params[k] = newparams[k]
        
        self.N = self.params['N']
        self.bi =self.params['bi']
        self.Qo = self.params['Qo']
        self.omega0 = self.params['omega0']
        self.omegaVector= self.params['omegaVector']
        self.Qi= self.params['Qi']
        self.Nprec = self.params['Nprec']
        
        
            
        
    def MakeMatrix(self,N,bi,bo,Qo,Qi,omegaVector):
    
        Zero = 0.0#Decimal(0)
    
        A11 = np.zeros([N+1,N+1])
        #A11 = np.array([[Decimal(x) for x in y ] for y in A11])
    
        A12 = np.identity(N+1)
        #A12 = np.array([[Decimal(x) for x in y ] for y in A12])

        B11 = np.zeros([N+1,N+1])
        B12 = np.zeros([N+1,N+1])
        
        k=bi

       
        ###### make B11

        for i in range(0,N+1):
        
            if i==0:
                row = np.append(-omega0**2 /(1-N*k**2),bo*omegaVector[0][j]**2/(1-N*k**2))
    
            else:
                row = np.append(bi*omega0**2/(1-N*k**2),-k**2*omegaVector[0][j]**2/(1-N*k**2))
                row[i]=-omegaVector[0][j]**2*(1 + ( k**2/(1-N*k**2) ) )
        
            B11.append(row)
    

        ###################################### make B12
        for i in range(0,N+1):
            row = []
    
            if i==0:
                row=np.append(-omega0/Qo/(1-N*k**2),bo*omegaVector[0][j]/Qi/(1-N*k**2))  
            else:
                row = np.append(bi*omega0/Qo/(1-N*k**2),-k**2*omegaVector[0][j]/Qi/(1-N*k**2))
                row[i]=-omegaVector[0][j]/Qi*(1 + ( k**2/(1-N*k**2) ) )       
        
            B12.append(row)

        M = np.vstack([np.hstack([A11,A12]),np.hstack([B11,B12])])
        
        #print "Matrix Sparse? ",sparse.issparse(M)

        return M
        
    def getMeanRates(self):
        bo = bi = self.bi
        Qo = self.Qo
        omega0= self.omega0
        omegaVector = self.omegaVector
        Nprec = self.Nprec
        Qi=self.Qi
        N =self.N
        
        bo= bi
        #k=np.sqrt(bo*bi)
        MM = self.MakeMatrix(N,bi,bo,Qo,Qi,omegaVector)
        w=LA.eigvals(MM)       
        vvN = np.sort(np.real(w))#[i.real for i in w ]
        
        self.eigenvalues = w

        self.cavRate = vvN[:2]
        self.eRate = vvN[2:]
        
        
        return [self.cavRate,self.eRate]        

        
def binFreqs(x,xbins):
    Ba=[]

    for i in range(0,len(xbins)-1):
        inds = np.where((x<xbins[i+1])*(x>xbins[i]))
        Ba.append(len(x[inds]))
        #db.append(np.std(B[inds]))
    
    Ba=np.array([Ba])
    x = x[1:]
      
    return Ba
    
    
Nprec = 16
bi=3.08e-6 
bi=5e-5
bo= bi

    
Qo = 1600 # TE121 simulated Q, nichrome, free space ends
Qi = 1.2e11 # free space Q for e at 30 GHz, 4 s decay time
omega0 = 2.2e10*2*np.pi

ModeBW = omega0/Qo
MGamma = ModeBW/2/np.pi
repeats = 1
kScan = 10.0**np.arange(-5,-4,0.1)

omegacScan = np.arange(-2,2,0.1)*ModeBW +omega0
sigma = ModeBW/100
omega_spec = omegacScan#np.arange(-2,2,0.5)*ModeBW +omega0

eBW = omega0/Qi
gBins= np.arange(0,5*ModeBW,ModeBW/100)/2/np.pi


Numbers = np.arange(500,501,1)

NeRate = np.zeros(len(Numbers))
NCRate =  np.zeros(len(Numbers))


n=Numbers[0]

print n," electrons"
Dm = np.zeros((len(omegacScan),2*(n+1),repeats))
Dm1 = np.zeros((len(omegacScan),len(omega_spec)-1))
Dm2 = np.zeros((len(omegacScan),len(gBins)-1))


PsF = np.zeros((repeats,2*(n+1)))
PsF2 = np.zeros((repeats,2*(n+1)))



for j in range(0,len(omegacScan)):
    print "solveing omega ",j," out of ",len(omegacScan)
    for m in range(0,repeats):
        #print "simga" ,sig
        omegaVector = [np.append([omega0],sigma*np.random.randn(n)+omegacScan[j])]
#        plt.hist((omegaVector[0][1:]-omega0)/ModeBW)
#        plt.show()
        Simulation = NElectronMatrix(n,bi,Qo,omega0,omegaVector,Qi,Nprec)
        CavRate,eRate = Simulation.getMeanRates()
        eigenW = np.abs(np.imag(Simulation.eigenvalues))
        eigenG = np.abs(np.real(Simulation.eigenvalues))
        #Dm[:,m] = eigenW
        PsF[m] =eigenW
        PsF2[m] =eigenG
        
        ws = np.ones(len(eigenW))*omegacScan[j]
        plt.subplot(211)
        plt.plot((ws-omega0)/ModeBW,(eigenW-omega0)/ModeBW,'o')
        plt.subplot(212)
        plt.plot((ws-omega0)/ModeBW,np.log10(eigenG/eBW),'o')
        
            
        F = PsF.reshape(-1)
        F2 = PsF2.reshape(-1)
    Dm1[j,:] = binFreqs(F,omega_spec)
    Dm2[j,:] = binFreqs(F2,gBins)
        #print np.sum(Dm1[i:,]) ," ==? ",repeats*2*(n+1)
        
        
plt.subplot(211)
plt.ylabel('(f - fo)/$\Delta$ fo (MHz)')
plt.subplot(212)
plt.ylabel('log10( $\Gamma / \Gamma_{freespace}$)')
plt.xlabel('(f_cyc - fo)/$\Delta$ fo (MHz)')
plt.show()            




matplotlib.rcParams.update({'font.size': 20, 'font.family': 'serif'})

detuning = (omegacScan[1:]-omega0)/2/np.pi
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#for c, z in zip(['r', 'g', 'b', 'y'], detuning):
#    xs =detuning#(omegacScan[1:]-omega0)/2/np.pi# np.arange(20)
#    ys = Dm1[i,:]#np.random.rand(20)
#    cs = [c] * len(xs)
#    cs[0] = 'c'
#    ax.bar(xs, ys, zs=z, zdir='y', color=cs, alpha=0.8)
#
#ax.set_xlabel('Detuning')
#ax.set_ylabel('k')
#ax.set_zlabel('PSF')
#
#plt.show()

#fig = plt.figure()
#
#ax = fig.gca(projection='3d')
#
#
#cc = lambda arg: colorConverter.to_rgba(arg, alpha=0.6)
#xs = detuning #np.arange(0, 10, 0.4)
#verts = []
#zs = detuning #[0.0, 1.0, 2.0, 3.0,4.0,5.0,6.0,7.0]
#for i in range(0,len(zs)):
#    ys = Dm1[i,:]#Gamma[i]#np.random.rand(len(xs))
#    xs = detuning#omegacScan[1:]#24.17/250.05*BN[i]
#    ys[0], ys[-1] = 0, 0
#    verts.append(list(zip(xs, ys)))
#    #plt.plot(xs,ys)
#    #plt.show()
#
##cs = [c] * len(xs)
#poly = PolyCollection(verts, facecolors = [cc('w')])#, cc('g'), cc('b'), cc('y'),cc('r'), cc('g'), cc('b'), cc('y')])
#poly.set_alpha(0.7)
#ax.add_collection3d(poly, zs=zs, zdir='y')
#
#ax.set_xlabel('Eigen frequency  (GHz)')
#ax.set_xlim3d(-100e6,100e6)
#ax.set_ylabel('cyclotron frequency')
#ax.set_ylim3d(-100e6,100e6)
#ax.set_zlabel('PSF')
#ax.set_zlim3d(0, 100)
#plt.show()
#
#fig = plt.figure()
#
#ax = fig.gca(projection='3d')
#cc = lambda arg: colorConverter.to_rgba(arg, alpha=0.6)
#xs = detuning #np.arange(0, 10, 0.4)
#verts = []
#zs = detuning #[0.0, 1.0, 2.0, 3.0,4.0,5.0,6.0,7.0]
#for i in range(0,len(zs)):
#    ys = Dm2[i,:]#Gamma[i]#np.random.rand(len(xs))
#    xs = gBins/MGamma#omegacScan[1:]#24.17/250.05*BN[i]
#    ys[0], ys[-1] = 0, 0
#    verts.append(list(zip(xs, ys)))
#    #plt.plot(xs,ys)
#    #plt.show()
#
##cs = [c] * len(xs)
#poly = PolyCollection(verts, facecolors = [cc('w')])#, cc('g'), cc('b'), cc('y'),cc('r'), cc('g'), cc('b'), cc('y')])
#poly.set_alpha(0.7)
#ax.add_collection3d(poly, zs=zs, zdir='y')
#
#ax.set_xlabel('Gamma  (1/s)')
#ax.set_xlim3d(0,5)
#ax.set_ylabel('cyclotron frequency')
#ax.set_ylim3d(-100e6,100e6)
#ax.set_zlabel('PSF')
#ax.set_zlim3d(0, 100)
#
#plt.show()