# Balanced fault analysis...
# Code written by Knapoj Chaimanekorn based on the building algorithm 
# presented in Hadi Saadat's book.
# 14 August 2022

import numpy as np
from tabulate import tabulate
from .faultutils import *

def symfault(zdata,zbus,vbus = None):

    if type(vbus) != np.ndarray: # Use unity prefault voltage if not specified.
        vbus = np.ones(len(zbus))
    else: pass
    
    fbus, zf, k = faultinput(zdata)
    
    ikf = vbus[k]/(zf + zbus[k,k])
    print('\nBalanced three-phase fault at bus No. {}'.format(fbus))
    print('Total fault current = {:.4f} Per unit'.format(np.abs(ikf)))

    print('\nBus Voltages during the fault in per unit')
    delv = zbus[:,k]*vbus[k]/(zbus[k,k] + zf)
    vf = vbus - delv

    vdata = np.array([range(1,1+len(vbus)), np.abs(vf), np.rad2deg(np.angle(vf))]).T
    print(tabulate(np.round(vdata,4), headers = ['Bus No.', 'Magnitude (pu.)', 'Angle (deg.)']))

    # Calculate line current using zdata...
    
    zdata = zdata[zdata[:,1].argsort()]
    zdata = zdata[zdata[:,0].argsort(kind = 'mergesort')]

    iline = np.array([], dtype = complex)
    frombus = np.array([], dtype = int)
    tobus = np.array([], dtype = int)
    for ide, e in enumerate(zdata):
        p,q = int(e[0]) - 1,  int(e[1]) - 1
        zpq = complex(e[2], e[3])
        if p < 0: continue
        
        if np.angle((vf[p] - vf[q])/zpq) > 0:
            iline = np.append(iline, (vf[q] - vf[p])/zpq)
            frombus, tobus = np.append(frombus, q +1), np.append(tobus, p +1)
        else: 
            iline = np.append(iline, (vf[p] - vf[q])/zpq)
            frombus, tobus = np.append(frombus, p +1), np.append(tobus, q +1)

    idata = np.array([frombus, tobus, np.abs(iline), np.rad2deg(np.angle(iline))]).T
    
    idata = np.append(idata, [[fbus , -1,np.abs(ikf), np.rad2deg(np.angle(ikf))]], axis = 0)
    idata = idata[idata[:,1].argsort()]
    idata = idata[idata[:,0].argsort(kind = 'mergesort')] # Sort data

    itbl = np.char.replace(
        tabulate(np.round(idata,4), headers = ['From', 'To', 'Magnitude (pu.)', 'Angle (deg.)'], numalign = 'right'), '-1 ', ' F ')

    print('\nLine currents for fault at bus No. {}'.format(fbus))
    print(itbl)