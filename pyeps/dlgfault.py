# Double Line to Ground fault analysis...
# Code written by Knapoj Chaimanekorn based on the building algorithm 
# presented in Hadi Saadat's book.
# 14 August 2022

import numpy as np
from tabulate import tabulate
from .faultutils import *

def dlgfault(zdata0,zbus0,zdata1,zbus1,zdata2,zbus2, vbus = None):
    if type(vbus) != np.ndarray: # Use unity prefault voltage if not specified.
        vbus = np.ones(len(zbus1))
    else: pass

    fbus, zf, k = faultinput(zdata1)
    ik1 = vbus[k]/(zbus1[k,k] + zbus2[k,k]*(zbus0[k,k] + 3*zf)/(zbus2[k,k] + zbus0[k,k] + 3*zf))
    ik2 = - (vbus[k] - zbus1[k,k]*ik1)/zbus2[k,k]
    ik0 = - (vbus[k] - zbus1[k,k]*ik1)/(zbus0[k,k] + 3*zf)
    ik012f = np.array([ik0,ik1,ik2])
    ikabcf = np.matmul(matA(),ik012f)

    print('\nDouble line to-ground fault at bus No. {}'.format(fbus))
    print('Total fault current = {:.4f} Per unit'.format(np.abs(np.sum(ikabcf))))

    vf012, vtbl = calculatevoltages(vbus,zbus0, zbus1, zbus2, ik012f, k)
    print('\nBus Voltages during the fault in per unit')
    print(vtbl)

    itbl = calculatecurrents(vf012,ikabcf, fbus, zdata0, zdata1, zdata2)
    print('\nLine currents for fault at bus No. {}'.format(fbus))
    print(itbl)