# Line to Line fault analysis...
# Code written by Knapoj Chaimanekorn based on the building algorithm 
# presented in Hadi Saadat's book.
# 14 August 2022

import numpy as np
from tabulate import tabulate
from .faultutils import *

def llfault(zdata0,zbus0,zdata1,zbus1,zdata2,zbus2, vbus = None):
    if type(vbus) != np.ndarray: # Use unity prefault voltage if not specified.
        vbus = np.ones(len(zbus1))
    else: pass

    fbus, zf, k = faultinput(zdata1)
    i1f = (vbus[k]/(zf + zbus1[k,k] + zbus2[k,k]))
    ik012f = np.array([0, i1f,-i1f])
    ikabcf = np.matmul(matA(),ik012f)

    print('\nLine to-line fault at bus No. {}'.format(fbus))
    print('Total fault current = {:.4f} Per unit'.format(np.abs(np.max(ikabcf))))

    vf012, vtbl = calculatevoltages(vbus,zbus0, zbus1, zbus2, ik012f, k)
    print('\nBus Voltages during the fault in per unit')
    print(vtbl)

    itbl = calculatecurrents(vf012,ikabcf, fbus, zdata0, zdata1, zdata2)
    print('\nLine currents for fault at bus No. {}'.format(fbus))
    print(itbl)