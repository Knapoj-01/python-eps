# Fault analysis utility functions...
# Code written by Knapoj Chaimanekorn based on the building algorithm 
# presented in Hadi Saadat's book.
# 14 August 2022

import numpy as np
from tabulate import tabulate

def faultinput(zdata):
    fbus = int(input("Enter Faulted Bus No. -> "))
    m = int(np.max(zdata[:,:2])) # Check 
    if fbus > m:
        raise TypeError("Error: No bus found.")
    zf = input("Enter Fault Impedance Zf = R + Xj in complex form (for bolted fault enter 0). Zf = ")
    zf = complex(zf.replace(" ", "")) # Casting
    k = fbus - 1 # Index
    return fbus, zf, k

def calculatecurrents(vf012, ikabcf, fbus, zdata0, zdata1, zdata2):
    
    # Calculate line current using zdata...
    iline0 = np.array([], dtype = complex)
    iline1 = np.array([], dtype = complex)
    iline2 = np.array([], dtype = complex)
    frombus = np.array([], dtype = int)
    tobus = np.array([], dtype = int)

    # Sort data
    zdata0 = zdata0[zdata0[:,1].argsort()]
    zdata0 = zdata0[zdata0[:,0].argsort(kind = 'mergesort')]
    zdata1 = zdata1[zdata1[:,1].argsort()]
    zdata1 = zdata1[zdata1[:,0].argsort(kind = 'mergesort')]
    zdata2 = zdata2[zdata2[:,1].argsort()]
    zdata2 = zdata2[zdata2[:,0].argsort(kind = 'mergesort')]

    for i in range(len(zdata0)):
        p,q = int(zdata0[i,0]) - 1,  int(zdata0[i,1]) - 1
        zpq0 = complex(zdata0[i,2], zdata0[i,3])
        zpq1 = complex(zdata1[i,2], zdata1[i,3])
        zpq2 = complex(zdata2[i,2], zdata2[i,3])

        if p < 0: continue
        if np.angle((vf012[p,0] - vf012[q,0])/zpq0) > 0:
            iline0 = np.append(iline0, (vf012[q,0] - vf012[p,0])/zpq0)
            iline1 = np.append(iline1, (vf012[q,1] - vf012[p,1])/zpq1)
            iline2 = np.append(iline2, (vf012[q,2] - vf012[p,2])/zpq2)
            frombus, tobus = np.append(frombus, q +1), np.append(tobus, p +1)
        else:
            iline0 = np.append(iline0, (vf012[p,0] - vf012[q,0])/zpq0)
            iline1 = np.append(iline1, (vf012[p,1] - vf012[q,1])/zpq1)
            iline2 = np.append(iline2, (vf012[p,2] - vf012[q,2])/zpq2)
            frombus, tobus = np.append(frombus, p +1), np.append(tobus, q +1)

    # Convert to abc.
    iline012 = np.array([iline0, iline1, iline2]).T
    ilineabc = np.copy(iline012)
    for i in range(len(iline012)):
        ilineabc[i] = np.matmul(matA(), iline012[i])

    idata = np.array([
        frombus, tobus, np.abs(ilineabc[:,0]), np.rad2deg(np.angle(ilineabc[:,0])),
        np.abs(ilineabc[:,1]), np.rad2deg(np.angle(ilineabc[:,1])),
        np.abs(ilineabc[:,2]), np.rad2deg(np.angle(ilineabc[:,2])),
    ]).T

    idata = np.append(idata, [[
        fbus , -1,np.abs(ikabcf[0]), np.rad2deg(np.angle(ikabcf[0])),
        np.abs(ikabcf[1]), np.rad2deg(np.angle(ikabcf[1])),
        np.abs(ikabcf[2]), np.rad2deg(np.angle(ikabcf[2]))
    ]], axis = 0)

    idata = idata[idata[:,1].argsort()]
    idata = idata[idata[:,0].argsort(kind = 'mergesort')] # Sort data

    itbl = np.char.replace(
        tabulate(np.round(idata,4), headers = [
            'From', 'To', 'Phase A\nMag. (pu.)', 'Phase A\nAng. (deg.)',
            'Phase B\nMag. (pu.)', 'Phase B\nAng. (deg.)',
            'Phase C\nMag. (pu.)', 'Phase C\nAng. (deg.)'
        ], numalign = 'right'), '-1 ', ' F ')
    return(itbl)

def calculatevoltages(vbus,zbus0, zbus1, zbus2, ik012f, k):
    v012 = np.zeros((3,len(vbus)), dtype = complex).T
    v012[:,1] = vbus
    delv012 = np.copy(v012)
    
    for i in range(len(delv012)):
        delv012[i] = np.multiply(np.array([zbus0[i, k], zbus1[i, k], zbus2[i, k]]), ik012f)
    
    vf012 = v012 - delv012
    vfabc = np.copy(vf012)
    for i in range(len(vf012)):
        vfabc[i] = np.matmul(matA(), vf012[i])

    vdata = np.array([
        range(1,1+len(vbus)), np.abs(vfabc[:,0]), np.rad2deg(np.angle(vfabc[:,0])),
        np.abs(vfabc[:,1]), np.rad2deg(np.angle(vfabc[:,1])),
        np.abs(vfabc[:,2]), np.rad2deg(np.angle(vfabc[:,2]))
    ]).T
    vtbl = tabulate(np.round(vdata,4), 
                   headers = [
                       'Bus No.', 'Phase A\nMag. (pu.)', 'Phase A\nAng. (deg.)',
                       'Phase B\nMag. (pu.)', 'Phase B\nAng. (deg.)',
                       'Phase C\nMag. (pu.)', 'Phase C\nAng. (deg.)'
                   ])
    return vf012, vtbl

def matA():
    a = np.exp(2j*np.pi/3)
    return np.array([[1,1,1],[1,a**2,a],[1,a,a**2]])

def matAinv():
    A = matA()
    return np.linalg.inv(A)