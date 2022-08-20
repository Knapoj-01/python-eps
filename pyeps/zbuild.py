# Bus impedance matrix building, manual procedure.
# Code written by Knapoj Chaimanekorn based on the building algorithm 
# presented in Hadi Saadat's book.
# 12 August 2022

import numpy as np
def zbuild(zdata):
    
    # 0.. Initialization.
    k = int(np.max(zdata[:,:2])) # Check for size.
    zbus = np.zeros([k,k], dtype = complex)
    zdata = zdata[zdata[:,1].argsort()]
    zdata = zdata[zdata[:,0].argsort(kind = 'mergesort')] # Sort data

    processed_row = np.array([]) # Processed row
    existing_bus = np.array([0]) # Existing bus 

    # 1.. Add main branchs to the reference.
    i = 0 # Index
    while zdata[i,0] == 0:
        q = int(zdata[i,1]) - 1
        zpq = complex(zdata[i,2], zdata[i,3]) 
        zbus[q,q] = zpq # Add diagonal elements.
        existing_bus = np.append(existing_bus, zdata[i,1]) # Add bus number to the existing list.
        processed_row = np.append(processed_row, i)
        i += 1

    # 2.. Add new branches.
    addlist = np.array([], dtype = int)

    for ide, e in enumerate(zdata): # Iterate through the zdata matrix.
        if e[1] not in existing_bus: # Find new bus.
            addlist = np.append(addlist, ide)
            existing_bus = np.append(existing_bus, int(e[1]))
            processed_row = np.append(processed_row, ide)

    newbus = zdata[addlist] # Add new branches.
    newbus = newbus[newbus[:,1].argsort()]
    for e in newbus:
        # Redefine variables.
        p, q = int(e[0])-1, int(e[1]) - 1
        zpq = complex(e[2], e[3]) 
        for j in range(0,q+1): 
            zbus[j,q] = zbus[j,p]
            zbus[q,j] = zbus[p,j]
        zbus[q,q] = zbus[p,p] + zpq

    # 3.. Add cotree link between existing buses.
    for ide, e in enumerate(zdata):
        if ide not in processed_row:
            p, q = int(e[0])-1, int(e[1])-1
            zpq = complex(e[2], e[3])
            z1 =  np.reshape(zbus[:, p], (1,-1))

            if q > 0: # If q is not the reference bus.
                z2 = np.reshape(zbus[:, q], (1,-1))
                delz = z2-z1
                zll = zpq + zbus[p,p] + zbus[q,q] - 2*zbus[p,q]

            else: 
                delz = -z1
                zll = zpq + zbus[p,p]
            zbus = zbus - (1/zll)*np.matmul(delz.T,delz)
    return zbus