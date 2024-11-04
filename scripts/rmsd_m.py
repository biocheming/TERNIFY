#!/usr/bin/env python


import numpy as np
import re
import sys


def Seq(f1, f2):
    s1 = {}
    s2 = {}
    ss = {}
    with open(f1) as f:
         for line in f:
             if line.startswith('ATOM') and line[13:15] == 'CA':
                x = line[17:20]+line[22:26]
                s1[x] = None
             elif line.startswith('END'): break
    with open(f2) as f:
         for line in f:
             if line.startswith('ATOM') and line[13:15] == 'CA':
                x = line[17:20]+line[22:26]
                s2[x] = None
             elif line.startswith('END'): break
    for x in s1:
        if x in s2: ss[x] = None
    return ss  

def ReadStructRef(filename,ss):
    Coor_ref = []
    with open(filename) as f:
         for line in f:
             if line.startswith('ATOM') and line[13:15] == 'CA':
                x = line[17:20]+line[22:26]
                if x in ss:
                   Coor_ref.append(np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])]))
    return np.array(Coor_ref)

def ReadStructs(filename,ss):
    coor = []
    Coors = []
    with open(filename) as f:
         for line in f:
             if line.startswith('ATOM') and line[13:15] == 'CA':
                x = line[17:20]+line[22:26]
                if x in ss:
                   coor.append(np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])]))
             if line.startswith('END'):
                Coors.append(coor)
                coor = []
    return Coors

def rmsd(A, B):
    Coord = len(A[0])
    NAtom = len(A)
    cum = 0.0
    for i in range(NAtom):
        for j in range(Coord):
            cum += (A[i][j] - B[i][j])**2.0
    return np.sqrt(cum / NAtom)


if __name__ == "__main__":

   frmsd = open('rmsd.txt', 'w')
   ss = Seq('TC_protein.pdb',sys.argv[1]) 
   coor_ref = ReadStructRef(sys.argv[1],ss) 
   coors = ReadStructs('TC_protein.pdb',ss)

   clusters =[[None,None] for _ in range(len(coors))]

   ic = 0
   for i in range(len(coors)):
       if clusters[i][0] is None:
          clusters[i][0] = ic
       else:
          continue   
       for j in range(i+1,len(coors)):
           # 10 angstrom for protein ca, jcim 2020, 60, 5234
           if clusters[j][0] is None and rmsd(coors[i],coors[j])<=10:
              clusters[j][0] = ic
       ic += 1
   
   for i in range(len(coors)):
       clusters[i][1] = rmsd(coors[i],coor_ref)
       frmsd.write('%d  %d  %.1f \n'%(i+1,clusters[i][0]+1,clusters[i][1]))
 
   idx = np.argmin(np.array(clusters)[:,1])
   nwithin10 = 0
   high = ic + 1
   for c in clusters:
       if c[1] <= 10: 
           nwithin10 +=1
           if c[0] < high: high = c[0]
   print('RMSD Minimum: %d  %d  %.1f'%(idx+1,clusters[idx][0]+1,clusters[idx][1]))
   print('Number of poses within 10 Angstrom: %d out of %d'%(nwithin10,len(clusters)))
   if high > ic: print('No cluster having a pose within 10 Angstrom out of %d\n'%ic)
   else: print('Near-native cluster rank %d out of %d\n'%(high+1,ic))

