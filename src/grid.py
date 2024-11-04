import numpy as np
from parameters import atomtypes, hbondtypes, paras
from multiprocessing import Pool

def CalcEn(atoms,coor,pocket,value):
    if value[0] is None:
       x = (coor[0]+0.5)*paras['grid_space'] + pocket[0][0]
       y = (coor[1]+0.5)*paras['grid_space'] + pocket[1][0]
       z = (coor[2]+0.5)*paras['grid_space'] + pocket[2][0]
       vdw = 0
       elec = 0
       for atom in atoms:   
           dist_square = (x-atom[0][0])**2 + (y-atom[0][1])**2 + (z-atom[0][2])**2
           if dist_square <= paras['nb_cutoff_2']:                 
              elec += atom[1]/dist_square
              temp = paras['rmin_6']/dist_square**3       
              vdw += paras['eps']*temp*(temp-2)
       value[0] = vdw
       value[1] = elec  
    return (coor,value)

def Grid(fpro,pocket,processes):
    eb_clash = int(paras['dist_clash']/paras['grid_space'])
    ub_hbond = int(paras['ub_hbond_dist']/paras['grid_space'])
    dist_clash_2 = paras['dist_clash']**2
    grid_space_2 = paras['grid_space']**2

    nbsize = np.copy(pocket)
    nbsize[:,0] -= paras['nb_cutoff']
    nbsize[:,1] += paras['nb_cutoff']

    atoms = []
    dim = (pocket[:,1] - pocket[:,0])/paras['grid_space']
    dim = np.array(dim).astype(int) + 1
    grid = [[[[None,0,0,0] for _ in range(dim[2])] for _ in range(dim[1])] for _ in range(dim[0])]

    with open(fpro) as f:
         for line in f:
             if line.startswith('ATOM'):
                if line[13] != 'H' and (line[16]==' ' or line[16]=='A'):
                   coor = np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])])
                   if coor[0]>= nbsize[0][0] and coor[0]<=nbsize[0][1] and coor[1]>nbsize[1][0] and coor[1]<=nbsize[1][1] and coor[2]>=nbsize[2][0] and coor[2]<=nbsize[2][1]:
                      res = line[17:20]
                      if res == 'HIS': res = 'HSD'
                      atp = res + '_' + line[13:16].strip()
                      if atp in atomtypes:
                         q = paras['elec_scaling']*atomtypes[atp] 
                      else:
                         q = 0
                      atoms.append((coor,q))
                   if coor[0]>= pocket[0][0] and coor[0]<=pocket[0][1] and coor[1]>pocket[1][0] and coor[1]<=pocket[1][1] and coor[2]>=pocket[2][0] and coor[2]<=pocket[2][1]:
                      loc = (coor - pocket[:,0])/paras['grid_space']
                      loc = loc.astype(int)
                      xmin = max(loc[0]-eb_clash, 0)
                      xmax = min(loc[0]+eb_clash, dim[0]-1)
                      ymin = max(loc[1]-eb_clash, 0)
                      ymax = min(loc[1]+eb_clash, dim[1]-1)
                      zmin = max(loc[2]-eb_clash, 0)
                      zmax = min(loc[2]+eb_clash, dim[2]-1)
                      for i in range(xmin,xmax+1):
                          for j in range(ymin,ymax+1):
                              for k in range(zmin,zmax+1):
                                  dist = ((i-loc[0])**2 + (j-loc[1])**2 + (k-loc[2])**2)*grid_space_2
                                  if dist <= dist_clash_2:
                                     grid[i][j][k][0] = paras['e_clash']
                      #hbond
                      if atp in hbondtypes:
                         xmin = max(loc[0]-ub_hbond, 0)
                         xmax = min(loc[0]+ub_hbond, dim[0]-1)
                         ymin = max(loc[1]-ub_hbond, 0)
                         ymax = min(loc[1]+ub_hbond, dim[1]-1)
                         zmin = max(loc[2]-ub_hbond, 0)
                         zmax = min(loc[2]+ub_hbond, dim[2]-1)
                         for i in range(xmin,xmax+1):
                             for j in range(ymin,ymax+1):
                                 for k in range(zmin,zmax+1):
                                     if grid[i][j][k][0] is None:
                                        dist = np.sqrt((i-loc[0])**2 + (j-loc[1])**2 + (k-loc[2])**2)*paras['grid_space']
                                        if dist < paras['ub_hbond_dist']:
                                           grid[i][j][k][hbondtypes[atp]] += (paras['ub_hbond_dist'] - dist)*paras['e_hbond']
                                           if grid[i][j][k][hbondtypes[atp]] < paras['e_hbond']:
                                              grid[i][j][k][hbondtypes[atp]] =  paras['e_hbond']                                  
    ###
    arguments = []
    for i in range(dim[0]):
        for j in range(dim[1]):
            for k in range(dim[2]):
                coor = [i,j,k]
                arguments.append([atoms,coor,pocket,grid[i][j][k]])
    with Pool(processes) as p:
            xyz = p.starmap(CalcEn, arguments)
    for xx in xyz:
        grid[xx[0][0]][xx[0][1]][xx[0][2]] = xx[1]
    #
    return (grid,pocket)


