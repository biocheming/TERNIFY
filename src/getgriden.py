from parameters import paras

def GetGridEn(grid,coors,q):
    vdw = 0
    ele = 0
    hb = 0
    for i, coor in enumerate(coors):
        if coor[0]>= grid[1][0][0] and coor[0]<=grid[1][0][1] and coor[1]>grid[1][1][0] and coor[1]<=grid[1][1][1] and coor[2]>=grid[1][2][0] and coor[2]<=grid[1][2][1]:
           loc = (coor - grid[1][:,0])/paras['grid_space']
           loc = loc.astype(int)
           voxel = grid[0][loc[0]][loc[1]][loc[2]]
           vdw += voxel[0]
           ele += voxel[1]*q[i][1]
           if q[i][2] is not None:
              hb += voxel[q[i][2]]
    #print(vdw,ele,hb)          
    return vdw+ele+hb
