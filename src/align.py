import numpy as np


def Kabsch(coord_var, coord_ref):
    covar = np.dot(coord_var.T, coord_ref)
    u, s, vh = np.linalg.svd(covar)
    d = (np.linalg.det(u) * np.linalg.det(vh)) < 0.0
    if d: u[:, -1] = -u[:, -1]
    R = np.dot(u, vh)
    return R

def Align(coords, coord_subs_var, coord_ref):
    center = coord_ref.mean(axis=0)
    coord_ref = coord_ref - center
    R = Kabsch(coord_subs_var,coord_ref)
    coords = np.dot(coords,R) + center    
    return coords

def Align_2(coords, coord_subs_var, coord_ref, translation):
    center = coord_subs_var.mean(axis=0)
    coord_subs_var = coord_subs_var - center
    coords = coords - center
    R = Kabsch(coord_subs_var,coord_ref)
    coords = np.dot(coords,R) + translation   
    return coords

