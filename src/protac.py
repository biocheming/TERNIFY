from scipy.optimize import minimize
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdmolops, rdMolTransforms
import random, math
from multiprocessing import Pool
from protein import PROTEIN
from align import Align, Align_2
from getgriden import GetGridEn
from minimize_h import Minimize_H
from parameters import paras

class PROTAC:
    def __init__(self):
        self.protac = None
        self.rot_dihe = []
        self.E_intra_ref = None
        self.solu = None
        self.processes = 1
        self.grid_anchor = None
        self.grid_flex = None
        self.protein = PROTEIN()
        self.coord_subs_var = None
        self.idx = []
        self.q_flex = []
        self.q_anchor = []
        self.translation = None
        self.list_vdw = []
        self.list_dihe = []

    def init(self, protac, w_anch, w_flex, fpro_flex, processes, grid_anchor, grid_flex):
        self.processes = processes
        self.grid_anchor = grid_anchor
        self.grid_flex = grid_flex

        protac_0 = Minimize_H(protac)
        self.protac = protac

        conf_prot = self.protac.GetConformer()
        conf_anch = w_anch.GetConformer()
        conf_flex = w_flex.GetConformer()

        AllChem.ComputeGasteigerCharges(self.protac)

        #w_flex coords moved to origin
        self.coord_subs_var = conf_flex.GetPositions()
        self.translation = self.coord_subs_var.mean(axis=0)
        self.coord_subs_var -= self.translation
        #protein
        self.protein.ReadProt(fpro_flex, self.translation)

        #align on the bound conformation of the flexible warhead
        match_01 = self.protac.GetSubstructMatch(w_flex)
        cp_01 = None
        for i in range(len(match_01)): 
            for atom in self.protac.GetAtomWithIdx(match_01[i]).GetNeighbors():
                if atom.GetIdx() not in match_01:
                   cp_01 = i
                   break
        pairs = []
        pairs.append((match_01[cp_01],cp_01))
        for i in range(len(match_01)):
            if w_flex.GetBondBetweenAtoms(cp_01,i):
               pairs.append((match_01[i],i))  
        rdMolAlign.AlignMol(self.protac,w_flex,atomMap=pairs)
        for i,j in enumerate(match_01):
            coor = conf_flex.GetAtomPosition(i)
            conf_prot.SetAtomPosition(j,coor)
            self.idx.append(j)
        for i in range(len(self.protac.GetAtoms())):
            if i not in match_01:
               atom = self.protac.GetAtomWithIdx(i)
               q = atom.GetDoubleProp('_GasteigerCharge') + atom.GetDoubleProp('_GasteigerHCharge')
               self.q_anchor.append((i,q,None))

        #align on the bound conformation of the anchor warhead
        match_02 = self.protac.GetSubstructMatch(w_anch)
        cp_02 = None
        for i in range(len(match_02)): 
            for atom in self.protac.GetAtomWithIdx(match_02[i]).GetNeighbors():
                if atom.GetIdx() not in match_02:
                   cp_02 = i
                   break
        pairs = []
        pairs.append((match_02[cp_02],cp_02))
        for i in range(len(match_02)):
            if w_anch.GetBondBetweenAtoms(cp_02,i):
               pairs.append((match_02[i],i))  
        rdMolAlign.AlignMol(self.protac,w_anch,atomMap=pairs)
        for i,j in enumerate(match_02):
            coor = conf_anch.GetAtomPosition(i)
            conf_prot.SetAtomPosition(j,coor)
        for i in range(len(self.protac.GetAtoms())):
            if i not in match_02:
               atom = self.protac.GetAtomWithIdx(i)
               q = atom.GetDoubleProp('_GasteigerCharge') + atom.GetDoubleProp('_GasteigerHCharge')
               self.q_flex.append((i,q,None))
        
        #find linker
        linker = []
        warheads = match_01 + match_02
        for i in range(len(self.protac.GetAtoms())):
            if i not in warheads:
               linker.append(i)

        #find rotable bonds
        rbond = []
        #patt = Chem.MolFromSmarts('[!$(C(=[N,O,S])-!@[#7H,O,S])&!$([#7H,O,S]-!@C=[N,O,S])&!D1&!$(*#*)]-&!@[!D1&!$(*#*)]')
        patt = Chem.MolFromSmarts('[!$(C(=[N,O,S])-!@[#7H,O,S])&!$([#7H,O,S]-!@C=[N,O,S])&!D1]-&!@[!D1]')
        bonds = self.protac.GetSubstructMatches(patt)
        for a,b in bonds:
            if a in linker or b in linker:
               d1 = len(rdmolops.GetShortestPath(self.protac,a,match_02[cp_02])) if a != match_02[cp_02] else 0
               d2 = len(rdmolops.GetShortestPath(self.protac,b,match_02[cp_02])) if b != match_02[cp_02] else 0
               if d1<d2: rbond.append((a,b))
               else: rbond.append((b,a))                 
        #print(rbond)
        #find dihedrals
        for bond in rbond:
            na = []
            for atom in self.protac.GetAtomWithIdx(bond[0]).GetNeighbors():
                if atom.GetIdx() != bond[1]:
                   na.append(atom.GetIdx())
                   break
            for atom in self.protac.GetAtomWithIdx(bond[1]).GetNeighbors():
                if atom.GetIdx() != bond[0]:
                   na.append(atom.GetIdx())
                   break
            self.rot_dihe.append([na[0],bond[0],bond[1],na[1]])  
        #
        self.protac = protac_0        
        self.list(warheads,rbond)
        self.E_intra_ref = self.e_intra()
        self.protac = Minimize_H(protac)
        #

    def list(self,warheads,rbond):
        wh_h = {}
        for i in warheads:
            wh_h[i] = None
            for atom in self.protac.GetAtomWithIdx(i).GetNeighbors():
                if atom.GetAtomicNum() == 1: wh_h[atom.GetIdx()] = None              

        ffps = AllChem.MMFFGetMoleculeProperties(self.protac, mmffVariant='MMFF94s')
        #vdw pairs
        dist = AllChem.GetDistanceMatrix(self.protac)
        n = len(self.protac.GetAtoms())
        for i in range(n):
            for j in range(i,n):
                if dist[i,j] > 2:
                   if (i in wh_h) and (j in wh_h): continue
                   x = ffps.GetMMFFVdWParams(i,j)
                   self.list_vdw.append([i,j,x[3],1.07*x[2],0.07*x[2],1.12*x[2]**7,0.12*x[2]**7])
        #dihedrals
        for bond in rbond:
            for a in self.protac.GetAtomWithIdx(bond[0]).GetNeighbors():
                if a.GetIdx() != bond[1]:
                   for b in self.protac.GetAtomWithIdx(bond[1]).GetNeighbors():
                       if b.GetIdx() != bond[0]:
                          x = ffps.GetMMFFTorsionParams(self.protac,a.GetIdx(),bond[0],bond[1],b.GetIdx())
                          if x is not None:
                             self.list_dihe.append([a.GetIdx(),bond[0],bond[1],b.GetIdx(),x[1],x[2],x[3]]) 
        ###

    def output(self, w, fpro, nKeep,fpro_w):
        self.solu.sort(key=lambda tup: tup[1]) 
        for i,x in enumerate(self.solu):
            if i == nKeep: break
            score = x[1] #-self.E_intra_ref
            #if score > 0: continue
            for j in range(len(self.rot_dihe)):
                rdMolTransforms.SetDihedralDeg(self.protac.GetConformer(),*self.rot_dihe[j],x[0][j])
            self.protac.SetProp('score','%.2f'%(score))
            w.write(self.protac)

            #output proteins
            if fpro_w is False: continue
            pos = self.protac.GetConformer().GetPositions()
            ref = []
            for ii in self.idx:
                ref.append(pos[ii])
            ref = np.array(ref)

            coords_pro = Align(self.protein.coords, self.coord_subs_var, ref)

            fpro.write('REMARK Score: '+'%.2f'%(score) + '\n')
            index = 0
            fpro.write('MODEL'+'%9d'%(i+1) + '\n')
            for line in self.protein.struct:
                if line.startswith('ATOM') or line.startswith('HETATM'): 
                   fpro.write(line[:30]+'%8.3f'%coords_pro[index][0]+'%8.3f'%coords_pro[index][1]+'%8.3f'%coords_pro[index][2]+'\n')
                   index +=1
                else:
                   fpro.write(line)
            fpro.write('ENDMDL\n')

    def sample_single(self):
        dihe = np.array([random.uniform(-180, 180) for _ in range(len(self.rot_dihe))])
        ener = self.score(dihe)
        return (dihe,ener)

    def sample(self,ntotal=100, nsolu=100):
        with Pool(self.processes) as p:
             solu = p.starmap(self.sample_single,[[]]*ntotal) 
        solu.sort(key=lambda tup: tup[1]) 
        if nsolu < ntotal:
           self.solu = solu[:nsolu]
        else:
            self.solu = solu     
        self.search()

    def search_single(self,solu):
        #phase 1
        Temp = [60, 50, 40, 30, 20]
        alpha = 1.1
        best = solu
        for T in Temp:
            alpha = alpha * 0.9
            res = minimize(self.score,best[0],method='Powell',tol=0.01)
            if res.fun < best[1]: best = (res.x,res.fun)
            prev = best
            for _ in range(200):  
                delta = np.array([random.uniform(-20, 20) for _ in range(len(self.rot_dihe))]) * alpha
                current = prev[0] + delta       
                score = self.score(current)  
                if score < prev[1]: 
                   prev = (current,score)
                   if score < best[1]: best = prev
                elif np.exp((prev[1]-score)/T) > random.random():
                   prev = (current,score)       
        #phase 2
        res = minimize(self.score,best[0],method='Powell',tol=0.01)
        if res.fun < best[1]: best = (res.x,res.fun) 
        prev = best
        T = 10
        for _ in range(300):  
            delta = np.array([random.uniform(-5, 5) for _ in range(len(self.rot_dihe))]) 
            current = prev[0] + delta       
            score = self.score(current)  
            if score < prev[1]: 
               prev = (current,score)
               if score < best[1]: best = prev
            elif np.exp((prev[1]-score)/T) > random.random():
               prev = (current,score) 
        #final
        res = minimize(self.score,best[0],method='Powell',tol=0.01)
        if res.fun < best[1]: best = (res.x,res.fun)
        return best
    
    def search(self):
        with Pool(self.processes) as p:
            self.solu = p.map(self.search_single, self.solu)

    def e_intra(self) -> float:
        vdw = 0
        conf = self.protac.GetConformer()
        coors = [conf.GetAtomPosition(i) for i in range(len(self.protac.GetAtoms()))]
        for p in self.list_vdw:
            dist = coors[p[0]].Distance(coors[p[1]])
            vdw += p[2]*(p[3]/(dist+p[4]))**7*(p[5]/(dist**7+p[6])-2)
        dihe = 0
        for d in self.list_dihe:
            a = rdMolTransforms.GetDihedralRad(conf,d[0],d[1],d[2],d[3])
            dihe += 0.5*(d[4]*(1+math.cos(a)) + d[5]*(1-math.cos(2*a)) + d[6]*(1+math.cos(3*a)))
        return vdw+dihe

    def score(self,dihe) -> float:
        for j in range(len(self.rot_dihe)):
            rdMolTransforms.SetDihedralDeg(self.protac.GetConformer(),*self.rot_dihe[j],dihe[j])

        pos = self.protac.GetConformer().GetPositions()

        #intramoelcular torison+vdw
        ener = self.e_intra() - self.E_intra_ref
        if ener > paras['ub_strain']: ener = paras['ub_strain']

        #get interaction energy with the anchroing protein
        coors = []
        for q in self.q_flex:
            coors.append(pos[q[0]])
        ener += GetGridEn(self.grid_anchor,coors,self.q_flex)
        
        #protein-protein
        ref = []
        for ii in self.idx:
            ref.append(pos[ii])
        ref = np.array(ref)
        coords_pro = Align(self.protein.coords, self.coord_subs_var, ref)
        ener += GetGridEn(self.grid_anchor,coords_pro,self.protein.para)  
        
        #get interaction energy with the flexible protein
        coors = []
        for q in self.q_anchor:
            coors.append(pos[q[0]])
        coors = np.array(coors)
        coors = Align_2(coors,ref,self.coord_subs_var,self.translation)
        ener += GetGridEn(self.grid_flex,coors,self.q_anchor)
        #
        
        return ener

