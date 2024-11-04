#!/usr/bin/env python

import numpy as np
import argparse, textwrap
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
from protac import PROTAC
from grid import Grid

def main(f_paras):
    paras = {
        'PROTACs': 'protac.sdf',
        'Warhead_anchor': 'e3.sdf',
        'Warhead_flex': 'poi.sdf',
        'Protein_anchor': 'e3.pdb',
        'Protein_flex': 'poi.pdb',
        'Output_protac': 'TC_protac.sdf',
        'Output_protein': 'TC_protein.pdb',
        'N_ini': 10000,
        'N_search': 900,
        'N_keep': 900,
        'N_processes': 1,
        'Interface': None
    }
    
    with open(f_paras) as f:
        for line in f:
            cells = line.strip().split(':')
            if cells[0] in paras:
               if cells[0].startswith('N'):
                  paras[cells[0]] = int(cells[1].split()[0])
               elif cells[0].startswith('Interface'):
                  paras[cells[0]] = np.reshape(np.array(cells[1].split(',')).astype(float), (-1,2)) 
               else:
                  paras[cells[0]] = cells[1].split()[0]
    
    w_anch = Chem.SDMolSupplier(paras['Warhead_anchor'])[0]
    w_flex = Chem.SDMolSupplier(paras['Warhead_flex'])[0]
    protacs = Chem.SDMolSupplier(paras['PROTACs'])
    fprotac = Chem.SDWriter(paras['Output_protac'])
    if paras['Output_protein'] == 'None':
       fprotein = None
       fpro_w = False
    else: 
       fprotein = open(paras['Output_protein'],'w')
       fpro_w = True

    #grid
    #anchoring protein
    grid_anchor = Grid(paras['Protein_anchor'],paras['Interface'],paras['N_processes'])
    #flexible protein
    coord_subs_var = w_flex.GetConformer().GetPositions()
    site_flex = []
    for i in range(3): 
          site_flex.append(np.min(coord_subs_var[:,i]) - 5)
          site_flex.append(np.max(coord_subs_var[:,i]) + 5)
    site_flex = np.reshape(site_flex, (-1,2))
    grid_flex = Grid(paras['Protein_flex'],site_flex,paras['N_processes'])
    ###
    
    for mol in protacs:
        if mol is not None and mol.HasSubstructMatch(w_anch) and mol.HasSubstructMatch(w_flex):
           protac = PROTAC()
           protac.init(mol,w_anch,w_flex,paras['Protein_flex'],paras['N_processes'],grid_anchor,grid_flex)
           protac.sample(paras['N_ini'],paras['N_search'])
           protac.output(fprotac,fprotein,paras['N_keep'],fpro_w)


if __name__ == "__main__":
   parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                    description=textwrap.dedent('''\
                                    TERNIFY: Efficient Sampling of PROTAC-Induced Ternary Complexes
                                    Hongtao Zhao, PhD'''))     

   parser.add_argument('-p', action='store', dest='f_paras', default='tcs.inp',
                        help=textwrap.dedent('''path to the parameter file\ndefault: tcs.inp'''))


   arg_dict = vars(parser.parse_args())
   
   main(**arg_dict)

