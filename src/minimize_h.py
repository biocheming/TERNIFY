from rdkit import Chem
from rdkit.Chem import AllChem

def Minimize_H(mol):
    mol = Chem.AddHs(mol, addCoords = True)
    ffps = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
    ffps.SetMMFFVdWTerm(False)  
    ffps.SetMMFFEleTerm(False)
    ff = AllChem.MMFFGetMoleculeForceField(mol, ffps)
    for a in mol.GetAtoms():
        if a.GetAtomicNum() > 1:
           ff.MMFFAddPositionConstraint(a.GetIdx(), 0.0, 1.e9)
    ff.Minimize(maxIts=10000)
    return mol

