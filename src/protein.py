import numpy as np
from parameters import atomtypes, hbondtypes

class PROTEIN:
    def __init__(self):
        self.coords = []
        self.struct = []
        self.para = []

    def ReadProt(self,filename, translation):
        with open(filename) as f:
             for line in f:
                 if line.startswith('ATOM'):
                    if line[13] != 'H' and (line[16]==' ' or line[16]=='A'):
                       self.struct.append(line)
                       self.coords.append(np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])]))
                       res = line[17:20]
                       if res == 'HIS': res = 'HSD'
                       atp = res + '_' + line[13:16].strip()
                       if atp in atomtypes:
                          if atp in hbondtypes:
                             if hbondtypes[atp] == 2:
                                self.para.append((None, atomtypes[atp], 3))
                             else:
                                self.para.append((None, atomtypes[atp], 2))  
                          else:
                             self.para.append((None, atomtypes[atp], None))
                       else:
                          self.para.append((None, 0, None))
                 elif line.startswith('END'):
                    break
                 elif line.startswith('TER'):
                    self.struct.append(line)
        self.coords = np.array(self.coords)
        self.coords -= translation
