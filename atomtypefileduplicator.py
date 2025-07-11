class Info:
    def __init__(self, index, type, x, y, z, parameter):
        self.index = index
        self.type = type
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.parameter = parameter
def parreferenceparse(filename):
    atoms = []
    with open(filename,'r') as f:
        lines = f.readlines()
        #totatoms = int(lines[0].strip())
        totatoms = 0
        for i, line in enumerate(lines[2:], start=1):
            parts = line.split()
            atom = Info(i, parts[0], parts[1], parts[2], parts[3], parts[4])
            atoms.append(atom)
            totatoms = totatoms + 1
    return atoms, totatoms
def xyzreferenceparse(filename, totatoms, atoms):
    with open(filename,'r') as f, open('output.xyz','w') as w:
        lines = f.readlines()
        for i, line in enumerate(lines[2:], start=0):
            cycle_value = (i) % (totatoms)
            line = line.strip()
            print(line)
            newline= f"{line}\t{atoms[cycle_value].parameter}\n"
            w.write(newline)


parname='atomtypes.xyz'
xyzname='molecules.xyz'
atomsdata= parreferenceparse(parname)
atoms=atomsdata[0]
totatoms = atomsdata[1]
print(totatoms)
xyzs= xyzreferenceparse(xyzname,totatoms,atoms)


