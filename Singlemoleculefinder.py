'''
This Python script finds all bonds within a given .xyz structure file and uses them to find all organic spacer molecules within the structure.
It requires a replicatable unitcell in .xyz along with its Unit cell information (cell angle, cell lengths) as inputs.
It outputs .xyz files containing all unique molecules, all organic spacer molecules, and a .txt file containing all bonds within the structure.
'''
import numpy as np
from collections import defaultdict
import itertools

#Class containing all atom info
class Info:
    def __init__(self, index, type, x, y, z):
        self.index = index
        self.type = type
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)


# Cutoff dictionary
# These are the cutoff values for specific bonds within the structure in Angstrom,
# Atom pairs within this cutoff distance are considered bonded
cutoff_dict = {
    ('H', 'N'): 1.5,
    ('C', 'N'): 2,
    ('C', 'H'): 1.5,
    ('C', 'C'): 1.8,
    ('C', 'O'): 1.8,
    # ('I', 'Pb'): 4,
}


# Function for loading .xyz structure data
# Assigns atom info to each atom within the structure
def read_atoms(path):
    atoms = []
    with open(path, 'r') as file:
        lines = file.readlines()
        totatoms = int(lines[0].strip())
        for i, line in enumerate(lines[2:], start=1):
            parts = line.split()
            atom = Info(i, parts[0], parts[1], parts[2], parts[3])
            atoms.append(atom)
    return totatoms, atoms


# Bond finder function, uses values in cutoff_dict to find bonds
# Contains a debug print for each bond that is detected, displaying the distance at which it was found
def calcbonds(atoms):
    bonds = []
    for i, atom in enumerate(atoms):
        for j in range(i + 1, len(atoms)):
            neighbor = atoms[j]
            pair = (atom.type, neighbor.type)
            sorted_pair = tuple(sorted(pair))  # checken
            if sorted_pair in cutoff_dict:
                cutoff = cutoff_dict[sorted_pair]
                atompos = np.array([atom.x, atom.y, atom.z])
                neighborpos = np.array([neighbor.x, neighbor.y, neighbor.z])
                distance = np.linalg.norm(atompos - neighborpos)
                print(
                    f'Checking bond between {atom.index + 1} ({atom.type}) and {neighbor.index + 1} ({neighbor.type}): Distance = {distance:.3f}, Cutoff = {cutoff}')
                if distance < cutoff_dict[sorted_pair]:
                    bonds.append((atom.index + 1, neighbor.index + 1, atom.type, neighbor.type))
                    print(f'BOND - IDS: {atom.index + 1}, {neighbor.index + 1}, ATOMS: {atom.type}, {neighbor.type}')
    return len(bonds), bonds


#Molecule finder algorithim using Depth First Search (DFS) for graph traversal
'''
-----------------------------
DFS CODE START
-----------------------------
'''

#Functions defining each molecule as a graph
def add_edge(graph, u, v):
    graph[u].append(v)
    graph[v].append(u)
    print(graph)


def dfs(graph, node, visited, molecule):
    visited.add(node)
    molecule.append(node)
    for neighbor in graph[node]:
        if neighbor not in visited:
            dfs(graph, neighbor, visited, molecule)

def find_molecules(bonds):
    graph = defaultdict(list)
    for bond in bonds:
        #print(bond)
        add_edge(graph, bond[0], bond[1])

    visited = set()
    molecules = []
    for atom in graph:
        if atom not in visited:
            molecule = []
            dfs(graph, atom, visited, molecule)
            # molecule=tuple(sorted(molecule))
            molecules.append(molecule)
    return molecules

'''
-----------------------------
DFS CODE END
-----------------------------
'''


# Unique molecule finder function, takes find_molecule output to find each unique molecule within the structure

def uniqmol(molecules, atoms, unique):
    seenmol = set()
    uniqmols = []
    uniqmolname = []
    for molecule in molecules:
        moleculeid = []
        moleculeindex = []
        for i in molecule:
            j = i - 1
            if j < len(atoms):
                atomid = atoms[j].type
                atomindex = atoms[j].index
                moleculeid.append(atomid)
                moleculeindex.append(atomindex)
        moleculeid = tuple(sorted(moleculeid))
        # moleculeindex = tuple(sorted(moleculeindex))
        if moleculeid not in seenmol and unique == True:
            seenmol.add(moleculeid)
            uniqmolname.append(moleculeid)
            uniqmols.append(moleculeindex)
        elif unique == False:
            seenmol.add(moleculeid)
            uniqmolname.append(moleculeid)
            uniqmols.append(moleculeindex)
    return uniqmols, uniqmolname

#Atom replication functions nescesarry for handling Triclinic unit cells
import itertools


def lattice_vectors(a, b, c, alpha, beta, gamma):
    alpha, beta, gamma = np.radians([alpha, beta, gamma])
    v_x = np.array([a, 0, 0])
    v_y = np.array([b * np.cos(gamma), b * np.sin(gamma), 0])
    cx = c * np.cos(beta)
    cy = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    cz = np.sqrt(c ** 2 - cx ** 2 - cy ** 2)
    v_z = np.array([cx, cy, cz])
    return np.array([v_x, v_y, v_z])


def replicate(atoms, a, b, c, alpha, beta, gamma, replication_factors):
    # Get lattice vectors
    lattice = lattice_vectors(a, b, c, alpha, beta, gamma)
    # Extract fractional coordinates and types
    typearray = [atom.type for atom in atoms]
    atomcoords = np.array([[float(atom.x) / a, float(atom.y) / b, float(atom.z) / c] for atom in atoms])

    # Replication
    supercell_atoms = []
    count = 0
    for i, j, k in itertools.product(*[range(n) for n in replication_factors]):
        shift = i * lattice[0] + j * lattice[1] + k * lattice[2]
        for idx, frac_coord in enumerate(atomcoords):
            cartesian = frac_coord @ lattice + shift
            atom_type = typearray[idx]
            supercell_atoms.append(Info(count, atom_type, *cartesian))
            count += 1

    return supercell_atoms


#Function for outputting atom data into a .xyz file.
def convtoxyz(totatom, atoms, molecules, name):
    with open(f'{name}.xyz', 'w') as file:
        # read bonds array
        # total atoms in molecules
        tot = 0
        for molecule in molecules:
            for atom in molecule:
                tot = tot + 1
        totatom = tot
        file.write(f'{totatom}\n')
        file.write(f'convertedtestfile\n')
        for molecule in molecules:
            for i in molecule:
                if i < len(atoms):
                    aindex = i
                    file.write(
                        f'{atoms[aindex].type}       {atoms[aindex].x}   {atoms[aindex].y}    {atoms[aindex].z}\n')
                else:
                    break
        # extract info out of bonds array
        # write info in xyz format


# Debug function, outputting each atom in the bondlist into a .xyz file, useful for diagnosis to see if each atom is accounted for.
def bondtoxyz(totatom, atoms, bonds):
    with open('bondoutput.xyz', 'w') as file:
        # read bonds array
        # total atoms in molecules
        tot = 0
        file.write(f'{totatom}\n')
        file.write(f'convertedtestfile\n')
        lowestindex = 0
        writeindex = []
        texts = list()
        for bond in bonds:
            aindex = bond[0] - 1
            nindex = bond[1] - 1
            if 0 <= aindex < len(atoms) and 0 <= nindex < len(atoms):
                writeindex.append(aindex)
                writeindex.append(nindex)

        # writing
        writeindex = sorted(writeindex)
        print(writeindex)
        for i in range(len(writeindex)):
            j = i - 1
            if j < len(atoms):
                text = f'{atoms[j].type}       {atoms[j].x:.6f}   {atoms[j].y:.6f}    {atoms[j].z:.6f}'
                parts = text.split()
                formatted_text = f"{parts[0]:<8}{parts[1]:>12}{parts[2]:>12}{parts[3]:>12}"
                file.write(f'{formatted_text}\n')

                # file.write(f'{formatted_texta}\n')
                # file.write(f'{formatted_textn}\n')


# Output function for replicated structure obtained using replicate function
def write_replicatedxyz(atoms, filename):
    with open(filename, 'w') as f:
        # Write the number of atoms
        f.write(f"{len(atoms)}\n")
        # Write a comment line (optional)
        f.write("Replicated unit cell\n")
        # Write atom data
        for atom in atoms:
            f.write(f"{atom.type} {atom.x:.6f} {atom.y:.6f} {atom.z:.6f}\n")


#######


def main():
    path = 'aj40b.xyz' #Main input file containing replicatable unit cell.
    totatoms, atoms = read_atoms(path)

    '''
    Following are unit cell information inputs for several commonly found organic spacers
    The replicate function requires the cell lengths, and cell angles to correctly replicate the structure.
    The unit cell data can be obtained through the Mercury visaulization software.
    This is a nescesarry step for finding the correct bonds between atoms.
    The replication factor 
    '''
    # Replicate cell
    # BUTYLAMMONIUM (BA) unitcell
    # atoms = replicate(atoms, a=8.4280, b=8.986, c=26.233, replication_factors=(1, 1, 1), alpha= 90, beta= 90, gamma= 90)
    # 1,8-Diamino-3,6-dioxaoctane (DA-DOE) unit cell
    atoms = replicate(atoms, a=6.494, b=29.461, c=9.2666, replication_factors=(1, 1, 1), alpha=90, beta=91.777,
                      gamma=90)
    # Butylammonium (BA) + Methylammonium (MA) unitcell
    # atoms= replicate(atoms, a=8.9470, b=8.947, c=20.1757, replication_factors=(1, 1, 1), alpha=90, beta= 90, gamma= 90)
    # phenylethylammonium (PEA) unitcell
    # atoms = replicate(atoms, a=8.7389, b=8.7403, c=32.9952, replication_factors=(1, 1, 1), alpha=84.646, beta=84.657, gamma=89.643)

    #Debug structure output for replication function
    write_replicatedxyz(atoms, 'replicated.xyz')

    bonds = calcbonds(atoms)[1]
    print(f'Total bonds: {calcbonds(atoms)[0]}')

    #Print for all organic spacer molecules within the structure
    molecules = find_molecules(bonds)

    #Debug prints, useful for diagnosis
    #print(molecules)
    #print(uniqmol(molecules, atoms, unique=True)[0])

    #Debug print
    print(f"Longest molecule consists of: {len(max(molecules, key=len))} atoms")
    print(len(molecules))


    ##############
    # output files#
    ##############

    # Bond list, used in next script to find angles and dihedrals
    with open('bondlist.txt', 'w') as file:
        for bond in bonds:
            file.write(f'{bond}\n')

    # xyz file outputs for both the uniquemolecules and non unique molecules within the structure
    convtoxyz(totatoms, atoms, uniqmol(molecules, atoms, unique=True)[0], 'uniquemolecules')
    convtoxyz(totatoms, atoms, uniqmol(molecules, atoms, unique=False)[0], 'molecules')
    bondtoxyz(totatoms, atoms, bonds)


if __name__ == '__main__':
    main()

