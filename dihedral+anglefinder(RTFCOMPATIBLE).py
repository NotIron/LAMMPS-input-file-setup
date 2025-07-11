import math
import numpy as np
from plotly.express import treemap
from collections import Counter as counter

def parse_xyz(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    atoms = []
    for line in lines[2:]:  # Skip the first two lines
        parts = line.split()
        atom = {
            'element': parts[0],
            'x': float(parts[1]),
            'y': float(parts[2]),
            'z': float(parts[3]),
        }
        atoms.append(atom)

    return atoms




def readbonds(file_path, atoms):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    bonded=[]
    for line in lines:
        parts = line.split()
        pair= parts[0], parts[1]
        if pair not in bonded:
            bonded.append(pair)
        intbonded = []
        for tup in bonded:
            intbonded.append((int(tup[0].strip('(),'))-1, int(tup[1].strip('(),'))-1))

    return intbonded

def find_angles(atoms, bonded,atomtypes):
    angles = []
    angletypes = []
    for i in range(len(atoms)):
        bonded_to_i = [j for j in range(len(atoms)) if (i, j) in bonded or (j, i) in bonded]
        '''
        if len(bonded_to_i)>0:
            print(atoms[bonded_to_i[0]]['element'])
        '''
        if len(bonded_to_i) > 1:
            for j in range(len(bonded_to_i)):
                for k in range(j+1, len(bonded_to_i)):
                    angles.append((bonded_to_i[j], i, bonded_to_i[k]))
                    angletype= atomtypes[bonded_to_i[j]], atomtypes[i], atomtypes[bonded_to_i[k]]
                    angletypes.append(angletype)

    # Count how many times each dihedral appears
    angle_counts = counter(angletypes)

    # Create final list with multiplicity
    anglesmultiplicity = [
        angletypes + (angle_counts[angletypes],)
        for angletypes in angle_counts
    ]
    return angles, anglesmultiplicity

def find_dihed(atoms, bonded,atomtypes):
    dihedrals = []
    num_atoms = len(atoms)
    dihedralids =[]

    # Create a bond lookup set for fast access
    bond_set = set(bonded) | {(j, i) for i, j in bonded}

    for i in range(num_atoms):
        for j in range(num_atoms):
            if (i, j) not in bond_set:
                continue
            for k in range(num_atoms):
                if k in {i, j} or (j, k) not in bond_set:
                    continue
                for l in range(num_atoms):
                    if l in {i, j, k} or (k, l) not in bond_set:
                        continue

                    dihedral = (
                        atomtypes[i],
                        atomtypes[j],
                        atomtypes[k],
                        atomtypes[l]
                    )
                    dihedralid = (i, j, k, l)

                    # Normalize the dihedral ID and type
                    dihedralid_rev = (l, k, j, i)
                    dihedral_rev = (atomtypes[l], atomtypes[k], atomtypes[j], atomtypes[i])

                    # Use the lexicographically smaller version
                    canonical_id = min(dihedralid, dihedralid_rev)
                    canonical_type = min(dihedral, dihedral_rev)

                    if canonical_id not in dihedralids:
                        dihedralids.append(canonical_id)
                        dihedrals.append(canonical_type)

    # Count how many times each dihedral appears
    dihedral_counts = counter(dihedrals)

    # Create final list with multiplicity
    dihedralsmultiplicity= [
    dihedral + (dihedral_counts[dihedral],)
    for dihedral in dihedral_counts
    ]

    return dihedralsmultiplicity, dihedralids

def write_dihedrals(dihedrals, dihedralids):
    with open('dihedralids.txt', 'w', newline='') as file:
        for dihedral in dihedralids:
            file.write(f"{atomtypes[dihedral[0]]} {atomtypes[dihedral[1]]} {atomtypes[dihedral[2]]} {atomtypes[dihedral[3]]}      {dihedral[0] + 1} {dihedral[1] + 1} {dihedral[2] + 1} {dihedral[3] + 1}\n")
    with open('dihedraltypes.txt', 'w', newline='') as file:
        for dihedral in dihedrals:
            file.write(f'{dihedral[0]} {dihedral[1]} {dihedral[2]} {dihedral[3]}\n')
def write_angles(angles, anglesmultiplicity):
    with open('angles.txt', 'w', newline='') as file:
        for angle in angles:
            #file.write(f"Angle between atoms {atomtypes[angle[0]]}({angle[0]}), {atomtypes[angle[1]]}({angle[1]}), {atomtypes[angle[2]]}({angle[2]}): {angle[3]:.2f} degrees\n")
            file.write(f"{atomtypes[angle[0]]} {atomtypes[angle[1]]} {atomtypes[angle[2]]}      {angle[0]+1} {angle[1]+1} {angle[2]+1}\n")
    with open('anglestypes.txt', 'w', newline='') as file:
        for angle in anglesmultiplicity:
            file.write(f"{angle[0]} {angle[1]} {angle[2]}\n")

def writebonds(bonds):
    seenbonds = set()
    with open('bondstype.txt', 'w', newline='') as file:
        for bond in bonds:
            print(bond)
            pair = atomtypes[bond[0]], atomtypes[bond[1]]
            #pair = tuple(sorted(pair))
            print(pair)
            if pair not in seenbonds:
                file.write(f"{atomtypes[bond[0]]}  {atomtypes[bond[1]]}\n")
                seenbonds.add(pair)
    with open('bondsall.txt', 'w', newline='') as file:
        for bond in bonds:
            file.write(f"{atomtypes[bond[0]]}  {atomtypes[bond[1]]}      {bond[0] + 1} {bond[1] + 1}\n")

def readrtf(path):
    with open(path, 'r', newline='') as file:
        lines = file.readlines()
    atomtypes=[]
    for line in lines:
        parts= line.split()
        atomtype=parts[4]
        atomtypes.append(atomtype)
    return atomtypes

# Example usage
file_path = 'molecules.xyz'
rtf_path='output.xyz'
bond_path='bondlistalligned.txt'
atoms = parse_xyz(file_path)
atomtypes=readrtf(rtf_path)
bonded = readbonds(bond_path, atoms)
#bonded= find_neigh_atoms(atoms)

#print (bonded2)
diheds = find_dihed(atoms, bonded,atomtypes)
angles = find_angles(atoms, bonded,atomtypes)
print(len(diheds[0]))
#print(diheds[1])
#print(diheds[0])
write_dihedrals(diheds[0], diheds[1])
write_angles(angles[0], angles[1])
print(bonded)
writebonds(bonded)


