import numpy as np
import time
import math

import urllib

currURL = "https://files.rcsb.org/download/1UBQ.pdb"
response = urllib.request.urlopen(currURL)


class Residue:
    def __init__(self,name,id):
        self.name = name
        self.id = id
        self.atoms = []


class Atom:
    pass


class ContactMap:
    pass


class ContactMapCalculator:
    pass

vdW = {
    "CA": 1.90,
    "C": 1.75,
    "CH": 2.01,
    "CH2": 1.92,
    "CH2b": 1.91,
    "CH2ch": 1.88,
    "CH3": 1.92,
    "CHar": 1.82,
    "Car": 1.74,
    "CHim": 1.74,
    "Cco": 1.81,
    "Ccoo": 1.76,
    "SH": 1.88,
    "S": 1.94,
    "N": 1.71,
    "NH": 1.66,
    "NH1": 1.65,
    "NH2": 1.62,
    "NH21": 1.67,
    "NH31":  1.67,
    "O":  1.49,
    "Oco":  1.52,
    "Ocoo":  1.49,
    "OH":  1.54,
    "H2O": 1.68,
}

def get_atom_radius(atom):
    return water_radii

def get_atom_radius(atom):
    if atom in vdW:
        return vdW[atom]
    elif atom[0] in vdW:
        return vdW[atom[0]]
    else:
        return water_radii

lines = response.readlines()
water_radii=1.42 #Sobolev
water_radii=1.68 #Nussinov
protein = {}
for line in lines:
    line = line.decode('ascii')
    if(line[:4]=='ATOM'):
        chain=line[21]
        res_name=line[17:20]
        res=line[22:26]
        atom_name=line[13:16]
        x=float(line[31:38])
        y=float(line[39:46])
        z=float(line[47:54])
        if chain not in protein:
            protein[chain]={}
        if res not in protein[chain]:
            protein[chain][res] = Residue(res_name,res)
        protein[chain][res].atoms.append([atom_name,x,y,x])
        #print(line,chain,res,atom_name,x,y,z)


def create_interchain_contact_map(chain):
    residues = list(chain.keys())
    nlen = len(residues)
    contacts = []
    for r1 in range(nlen):
        residue1 = chain[residues[r1]]
        for r2 in range(r1,nlen):
            residue2 = chain[residues[r2]]
            in_contact = False
            for atom1 in residue1.atoms:
                r1 = get_atom_radius(atom1[0])
                for atom2 in residue2.atoms:
                    if atom1 == atom2:
                        continue
                    dx = atom1[1]-atom2[1]
                    dy = atom1[2]-atom2[2]
                    dz = atom1[3]-atom2[3]
                    dr = math.sqrt(dx*dx+dy*dy+dz*dz)
                    r2 = get_atom_radius(atom2[0])
                    if dr < water_radii+r1+r2:
                        in_contact = True
            if in_contact:
                contacts.append([residue1.id,residue2.id])
    return contacts



def make_surface(f):
    fiba = 0
    fibb = 1
    for i in range(f):
        fiba, fibb = fibb, fiba+fibb
    s = np.zeros([fibb, 3])
    k = np.array(range(fibb))
    theta = np.arccos(1.0-2.0*k/fibb)
    phi = np.mod(2.0*k*fiba, fibb)
    sin_theta = np.sin(theta)
    s[:, 0] = sin_theta*np.cos(phi)
    s[:, 1] = sin_theta*np.sin(phi)
    s[:, 2] = np.cos(theta)

    return s


surface = make_surface(10)


def create_interchain_contact_map(chain):
    residues = list(chain.keys())
    nlen = len(residues)
    contacts = []
    for r1 in range(nlen):
        residue1 = chain[residues[r1]]
        for r2 in range(nlen):
            residue2 = chain[residues[r2]]
            in_contact = False
            for atom1 in residue1.atoms:
                r1 = get_atom_radius(atom1[0])
                for atom2 in residue2.atoms:
                    if atom1 == atom2:
                        continue
                    dx = atom1[1]-atom2[1]
                    dy = atom1[2]-atom2[2]
                    dz = atom1[3]-atom2[3]
                    dr2 = dx*dx+dy*dy+dz*dz
                    r2 = get_atom_radius(atom2[0])
                    if dr2 < (water_radii+r1+r2)*(water_radii+r1+r2):
                        in_contact = True
            if in_contact:
                contacts.append([residue1.id,residue2.id])
    return contacts




start_time = time.time()
for chain in protein:
    create_interchain_contact_map(protein[chain])

print(time.time() - start_time)