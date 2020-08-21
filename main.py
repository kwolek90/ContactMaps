import numpy as np
import time
import urllib.request

from ContactMapCalculator import *

currURL = "https://files.rcsb.org/download/1UBQ.pdb"
response = urllib.request.urlopen(currURL)


class Residue:
    def __init__(self, name, id, chain):
        self.name = name
        self.id = id
        self.chain = chain
        self.atoms = []

    def add_atom(self, atom):
        self.atoms.append(atom)


class Atom:
    def __init__(self, name, id, x, y, z):
        self.name = name
        self.id = id
        self.x = x
        self.y = y
        self.z = z


class ContactMap:
    pass


lines = response.readlines()
protein = {}
for line in lines:
    line = line.decode('ascii')
    if line[:4]=='ATOM':
        chain = line[21]
        res_name = line[17:20]
        res = line[22:26]
        atom_name = line[13:16]
        x = float(line[31:38])
        y = float(line[39:46])
        z = float(line[47:54])
        if res not in protein:
            protein[res] = Residue(res_name, res, chain)

        protein[res].add_atom(Atom(atom_name, 0, x, y, x))


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


start_time = time.time()
for calculator in ContactMapCalculator.__subclasses__():
    calculator().calculate_contact_map(protein)
    print(calculator.__name__,time.time() - start_time)

print(time.time() - start_time)
