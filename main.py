
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







start_time = time.time()
for calculator in ContactMapCalculator.__subclasses__():
    calculator().calculate_contact_map(protein)
    print(calculator.__name__,time.time() - start_time)

print(time.time() - start_time)
