
import math

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

water_radii=1.42 #Sobolev
water_radii=1.68 #Nussinov


def get_atom_radius(atom):
    return water_radii


def get_atom_radius(atom):
    if atom in vdW:
        return vdW[atom]
    elif atom[0] in vdW:
        return vdW[atom[0]]
    else:
        return water_radii


class ContactMapCalculator:
    def calculate_contact_map(self, chain):
        raise NotImplementedError


class OverlapContactMapCalculator(ContactMapCalculator):
    def calculate_contact_map(self, chain):
        residues = list(chain.keys())
        nlen = len(residues)
        contacts = []
        for r1 in range(nlen):
            residue1 = chain[residues[r1]]
            for r2 in range(r1, nlen):
                residue2 = chain[residues[r2]]
                in_contact = False
                for atom1 in residue1.atoms:
                    r1 = get_atom_radius(atom1.name)
                    for atom2 in residue2.atoms:
                        if atom1 == atom2:
                            continue
                        dx = atom1.x-atom2.x
                        dy = atom1.y-atom2.y
                        dz = atom1.z-atom2.z
                        dr = math.sqrt(dx*dx+dy*dy+dz*dz)
                        r2 = get_atom_radius(atom2.name)
                        if dr < water_radii+r1+r2:
                            in_contact = True
                if in_contact:
                    contacts.append([residue1.id,residue2.id])
        return contacts


class CSUContactMapCalculator(ContactMapCalculator):
    def calculate_contact_map(self, chain):
        residues = list(chain.keys())
        nlen = len(residues)
        contacts = []
        for r1 in range(nlen):
            residue1 = chain[residues[r1]]
            for r2 in range(nlen):
                residue2 = chain[residues[r2]]
                in_contact = False
                for atom1 in residue1.atoms:
                    r1 = get_atom_radius(atom1.name)
                    for atom2 in residue2.atoms:
                        if atom1 == atom2:
                            continue
                        dx = atom1.x-atom2.x
                        dy = atom1.y-atom2.y
                        dz = atom1.z-atom2.z
                        dr2 = dx*dx+dy*dy+dz*dz
                        r2 = get_atom_radius(atom2.name)
                        if dr2 < (water_radii+r1+r2)*(water_radii+r1+r2):
                            in_contact = True
                if in_contact:
                    contacts.append([residue1.id, residue2.id])
        return contacts
