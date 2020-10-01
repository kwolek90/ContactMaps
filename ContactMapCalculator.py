import numpy as np
import math
import sys

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
                    r1 = get_atom_radius(atom1.name)+water_radii
                    for atom2 in residue2.atoms:
                        if atom1 == atom2:
                            continue
                        dx = atom1.x-atom2.x
                        dy = atom1.y-atom2.y
                        dz = atom1.z-atom2.z
                        dr = math.sqrt(dx*dx+dy*dy+dz*dz)
                        r2 = get_atom_radius(atom2.name)+water_radii
                        if dr < r1+r2:
                            in_contact = True
                if in_contact:
                    contacts.append([residue1.id,residue2.id])
        return contacts


class CSUContactMapCalculator(ContactMapCalculator):
    def calculate_contact_map(self, chain):
        residues = list(chain.keys())
        nlen = len(residues)
        contacts = []
        basic_surface = make_surface(14)
        slen = len(basic_surface)   
        res_box = []
        for r in range(nlen):
            residue = chain[residues[r]]
            xmin = sys.float_info.max
            xmax = -sys.float_info.max
            ymin = sys.float_info.max
            ymax = -sys.float_info.max
            zmin = sys.float_info.max
            zmax = -sys.float_info.max
            for atom in residue.atoms:
                xmin = min(xmin,atom.x)
                xmax = max(xmax,atom.x)
                ymin = min(ymin,atom.y)
                ymax = max(ymax,atom.y)
                zmin = min(zmin,atom.z)
                zmax = max(zmax,atom.z)
            res_box.append([(xmax+xmin)/2,(ymax+ymin)/2,(zmax+zmin)/2,(xmax-xmin)/2+7,(ymax-ymin)/2+7,(zmax-zmin)/2+7])

        for r1 in range(nlen):
            residue1 = chain[residues[r1]]
            b1=res_box[r1]
            unique_residue_contacts = set()
            for atom1 in residue1.atoms:
                r1 = get_atom_radius(atom1.name)+water_radii
                surf = basic_surface*r1+np.array([atom1.x,atom1.y,atom1.z])
                surf_closest_id=np.zeros(slen)-1
                surf_closest_dist=np.zeros(slen)+sys.float_info.max
                for r2 in range(nlen):
                    residue2 = chain[residues[r2]]
                    b2=res_box[r2]
                    if abs(b1[0]-b2[0]) > b1[3]+b2[3]:
                        continue
                    if abs(b1[1]-b2[1]) > b1[4]+b2[4]:
                        continue
                    if abs(b1[2]-b2[2]) > b1[5]+b2[5]:
                        continue

                    for atom2 in residue2.atoms:
                        if atom1 == atom2:
                            continue
                        dx = atom1.x-atom2.x
                        dy = atom1.y-atom2.y
                        dz = atom1.z-atom2.z
                        dr2 = dx*dx+dy*dy+dz*dz
                        r2 = get_atom_radius(atom2.name)+water_radii
                        if dr2<(r1+r2)*(r1+r2):
                            sdx = surf[:,0]-atom2.x
                            sdy = surf[:,1]-atom2.y
                            sdz = surf[:,2]-atom2.z
                            sdr2 = sdx*sdx+sdy*sdy+sdz*sdz
                            closer = sdr2 < surf_closest_dist
                            surf_closest_dist[closer] = sdr2[closer]
                            surf_closest_id[closer] = r2
                unique_residue_contacts.update(set(surf_closest_id))
            contacts += [[r1,r2] for r2 in unique_residue_contacts if r2 != -1 and r2 != r1]



                        
        return contacts
