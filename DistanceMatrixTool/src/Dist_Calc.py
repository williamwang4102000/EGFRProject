import Bio.PDB
import numpy

pdb_code = "1m14"
pdb_filename = "1m14.pdb" #not the full cage!

def calc_residue_dist(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""
    diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
    return numpy.sqrt(numpy.sum(diff_vector * diff_vector))

def calc_dist_matrix(chain_one, chain_two) :
    """Returns a matrix of C-alpha distances between two chains"""
    answer = numpy.zeros((len(chain_one), len(chain_two)), numpy.float)
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            answer[row, col] = calc_residue_dist(residue_one, residue_two)
    return answer
structure = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filename)
model = structure[0]
from Bio.PDB import *
chain = model["A"]

def aa_residues(chain):
    aa_only = []
    for i in chain:
        if i.get_resname() in standard_aa_names:
            aa_only.append(i)
    return aa_only
AA_1 = aa_residues(model["A"])
dist_matrix = calc_dist_matrix(AA_1, AA_1)

import pylab
pylab.matshow(numpy.transpose(dist_matrix))
pylab.colorbar()
pylab.show()