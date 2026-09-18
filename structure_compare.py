#compares the 3d structures using sequence-guided superposition
#aka use the result from the sequence alignement to find which residues correspond between the two portiens
#then use the biopythons superimposer to fidn the best fit translation that overlays one set of atoms onto the other
#reports the RMSD (root mean square deviation) in Angstroms
#PSA low RMSD means that the structure are more similar in the regiosn that align

from Bio.PDB import PDBParser, Superimposer

from sequence_compare import align_sequences, get_residue_correspondence, percent_identity

#parses a PDB file and returns list of Cα atoms in their reisdue order
def get_ca_atoms(pdb_path: str):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    chain = next(structure[0].get_chains())

    ca_atoms = []
    for residue in chain:
        if "CA" in residue:
            ca_atoms.append(residue["CA"])
    return ca_atoms


#aligns the two sequences to find the corresponding residues
#then pulls out the matching Cα atom for each corresponding
#then superimpose one set onto the other and compute the RMSD
#retuns a dict with RMSD, redidue pairs computed, and underlying sequence identity
def compare_structures(pdb_path_a: str, seq_a: str, pdb_path_b: str, seq_b: str) -> dict:
    alignment = align_sequences(seq_a, seq_b)
    correspondence = get_residue_correspondence(alignment)

    ca_atoms_a = get_ca_atoms(pdb_path_a)
    ca_atoms_b = get_ca_atoms(pdb_path_b)

    matched_a = []
    matched_b = []
    for idx_a, idx_b in correspondence: 
        #prevents any off by one mismatch between sequences 
        #used for aligment and the residues present in the PDB file
        if idx_a < len(ca_atoms_a) and idx_b < len(ca_atoms_b):
            matched_a.append(ca_atoms_a[idx_a])
            matched_b.append(ca_atoms_b[idx_b])

    if len(matched_a) < 3:
        #superimposer needs at;east 3 points to define a meaningful 3d rotation
        #
        raise ValueError(
            "Not enough corresponding residues to compare structures "
            f"(found {len(matched_a)}, need at least 3)"
        )

    superimposer = Superimposer()
    superimposer.set_atoms(matched_a, matched_b)

    return {
        "rmsd_angstrom": round(superimposer.rms, 2),
        "compared_residue_pairs": len(matched_a),
        "sequence_identity_percent": round(percent_identity(alignment), 1),
    }


#testingggggggg
if __name__ == "__main__":
    import sys
    if len(sys.argv) != 5:
        print("Usage: python structure_compare.py <pdb_a> <seq_a> <pdb_b> <seq_b>")
        sys.exit(1)
    result = compare_structures(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    print(result)