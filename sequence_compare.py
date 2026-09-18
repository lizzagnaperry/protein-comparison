#compares the two protein sequences using the pairwise sequence alignment thingy,
#aka the satndard way to quantify how similar two sequences are

from Bio.Align import PairwiseAligner, substitution_matrices


#uses the BLOSUM62 substitution matrix aka it scores amino acid substitutions #
#by hwo often they occur betqeen evolutionarily related proteins and then returns best alignment
def align_sequences(seq_a: str, seq_b: str):
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.mode = "global"

    alignments = aligner.align(seq_a, seq_b)
    return alignments[0]


#how much/fractionw-wise are the exact amino acid matching for the aligned positions 
def percent_identity(alignment) -> float:
    aligned_a, aligned_b = alignment[0], alignment[1]
    matches = 0
    compared = 0
    for a, b in zip(aligned_a, aligned_b):
        if a == "-" or b == "-":
            continue  
        compared += 1
        if a == b:
            matches += 1
    if compared == 0:
        return 0.0
    return 100 * matches / compared


#returns a list of pairs for every postions where both sequences have a real residue
# (a residue is when theres no gap on either side)
def get_residue_correspondence(alignment) -> list[tuple[int, int]]:
    aligned_a, aligned_b = alignment[0], alignment[1]
    correspondence = []
    idx_a = 0
    idx_b = 0
    for a, b in zip(aligned_a, aligned_b):
        if a != "-" and b != "-":
            correspondence.append((idx_a, idx_b))
        if a != "-":
            idx_a += 1
        if b != "-":
            idx_b += 1
    return correspondence


#testing with myoglobin-is vs hemoglobin-alpha-ish toy sequences
if __name__ == "__main__":
    seq1 = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF"
    seq2 = "MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRF"

    alignment = align_sequences(seq1, seq2)
    print(alignment)
    print(f"Percent identity: {percent_identity(alignment):.1f}%")
    print(f"Aligned residue pairs: {len(get_residue_correspondence(alignment))}")