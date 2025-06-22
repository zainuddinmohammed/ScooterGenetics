from Bio import SeqIO, Phylo
from Bio.Align import PairwiseAligner
from Bio.Align import MultipleSeqAlignment
import subprocess
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from ete3 import Tree, TreeStyle
import re


def CompareSequences(seq1path, seq2path, file_type="fasta", show_sequences=True, show_percent=True, use_local=False):
    # Load sequences
    seq1_record = list(SeqIO.parse(seq1path, file_type))[0]
    seq2_record = list(SeqIO.parse(seq2path, file_type))[0]
    seq1 = seq1_record.seq
    seq2 = seq2_record.seq

    # Homology warning
    if "(COI)" not in seq1_record.description or "(COI)" not in seq2_record.description:
        print("⚠️ Warning: These sequences may not be homologous (e.g., COI vs COII)\n")

    # Aligner setup
    aligner = PairwiseAligner()
    aligner.mode = 'local' if use_local else 'global'
    alignments = aligner.align(seq1, seq2)
    top_alignment = alignments[0]

    # Print alignment if requested
    if show_sequences:
        print(top_alignment)
        print(f"\nAlignment Score: {top_alignment.score}")

    # Combine aligned blocks into full sequences
    aligned_str1 = ''.join([str(seq1[start:end]) for start, end in top_alignment.aligned[0]])
    aligned_str2 = ''.join([str(seq2[start:end]) for start, end in top_alignment.aligned[1]])

    # Alignment coverage
    aligned_region_length = max(len(aligned_str1), len(aligned_str2))
    longest_seq_length = max(len(seq1), len(seq2))
    coverage = (aligned_region_length / longest_seq_length) * 100 if longest_seq_length > 0 else 0

    if show_percent:
        print(f"\nPercent Related: {coverage:.2f}%")

    # Genetic relationship suggestion (thresholds can be tuned)
    if coverage > 98:
        relationship = "✅ Likely same species"
    elif coverage > 93:
        relationship = "🔷 Likely same genus"
    elif coverage > 85:
        relationship = "🟡 Possibly same family"
    else:
        relationship = "❌ Likely different family"

    print(f"\nRelationship Inference: {relationship}")


def BuildPhylogeneticTree(fasta_files, tree_file="tree.dnd", clustalw_path="clustalw2"):
    """
    Build a phylogenetic tree from a list of fasta files.

    Parameters:
    - fasta_files: List of paths to FASTA files
    - alignment_file: Output path for the alignment file
    - tree_file: Output path for the resulting phylogenetic tree in Newick or XML format
    - clustalw_path: Path to ClustalW executable (e.g., 'clustalw2')

    Returns:
    - tree: The generated phylogenetic tree
    """
    merged_fasta = "merged_sequences.fasta"
    with open(merged_fasta, "w") as outfile:
        for f in fasta_files:
            for record in SeqIO.parse(f, "fasta"):

                pretty_name, safe_id = NameVermin(f)
                record.id = safe_id
                record.name = pretty_name
                record.description = pretty_name

                SeqIO.write(record, outfile, "fasta")

    subprocess.run([clustalw_path, "-infile=" + merged_fasta], check=True)

    alignment = SeqIO.parse("merged_sequences.aln", "clustal")
    alignment = MultipleSeqAlignment(alignment)

    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment)

    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)

    Phylo.write(tree, tree_file, "newick")
    Phylo.draw_ascii(tree)

    return tree


def NameVermin(path):
    parts = re.split(r'[/.]', path)
    capitalized = [part.replace('_', ' ').capitalize() for part in parts[:-1]]
    return ' '.join(capitalized), '_'.join(part.capitalize() for part in parts[:-1])


def sanitize_newick_for_ete(path):
    with open(path, "r") as f:
        tree_str = f.read()

    import re
    sanitized = re.sub(r'\bInner\d+\b(?=:\d+)', '', tree_str)

    return sanitized


def export_tree_ete(newick_path="tree.dnd", html_file="ete_tree.html"):
    tree_str = sanitize_newick_for_ete(newick_path)
    t = Tree(tree_str)
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.mode = "r"

    t.render("tree.png", w=600, tree_style=ts)
    with open(html_file, "w") as f:
        f.write(f"<html><body><h2>Phylogenetic Tree</h2><img src='tree.png'></body></html>")
