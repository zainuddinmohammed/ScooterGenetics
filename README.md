# 🧬 ScooterGenetics

**ScooterGenetics** is a lightweight Python library for comparing DNA sequences and building phylogenetic trees, tailored for genetic exploration and evolutionary visualization. Inspired by vermin studies (and scooters), this tool provides simple, readable outputs and educational visualizations for students and researchers alike.

---

## 🚀 Features

- ✅ **Compare DNA Sequences**
  - Align FASTA files using global or local alignment.
  - Get sequence alignment score, percent relatedness, and species/genus/family relationship inference.

- 🌳 **Build Phylogenetic Trees**
  - Merge multiple FASTA files.
  - Use ClustalW to create sequence alignments.
  - Generate and render Newick-format phylogenetic trees.

- 🖼 **Export Tree as Image + HTML**
  - Render circular or rectangular trees using [`ete3`](http://etetoolkit.org/).
  - Output a static `.png` plus a clickable `.html` preview.

---

## 🧪 Example Usage

```python
from ScooterGenetics import CompareSequences, BuildPhylogeneticTree, export_tree_ete

# Compare two sequences
CompareSequences("samples/A.fasta", "samples/B.fasta")

# Build tree from multiple FASTA files
tree = BuildPhylogeneticTree([
    "blattella/germanica/vermin_1.fasta",
    "parcoblatta/uhleriana/vermin_1.fasta",
    "parcoblatta/notha/vermin_1.fasta"
])

# Export visualized tree to HTML
export_tree_ete("tree.dnd", "phylo_tree.html")
```

---

## 📦 Requirements

- Python ≥ 3.8
- [Biopython](https://biopython.org/)
- [ete3](http://etetoolkit.org/)
- ClustalW2 (CLI tool must be in your system PATH)

Install Python dependencies:

```bash
pip install biopython ete3
```

---

## 🧠 Concepts Behind the Library

Want to learn more? Here are some relevant topics:
- FASTA Format & COI Barcoding
- Global vs Local Alignment
- Percent Identity vs Percent Coverage
- Homology & Phylogenetic Inference
- UPGMA & Distance Matrices
- Newick Tree Format
- Node Label Sanitization

---

## 🔧 Functions Overview

| Function | Description |
|----------|-------------|
| `CompareSequences(seq1, seq2)` | Aligns and compares two sequences |
| `BuildPhylogeneticTree(fasta_list)` | Merges sequences and builds a UPGMA tree |
| `NameVermin(path)` | Formats FASTA file path into readable label |
| `sanitize_newick_for_ete(path)` | Fixes internal node names for ETE |
| `export_tree_ete(newick_path)` | Creates a `.png` and `.html` visualization |

---

## 📄 License

MIT License — free to use, modify, and distribute.

![The Scooter in Question](https://i.imgur.com/CT8sDui.gif)
