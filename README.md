# Center Star Multiple Sequence Alignment

Python implementation of the Center Star heuristic for multiple sequence alignment (MSA), with statistical analysis, matrix visualizations, and a PyQt5 GUI.

## Features

- DNA and protein sequence support
- Configurable Needleman-Wunsch scoring parameters
- Consensus sequence generation with configurable gap threshold
- MSA statistics: average pairwise identity, fully conserved columns, sum-of-pairs score
- Matrix visualizations: center election matrix, pairwise identity matrix, pairwise distance matrix
- FASTA file import and export of all results

## Usage

```bash
pip install -r requirements.txt
python main.py
```

## Tech Stack

Python, PyQt5, NumPy, Matplotlib

## Documentation

Full technical report including algorithm details, complexity analysis, and biological example results (Cytochrome C alignment across vertebrate species) is available in [`docs/report.pdf`](docs/report.pdf).
