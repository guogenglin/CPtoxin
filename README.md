# CPtoxin
A tool for Clostridium perfringens toxintyping

# External Dependencies
BLAST+

# Introduction
CPtoxin will fastly do the toxintyping of your input sequence file and generate a concisely result in both terminal and output file.

# Usage
Put the python script and CPtoxin file into the folder contains your sequence file

```
CPtoxin [-i] [-r] [-o] [-t] [--min_gene_cov] [--min_gene_id ]
Input and Output:
  -i, --input             Input FASTA file
  -r, --reference         Reference CPtoxin sequence file
  -o, --output            Output file
Parameters:
  -t, --threads           Threads to use for BLAST searches
  --min_gene_cov          Minimum percentage coverage to consider a single gene complete. [default: 80.0%]
  --min_gene_id           Minimum percentage identity to consider a single gene complete. [default: 70.0%]
```
# Quick usage
``` Python
python CPtoxin.py -i *.fasta 
