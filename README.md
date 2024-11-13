# FMIndex: A Python Tool for FM Index Construction and Pattern Matching

## Overview
`FMIndex` is a Python tool designed to build FM Index components from a DNA sequence, stored in a FASTA format file.
The **FM Index** is a space-efficient data structure, designed for searching a specific substring within large texts, making it highly applicable in genomics for DNA sequence analysis. The different components that constitute it allow for efficient pattern matching, by supporting backward searching, where queries go across the index from the pattern's end. Additionally, the FM Index can be compressed, in order to reduce storage space. (Ferragina and Manzini, 2000), (Langmead, 2013)
The **different components** of the tool are:
- Burrows-Wheeler Transform (BWT)
- Suffix Array
- Occurrence Table
- LF Mapping
- Optional Run-Length Encoding (RLE) for BWT compression

## Features
- **FM Index Components**: Generates BWT, suffix array, occurrence table, and LF mapping.
- **Pattern Search**: Supports backward search within the FM Index for pattern matching.
- **Run-Length Encoding (RLE)**: Optional compression for BWT to save storage space.

## Requirements
- Python 3.6 or higher.

## Installation
Clone this repository:
```bash
git clone https://github.com/username/FMIndex.git
cd FMIndex
```

## Usage
The script requires a DNA sequence in FASTA format. You can also provide an optional pattern to search within the sequence and an option to apply run-length encoding.
**Command line arguments**:
- `fasta_file`: is the path to the FASTA file containing the DNA sequence;
- `--pattern`: (optional) is the pattern to search within the sequence;
- `--rle`: (optional) to apply run-length encoding (RLE) to the BWT for compression.
In order to construct **FM Index with both RLE and pattern search**:
```bash
python fm_index_script.py sample_dna.fasta --rle --pattern AGC
```

## Output
The script will create an `output` directory with the following files:
- `bwt.txt`: contains the Burrows-Wheeler Transform of the sequence;
- `suffix_array.txt`: contains the suffix array;
- `occurrence_table.txt`: contains the cumulative counts of each character in the BWT;
- `lf_map.txt`: contains the lf mapping array;
- `bwt_rle.txt`: contains the run-length encoded BWT (only if `--rle` is specified).
If a pattern is provided, the script will also display the positions in the sequence where the pattern is found.

## Version
Code created and tested on:
Python 3.11.7

## License
This project is licensed under the MIT License.






