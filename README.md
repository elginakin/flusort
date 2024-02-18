# flusort

## Description
Influenza segment annotation and cleaning script for Pekosz Lab Influenza Surveillance. The current release is targed towards maintainers of influenza genome sequences received from Dr. Heba Mostafa's team. `flusort` is comprised of 2 python scripts designed to perform BLAST searches on FASTA sequences and append subtype information to their headers based on the BLAST results. It also groups sequences based on their subtype and completeness, providing a summary of the input sequences.

`flusort` accepts multifasta files containing H1, H3 and IBV (tested on B/Victoria) segment sequences with a strict header following >XXXXXXX_Segment# example: `JH11111_1` which would be Sample JH11111 segment 1 (PB2)

## Dependencies
flusort is written in python v3.X.X and requires [biopython](https://biopython.org/wiki/Download) (validated on v1.83) and [numpy](https://pypi.org/project/numpy/) (validated on v1.23.4). 

example installation using [`pip`](https://pip.pypa.io/en/stable/installation/)
```
pip install biopython numpy
```

## Arguments

### `flusort.py`

- -i,  --input_file: Path to the input FASTA file. (Required)
- -db,  --blast_db: Path to the custom BLAST database file. (Optional)
- -o,  --output_directory: Path to the output directory. Defaults to 'flusort_out' in the script's directory. (Optional)

#### Output

- output.fasta: Filtered sequences written to a single file if not split.
- segment_{segment_number}_sorted.fasta: Filtered sequences for each segment written to separate files if the split flag is used (-s, --split).

### `flusort_split.py`

- -i, --input: Path to the input FASTA file. (Required)
- -o, --output: Path to the output FASTA file. (Optional)
- --sequence_id: Sequence ID to filter on. (Optional)
- --segment_number: Segment number to filter on. (Optional)
- --influenza_type: Influenza type to filter on. (Optional)
- --ha_subtype: HA subtype to filter on. (Optional)
- --na_subtype: NA subtype to filter on. (Optional)
- --genome_completeness: Genome completeness to filter on. (Optional)
- --sequence: Sequence to filter on. (Optional)
- -s, --split: Split filtered sequences into separate files based on segment number. (Optional)

#### flusort_split output 

- blast_output.txt: BLAST search results.
- sequences_without_hits.fasta: Sequences without BLAST hits.
- original_dataframe.csv: Original DataFrame summarizing input sequences.
- grouped_dataframe.csv: Grouped DataFrame summarizing sequences based on subtype and completeness.
- fasta_with_subtype.fasta: FASTA file with subtype information appended to headers.

## Example Usage

first running `flusort.py`

```
python flusort.py -i ./example_input/flusort_example.fasta
```
All outputfiles will be  directory named `flusort_output` will be created unless specified.

Second (if required), run  `flusort_split.py` to filter your sequences by their desired traits. Here, I only want the complete genomes of H3s split by segment into differen fasta files. 

```
python3 ./scripts/flusort_split.py -i ./flusort_output/fasta_with_subtype -o output.fasta --ha_subtype H3 --genome_completeness Complete --split
```

