

## Code for DIRECT: Digital Microfluidics for Isolation-Free Shared Library Construction of Single-Cell DNA Methylome and Transcriptome

## Running the Code
### Input files
1. The fq.gz file from DIRECT dual-omics sequencing
2. Bismark reference file
### Processing steps
1. Specify the input sequence file and reference in process_bs.py
2. run `python -u process_bs.py`
3. run `python -u split_dual.py`
4. The transcriptome and methylationome were obtained separately.


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.