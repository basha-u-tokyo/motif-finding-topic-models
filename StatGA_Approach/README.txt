
This program is coded in Perl and R and the installation of the following packages is required:

Perl:
- bioperl
- PDL

R:
- NLP
- tm
- lda
- topicmodels

The Perl code uses also some modules that have some parts coded in C, so a C compiler is also required.

To run the program, the following command is used:

./statGA.sh {input_file} {stats_file} {output_file}

where the parameters are as follows:

- input_file: The input file with a set of sequences in Multi-FASTA format.
- stats_file: The name of the an output file for statistical information regarding the execution (if it does not exist, it will be created)
- output_file: The name of the output file (if it does not exist, it will be created)

Example:

To test the program, you can run the following command for the example provided:

./statGA.sh example.fasta test_stats.txt test_output.txt

The output will be a list of motifs ordered by perplexity, each one of them with the format of the assessment created by Tompa et al at the University of Washington:
http://bio.cs.washington.edu/assessment/format.txt