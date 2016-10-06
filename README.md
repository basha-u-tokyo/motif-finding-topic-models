# A study on the application of topic models to motif finding algorithms

This repository contains the code of two motif finding methods based on topic models.

## CTM Approach

This program is coded in Perl and R and the installation of the following packages is required:

Perl:
* bioperl

R:
* NLP
* tm
* lda
* topicmodels

To run the program, the following command is used:

```perl ctm_motif_search.pl -input_file={input_file} -output_file={output_file} -stats_file={stats_file} -min_length={min_length} -max_length={max_length} -population_size={population_size} -number_of_motifs={number_of_motifs} -generations={generations} -instances_per_indiv={instances_per_indiv}
```

where the parameters are as follows:

- input_file: The input file with a set of sequences in Multi-FASTA format.
- output_file: The name of the output file (if it does not exist, it will be created)
- stats_file: The name of the an output file for statistical information regarding the execution (if it does not exist, it will be created)
- min_length: The minimum length of the k-mers searched
- max_length: The maximum length of the k-mers searched
- population_size: The size of the population for the Genetic Algorithm
- number_of_motifs: The number of motifs expected to be found in the input sequences
- generations: The number of generations for the Genetic Algorithm
- instances_per_indiv: The number of instances for each individual in the population

Example:

To test the program, you can run the following command for the example provided:

```perl ctm_motif_search.pl -input_file=example.fasta -output_file=test_result.txt -stats_file=test_stats.txt -min_length=6 -max_length=15 -population_size=50 -number_of_motifs=10 -generations=20 -instances_per_indiv=100
```

The output will be a list of motifs ordered by perplexity, each one of them with the format of the assessment created by Tompa et al at the [University of Washington](http://bio.cs.washington.edu/assessment/format.txt)

## Statistical GA Approach

This program is coded in Perl and R and the installation of the following packages is required:

Perl:
* bioperl
* PDL

R:
* NLP
* tm
* lda
* topicmodels

The Perl code uses also some modules that have some parts coded in C, so a C compiler is also required.

To run the program, the following command is used:

```./statGA.sh {input_file} {stats_file} {output_file}
```

where the parameters are as follows:

- input_file: The input file with a set of sequences in Multi-FASTA format.
- stats_file: The name of the an output file for statistical information regarding the execution (if it does not exist, it will be created)
- output_file: The name of the output file (if it does not exist, it will be created)

Example:

To test the program, you can run the following command for the example provided:

```./statGA.sh example.fasta test_stats.txt test_output.txt
```

The output will be a list of motifs ordered by perplexity, each one of them with the format of the assessment created by Tompa et al at the [University of Washington](http://bio.cs.washington.edu/assessment/format.txt)
