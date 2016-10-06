#!/bin/bash

FIRST_ARGUMENT=$1
SECOND_ARGUMENT=$2
THIRD_ARGUMENT=$3


perl project_motif_MW.pl -input_file=$FIRST_ARGUMENT -stats_file=$SECOND_ARGUMENT -motif_length=8 -max_mismatch=2 -population_size=200 -generations=100 -solutions=100


perl project_motif_MW.pl -input_file=$FIRST_ARGUMENT -stats_file=$SECOND_ARGUMENT -motif_length=10 -max_mismatch=3 -population_size=200 -generations=100 -solutions=100


perl project_motif_MW.pl -input_file=$FIRST_ARGUMENT -stats_file=$SECOND_ARGUMENT -motif_length=12 -max_mismatch=4 -population_size=200 -generations=100 -solutions=100



perl project_clustering.pl -output_instances=$THIRD_ARGUMENT -stats_file=$SECOND_ARGUMENT -input_sequences=$FIRST_ARGUMENT -max_similarity=0.5
