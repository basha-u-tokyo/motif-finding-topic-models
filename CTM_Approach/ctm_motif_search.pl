
#!/usr/bin/perl -w
use Bio::Perl;

use strict;
use warnings;
use Data::Dumper;
use List::Util qw(shuffle max min sum);
use Statistics::Basic qw(:all);
use POSIX qw(pow floor);
use Getopt::Long;
use List::Compare;
use Switch;
use Benchmark qw(:all);
use Math::BigInt;

use constant MAXINT => ~0;

# Measuring computation time
my $t0 = Benchmark->new;

# Check the validity of the parameters
my %args;
GetOptions(\%args,
           "input_file=s",
           "output_file=s",
           "stats_file=s",
           "min_length=s",
           "max_length=s",
           "population_size=s",
           "number_of_motifs=s",
	   "generations=s",
	   "instances_per_indiv=s",
) or die "Invalid arguments!";
die "Missing -input_file!" unless $args{input_file};
die "Missing -output_file!" unless $args{output_file};
die "Missing -stats_file!" unless $args{stats_file};
die "Missing -min_length!" unless $args{min_length};
die "Missing -max_length!" unless $args{max_length};
die "Missing -population_size!" unless $args{population_size};
die "Missing -number_of_motifs!" unless $args{number_of_motifs};
die "Missing -generations!" unless $args{generations};
die "Missing -instances_per_indiv!" unless $args{instances_per_indiv};

# Open the input file (FASTA FORMAT)
my $file = $args{input_file};
open (INPUTFILE, $file) or die "Could not open $args{input_file}. Please, use a valid file name";

# Getting the input file into an string an splitting it sequence by sequence
my @input = <INPUTFILE>;
my $input_string = join("", @input);
$input_string = substr($input_string, 1);
my $output = $args{output_file};
my $stats = $args{stats_file};

my @fasta_sequences = split("\n>", $input_string);
my @original_sequences;
foreach my $fasta (@fasta_sequences)
{
    my @lines = split("\n", $fasta);
    # Reading and removing the description from the list
    @lines = reverse(@lines);
    my $description = pop(@lines);
    @lines = reverse(@lines);
    my $seq = join("", @lines);
    $seq =~ s/\s//g;
    $seq = uc $seq;
    push @original_sequences, $seq;
}

# Close the input file
close (INPUTFILE);

my $init_pop_size = $args{population_size}; # individuals in the population
my $number_of_solutions = $args{number_of_motifs}; # number of solutions to report
my $number_of_motifs = $args{instances_per_indiv}; # number of instances for each individual in the population
my $motif_min_width = $args{min_length}; # minimum width of the motifs searched
my $motif_max_width = $args{max_length}; # maximum width of the motifs searched
my $mut_rate = 0.1;			# the mutation rate
my $generation_count = $args{generations};	        # run for this many generations
my $min_overrepresentation = roundup(length($original_sequences[0]) / 5); # minimum number of instances for a candidate motif

my $generation = 0;			# generation counter
my $pop_ref = [];                       # a reference to a population array

init_population($pop_ref, $init_pop_size);

do
{
 evaluate_fitness($pop_ref, \&fitness);

 # print out a generation summary line
 my @sorted_population = sort { $a->{fitness} <=> $b->{fitness} } @$pop_ref;

 printf "generation %d: size %d, least fit individual [%f], most fit individual [%f]\n",
     $generation,
     scalar @sorted_population,
     $sorted_population[-1]->{fitness},
     $sorted_population[0]->{fitness};

 open (STATS, ">>$stats") or die "Could not open. Please, use a valid file name";

 my @fitness_list = map { $_->{fitness} } @sorted_population;
 my $f_mean = mean @fitness_list;
 my $f_stddev = stddev @fitness_list;
 printf STATS "- Generation %d: Mean %f - Standard deviation %f\n", $generation, $f_mean, $f_stddev;

 close(STATS);

 survive($pop_ref);       # select survivors from the population
 select_parents($pop_ref);
 $pop_ref = recombine($pop_ref);        # returns a new population array reference

 # from this point on, we are working with a new generation in $pop_ref
 mutate($pop_ref, $mut_rate);	        # apply mutation to the individuals

} while ($generation++ < $generation_count); # run until we are out of generations


# Getting the final set of motifs
evaluate_fitness($pop_ref, \&fitness);

 my @fitness_list = map { $_->{fitness} } @$pop_ref;
 my $f_mean = mean @fitness_list;
 my $f_stddev = stddev @fitness_list;

 my @sorted_survivors = sort { $a->{fitness} <=> $b->{fitness} } (grep { ($_->{survived}) } @$pop_ref);

 my @best_individuals = @sorted_survivors[0..$number_of_solutions - 1];

my $index = 1;
foreach my $best_individual (@best_individuals)
{
    
    my @words;
    foreach my $instance (@{$best_individual->{motifs}})
    {
        push @words, find_motif($original_sequences[$instance->{seq}], $instance->{position}, $instance->{motif_width});
    }

    my %comb;
    @comb{@words} = ();
    my @uniqueWords = sort keys %comb;

    # Vocabulary file for CTM
    open (OUTPUT_CTM_VOCAB, ">./ctm_vocab.txt") or die "Could not open ctm_vocab.txt.";

    foreach my $instance (@uniqueWords)
    {
	printf OUTPUT_CTM_VOCAB "%s\n", $instance;
    }

    close(OUTPUT_CTM_VOCAB);

    open (OUTPUT_CTM_SEQS, ">./ctm_seqs.txt") or die "Could not open ctm_seqs.txt.";

    # Documents file for CTM
    foreach my $seq (@original_sequences)
    {
	my @words;
	my $i = 0;
	foreach my $word (@uniqueWords)
	{
	    my $count = approxPatternMatchCount($word, $seq, roundup(length($word)/4));
	    if ($count > 0) {
		push @words, { index => $i, count => $count };
	    }
	    $i++;
	}
	printf OUTPUT_CTM_SEQS "%d ", scalar(@words);
	foreach my $word (@words)
	{
	    printf OUTPUT_CTM_SEQS "%d:%d ", $word->{index}, $word->{count};
	}
	printf OUTPUT_CTM_SEQS "\n";
    }

    close(OUTPUT_CTM_SEQS);

    # Get the perplexity of the CTM
    my $number_of_instances = roundup(@uniqueWords / 5);
    my $perp = `Rscript ./finalMotif.R $number_of_instances`;

    my $perp_str = substr($perp, 4, length($perp) - 4);

    my $perplexity = sprintf("%.6f", $perp_str);


    open (INPUT_SOL, "./result_CTM.txt") or die "Could not open result_CTM.txt";

    # Getting the result file into a string and splitting it instance by instance
    my @results_file = <INPUT_SOL>;
    my $results_string = join("", @results_file);

    my @final_instances = split("\n", $results_string);

    close(INPUT_SOL);

    # Printing the results in an output file
    open (OUTPUT, ">>$output") or die "Could not open $args{output_file}. Please, use a valid file name";

    printf OUTPUT "> Motif %d - Perplexity %f\n", $index, $perplexity;
    
    my @final_positions;
    foreach my $instance (@final_instances)
    {
	foreach my $candidate (@{$best_individual->{motifs}})
	{
	    my $word = find_motif($original_sequences[$candidate->{seq}], $candidate->{position}, $candidate->{motif_width});
	    if ($instance eq $word)
	    {
		my $final_pos = $candidate->{position} - length($original_sequences[$candidate->{seq}]);
		my $complete_pos = $candidate->{seq}."-".$final_pos;
		if (!(grep { $_ eq $complete_pos } @final_positions))
		{
		    printf OUTPUT "%d,%d,%s\n", $candidate->{seq}, $final_pos, $word;
		    push @final_positions, $complete_pos;
		}
	    }
	}
    }

    printf OUTPUT "\n";

    close(OUTPUT);

    $index++;
}

my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);
open (STATS, ">>$stats") or die "Could not open. Please, use a valid file name";
print STATS "Elapsed time:",timestr($td),"\n";
close(STATS);

# Delete temporary files
unlink "result_CTM.txt", "ctm_seqs.txt", "ctm_vocab.txt";

sub init_population
{
 my $population = shift @_;
 my $pop_size = shift @_;

 # for each individual
 foreach my $id (1 .. $pop_size)
 {
     my @motifs;
     while (scalar(@motifs) < $number_of_motifs)
     {
	 my $seq = int(rand(scalar(@original_sequences)));
	 my $motif_width = int(rand($motif_max_width - $motif_min_width)) + $motif_min_width;
	 my $position = int(rand(length($original_sequences[$seq]) - $motif_width));
	 my $motif = find_motif($original_sequences[$seq], $position, $motif_width);
	 # Reverse Complement
	 my $rev_com = rev_complement($motif);
	 my $count = approxPatternMatchCount($motif, $original_sequences[$seq], roundup($motif_width/4)) +
	     approxPatternMatchCount($rev_com, $original_sequences[$seq], roundup($motif_width/4));
         # Generate shuffled motifs and calculate the count for them
	 my @motif_array = split("", $motif);
	 my @shuf_motif_array = shuffle(@motif_array);
	 foreach my $i (1..50) {
	     @shuf_motif_array = shuffle(@shuf_motif_array);
	 }
	 my $shuf_motif = join("", @shuf_motif_array);
         # Reverse Complement
	 my $rev_com_shuf = rev_complement($shuf_motif);
	 my $shuf_count = approxPatternMatchCount($shuf_motif, $original_sequences[$seq], roundup($motif_width/4)) +
	     approxPatternMatchCount($rev_com_shuf, $original_sequences[$seq], roundup($motif_width/4));
	 # Add to the population only if the difference between both scores is bigger or equal than $min_overrepresentation
	 if ($count - $shuf_count >= roundup($min_overrepresentation / ($motif_width ** 2)))  {
	     push @motifs, { seq => $seq, position => $position, motif_width => $motif_width };
	 }
     }
     push @$population, { motifs => \@motifs, survived => 1, parent => 0, fitness => MAXINT };
 }
}

sub evaluate_fitness
{
 my $population = shift @_;
 my $fitness_function = shift @_;

 foreach my $individual (@$population)
 {

     # set the fitness to the result of invoking the fitness function
     # on the individual
     if ($individual->{fitness} == MAXINT)
     {
	 $individual->{fitness} = $fitness_function->($individual);
     }
 }
}

sub survive
{
 my $population = shift @_;

 my @shuffled_population = shuffle(@$population);

 foreach my $i (1..(scalar(@shuffled_population)/2))
 {
     # Compare individuals in random pairs
     my $individual1 = $shuffled_population[$i-1];
     my $individual2 = $shuffled_population[$i];

     # The individual with better fitness survives
     if ($individual1->{fitness} <= $individual2->{fitness})
     {
	 $individual1->{survived} = 1;
	 $individual2->{survived} = 0;
	 # set the fitness to MAXINT for unfit individuals (so they won't procreate)
	 $individual2->{fitness} = MAXINT;
     } else
     {
	 # set the fitness to the result of invoking the fitness function
	 # on the individual's DNA
	 $individual2->{survived} = 1;
	 $individual1->{survived} = 0;
	 # set the fitness to MAXINT for unfit individuals (so they won't procreate)
	 $individual1->{fitness} = MAXINT;
     }
 }
}

sub select_parents
{
 my $population = shift @_;
 my $pop_size = scalar @$population;	# population size

 # create the weights array: select only survivors from the population,
 # then use map to have only the fitness come through
 my @weights = map { $_->{fitness} } grep { $_->{survived} } @$population;

 # if we have less than 2 survivors, the execution stops
 die "Population size $pop_size is too small" if $pop_size < 2;

 # we need to fill $pop_size parenting slots, to preserve the population size
 foreach my $slot (1..$pop_size)
 {
  my $index = sample(\@weights); # we pass a reference to the weights array here

  # do sanity checking on $index
  die "Undefined index returned by sample()"
   unless defined $index;
  die "Invalid index $index returned by sample()"
   unless $index >= 0 && $index < $pop_size;

  # increase the parenting slots for this population member
  $population->[$index]->{parent}++;

 }

}

sub recombine
{
    my $population = shift @_;
    my $pop_size = scalar @$population;	# population size
    my @parent_population;
    my @new_population = grep { $_->{survived} } @$population;

    my $total_parent_slots = 1;

    while ($total_parent_slots)
    {
	# find out how many parent slots are left
	$total_parent_slots = 0;
	$total_parent_slots += $_->{parent}>0 foreach @$population;

	last unless $total_parent_slots;

	# if we are here, we're sure we have at least one individual with parent > 0
	my $individual = undef;		# start with an undefined individual
	do
	{
	    # select a random individual
	    $individual = $population->[int(rand($pop_size))];
	    # individual is acceptable only if he can be a parent
	    undef($individual) unless $individual->{parent};
	} while (not defined $individual);

	push @parent_population, $individual;	# insert the individual in the parent population
	$individual->{parent}--;		# decrease the parenting slots by 1
    }


    while ($pop_size > scalar @new_population)
    {
	# select a random individual from the parent population (parent #2)
	my $parent1 = @parent_population[int(rand($pop_size))];
	my $parent2 = @parent_population[int(rand($pop_size))];
    
	my $child1 = { survived => 1, parent => 0, fitness => MAXINT };
	my $child2 = { survived => 1, parent => 0, fitness => MAXINT };

	my $point = int(rand($number_of_motifs));
	
	# Start crossover
	my @parent1_motifs = shuffle(@{$parent1->{motifs}}); # Shuffle the parent's motifs to combine them for the child
	my @parent2_motifs = shuffle(@{$parent2->{motifs}});

	my @child1_motifs = @parent1_motifs[0..($point - 1)];
	push @child1_motifs, @parent2_motifs[$point..($number_of_motifs - 1)];

	my @child2_motifs = @parent2_motifs[0..($point - 1)];
	push @child2_motifs, @parent1_motifs[$point..($number_of_motifs - 1)];

	$child1->{motifs} = \@child1_motifs;
	$child2->{motifs} = \@child2_motifs;

	push @new_population, $child1;		# the childs are now a part of the new generation
	push @new_population, $child2;
    }

    return \@new_population;
}

sub mutate
{
 my $population = shift @_;
 my $mut_rate   = shift @_;

 foreach my $individual (@$population)
 {
     # only mutate individuals if rand() returns more than mut_rate
     next if rand > $mut_rate;
     # mutate the individual by shifting a random number of elements
     my $motifs_mutated = int(rand($number_of_motifs));
     my @mutated_motifs = shuffle(@{$individual->{motifs}});
     foreach my $i (1..$motifs_mutated)
     {
	 my $new_position = $mutated_motifs[$i-1]->{position} + 
	     int(rand(2 * $mutated_motifs[$i-1]->{motif_width})) - $mutated_motifs[$i-1]->{motif_width};
	 while ($new_position > length($original_sequences[$mutated_motifs[$i-1]->{seq}]) - $mutated_motifs[$i-1]->{motif_width})
	 {
	     $new_position -= int(rand($mutated_motifs[$i-1]->{motif_width}));
	 }
	 while ($new_position < 0)
	 {
	     $new_position += int(rand($mutated_motifs[$i-1]->{motif_width}));
	 }

	 $mutated_motifs[$i-1]->{position} = $new_position; 
     }
     $individual->{motifs} = \@mutated_motifs;
     $individual->{fitness} = MAXINT;
 }
}

sub fitness
{
    my $individual = shift @_;

    # Vocabulary file for CTM
    my @words;
    foreach my $instance (@{$individual->{motifs}})
    {
        push @words, find_motif($original_sequences[$instance->{seq}], $instance->{position}, $instance->{motif_width});
    }

    my %comb;
    @comb{@words} = ();
    my @uniqueWords = sort keys %comb;

    open (OUTPUT_CTM_VOCAB, ">./ctm_vocab.txt") or die "Could not open ctm_vocab.txt.";

    foreach my $instance (@uniqueWords)
    {
	printf OUTPUT_CTM_VOCAB "%s\n", $instance;
    }

    close(OUTPUT_CTM_VOCAB);

    my $num_subsets = 1;
    my $length_subset = scalar(@original_sequences);
    my $length_last = scalar(@original_sequences);

    if (scalar(@original_sequences) > 10)
    {
	$num_subsets = floor(scalar(@original_sequences) / 10) + 1;
	$length_subset = roundup(scalar(@original_sequences) / $num_subsets);
	$length_last = scalar(@original_sequences) - ($num_subsets - 1) * $length_subset;
    }

    my $total_fitness = 0;

    for (my $k = 0; $k < $num_subsets; $k++)
    {
	open (OUTPUT_CTM_SEQS, ">./ctm_seqs.txt") or die "Could not open ctm_seqs.txt.";

	my @sequences;

	if ($k < $num_subsets - 1)
	{
	    @sequences = @original_sequences[$k * $length_subset..($k + 1) * $length_subset - 1];
	} else {
	    @sequences = @original_sequences[$k * $length_subset..($k * $length_subset) + $length_last - 1];
	}

	# Documents file for CTM
	foreach my $seq (@sequences)
	{
	    my @words;
	    my $i = 0;
	    foreach my $word (@uniqueWords)
	    {
		my $count = approxPatternMatchCount($word, $seq, roundup(length($word)/4));
		if ($count > 0) {
		    push @words, { index => $i, count => $count };
		}
		$i++;
	    }
	    printf OUTPUT_CTM_SEQS "%d ", scalar(@words);
	    foreach my $word (@words)
	    {
		printf OUTPUT_CTM_SEQS "%d:%d ", $word->{index}, $word->{count};
	    }
	    printf OUTPUT_CTM_SEQS "\n";
	}

	close(OUTPUT_CTM_SEQS);

	# Get the perplexity of the CTM
	my $perplexity = `Rscript ./getPerplexity.R`;

	my $fitness_str = substr($perplexity, 4, length($perplexity) - 4);

	my $fitness = sprintf("%.6f", $fitness_str);

	$total_fitness += $fitness;
    }

    # The fitness is the perplexity of the CTM
    return $total_fitness;

}


# Function to get the overrepresentation score
sub score
{
    my $motif = shift @_;

    # Generate shuffled motif and calculate the scores for it
    my @motif_array = split("", $motif);
    my @shuf_motif_array = shuffle(@motif_array);
    foreach my $i (1..50) {
	@shuf_motif_array = shuffle(@shuf_motif_array);
    }
    my $shuf_motif = join("", @shuf_motif_array);

    my $mismatches = roundup(length($motif)/4);

    my $score = 0;
    my $shuf_score = 0;

    foreach my $seq (@original_sequences)
    {
	$score += approxPatternMatchCount($motif, $seq, $mismatches);
	$shuf_score += approxPatternMatchCount($shuf_motif, $seq, $mismatches);
    }

    return $score - $shuf_score;

}

# Function to sample from an array of weighted elements
sub sample
{
 # get the reference to the weights array
 my $weights = shift @_ or return undef;
 # internal counting variables
 my ($count, $sample);

 for (my $i  = 0; $i < scalar @$weights; $i ++)
 {
  $count += $weights->[$i];
  $sample = $i if rand $count < $weights->[$i];
 }

 # return an index into the weights array
 return $sample;
}

# Function to round up a number
sub roundup {
    my $n = shift;
    return(($n == int($n)) ? $n : int($n + 1));
}

# Function to get the reverse complement
sub rev_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}

# Function that returns the motif given the absolute position
sub find_motif {
    my $seq = shift @_;
    my $pos = shift @_;
    my $width = shift @_;

    my $motif;

    $motif = substr ($seq, $pos, $width);
    return $motif;
}

sub approxPatternMatchCount
{
    my $pattern = shift @_;
    my $text = shift @_;
    my $d = shift @_;

    my $textLength = length($text);
    my $patternLength = length($pattern);

    my $matches = 0;

    for (my $i = 0; $i < ($textLength - $patternLength) + 1; $i++)
    {
	my $word = substr($text, $i, $patternLength);
	if (hammingDistance($pattern, $word) <= $d)
	{
	    $matches++;
	}
    }

    return $matches;
}

# Function to get the hamming distance between two motifs
sub hammingDistance
{ 

    return length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] );

}
