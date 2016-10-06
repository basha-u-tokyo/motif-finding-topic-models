
#!/usr/bin/perl -w
use Bio::Perl;

use strict;
use warnings;
use Data::Dumper;
use List::Util qw(shuffle max min sum);
use Statistics::Basic qw(:all);
use String::Approx;
use POSIX qw(pow);
use Bio::SeqIO;
use Bio::Matrix::Scoring;
use Bio::SimpleAlign;
use PDL::Basic;
use Getopt::Long;
use Statistics::Test::WilcoxonRankSum;
use List::Compare;
use Switch;
use Benchmark qw(:all);

# Measuring computation time
my $t0 = Benchmark->new;

# Check the validity of the parameters
my %args;
GetOptions(\%args,
           "input_file=s",
	   "stats_file=s",
           "motif_length=s",
           "max_mismatch=s",
           "population_size=s",
	   "generations=s",
	   "solutions=s",
) or die "Invalid arguments!";
die "Missing -input_file!" unless $args{input_file};
die "Missing -stats_file!" unless $args{stats_file};
die "Missing -motif_length!" unless $args{motif_length};
die "Missing -max_mismatch!" unless $args{max_mismatch};
die "Missing -population_size!" unless $args{population_size};
die "Missing -generations!" unless $args{generations};
die "Missing -solutions!" unless $args{solutions};

# Open the input file (FASTA FORMAT)
my $file = $args{input_file};
open (INPUTFILE, $file) or die "Could not open $args{input_file}. Please, use a valid file name";

# Getting the input file into an string an splitting it sequence by sequence
my @input = <INPUTFILE>;
my $input_string = join("", @input);
$input_string = substr($input_string, 1);
my $output = "inner_output.txt";
my @fasta_sequences = split("\n>", $input_string);
my $total_sequence_length = 0;
my $super_seq = "";
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
    $total_sequence_length += length($seq);
    $super_seq .= $seq;
    push @original_sequences, $seq;
}

# Close the input file
close (INPUTFILE);
# individuals in the population
my $init_pop_size = $args{population_size};
my $motif_min_width = $args{motif_length};
my $motif_max_width = 12;
my $motif_max_mismatch = $args{max_mismatch};
my $mut_rate = 0.1;			# the mutation rate
my $min_fitness = 0.1;			# the minimum fitness for survival
my $max_thinness = 0.6;
my $max_mw = 0.1;
my $generation_count = $args{generations};	        # run for this many generations
my $generation = 0;			# generation counter
my $pop_ref = [];                       # a reference to a population array
my $subsequence_size = 500;
my $max_number_of_solutions = $args{solutions};
my $solutions = 0;

# Creating the subsequences
my @shuffled_sequences = shuffle(@original_sequences);
my $shuffled_super_seq = join("", @shuffled_sequences);
my @sequence_list = split(/(?(?{ pos() % $subsequence_size })(?!))/, $shuffled_super_seq);

# Positions used
my @positions_used;

# Number of sequences
my $seq_number = scalar (@sequence_list);

my $min_sequence_length = length ($sequence_list[0]);


init_population($pop_ref, $init_pop_size);

do
{
 evaluate_fitness($pop_ref, \&fitness);

 # print out a generation summary line
 my @sorted_population = sort { $a->{fitness} <=> $b->{fitness} } @$pop_ref;

 printf "generation %d: size %d, least fit individual [%f], most fit individual [%f]\n",
     $generation,
     scalar @sorted_population,
     $sorted_population[0]->{fitness},
     $sorted_population[-1]->{fitness};

 survive($pop_ref, $min_fitness);       # select survivors from the population
 select_parents($pop_ref);
 $pop_ref = recombine($pop_ref);        # returns a new population array reference

 # from this point on, we are working with a new generation in $pop_ref
 mutate($pop_ref, $mut_rate);	        # apply mutation to the individuals

 # Shuffling order for new subsequences
 if ($generation % $seq_number == $seq_number - 1)
 {
     my @shuffled_sequences = shuffle(@original_sequences);
     my $shuffled_super_seq = join("", @shuffled_sequences);
     @sequence_list = split(/(?(?{ pos() % $subsequence_size })(?!))/, $shuffled_super_seq);
 }
} while (($generation++ < $generation_count) && ($solutions < $max_number_of_solutions)); # run until we are out of generations

 # Printing the results in an output file
open (OUTPUT, ">>$output") or die "Could not open the output file";

evaluate_fitness($pop_ref, \&fitness);

 my @fitness_list = map { $_->{fitness} } @$pop_ref;
 my $f_mean = mean @fitness_list;
 my $f_stddev = stddev @fitness_list;
 if ($f_stddev == 0)
 {
     $f_stddev = 1;
 }

 foreach my $individual (@$pop_ref)
 {
  # set the fitness to the result of invoking the fitness function
  # on the individual's DNA

  $individual->{survived} = ($individual->{fitness} - $f_mean / $f_stddev >= $min_fitness) && 
      (thinness($individual, $pop_ref) < $max_thinness);

  # set the fitness to 0 for unfit individuals (so they won't procreate)
  $individual->{fitness} = 0 if ($individual->{fitness} - $f_mean / $f_stddev < $min_fitness);

 }

foreach my $individual (grep { ($_->{survived}) } @$pop_ref)
{
    my $position = $individual->{position};

    my $motif = find_motif($position);

    # Performing Mann-Whitney test
    my @scores = @{$individual->{mw_original}};
    my @shuf_scores = @{$individual->{mw_shuffled}};
    my $prob = 1;
    if (sum(@scores) == 0)
    {
	$prob = 1;
    }
    elsif (sum(@shuf_scores) == 0)
    {
	$prob = 0.05;
    }
    else
    {
	my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
	
	$wilcox_test->load_data(\@scores, \@shuf_scores);
	$prob = $wilcox_test->probability();
    }

    if ($prob <= ($max_mw/2))
    {
	# For consensus sequence
	my $aln = Bio::SimpleAlign->new();

	# Get all the motifs for the matrix
        my $id_counter = 0;
	foreach my $sequence (@original_sequences)
	{
	    foreach my $pos (0..(length($sequence) - $motif_min_width))
	    {
		my $new_motif = substr ($sequence, $pos, $motif_min_width);
		if (String::Approx::amatch($motif, ["i35%I0D0S$motif_max_mismatch"], $new_motif))
		{
		    my $loc_seq = new Bio::LocatableSeq(-seq =>$new_motif, -id => $id_counter);
		    $aln->add_seq($loc_seq);
		    $id_counter++;
		}
		elsif (String::Approx::amatch($motif, ["i35%I0D0S$motif_max_mismatch"], rev_complement($new_motif)))
		{
		    my $loc_seq = new Bio::LocatableSeq(-seq =>rev_complement($new_motif), -id => $id_counter);
		    $aln->add_seq($loc_seq);
		    $id_counter++;
		}
		if ($id_counter > 100)
		{
		    last;
		}
	    }   
	    if ($id_counter > 100)
	    {
		last;
	    }
	}

	# Print consensus sequences for filtering
	printf OUTPUT "%d %d %s %f\n", $position, $motif_min_width, $aln->consensus_string, $prob / ($motif_max_width/10);
    }

}

close(OUTPUT);

my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);
my $stats = $args{stats_file};
open (STATS, ">>$stats") or die "Could not open $args{stats_file}. Please, use a valid file name";
print STATS "Elapsed time:",timestr($td),"\n";
close(STATS);

sub init_population
{
 my $population = shift @_;
 my $pop_size = shift @_;

 # for each individual
 foreach my $id (1 .. $pop_size)
 {
     # insert an anonymous hash reference in the population array with the individual's data
     # In the end it must be an array of positions
     push @$population, { position => make_unique(int(rand($total_sequence_length - $motif_min_width - 1)), $population), survived => 1, parent => 0, fitness => 0, mw_original => [], mw_shuffled => [] };
 }
}

sub evaluate_fitness
{
 my $population = shift @_;
 my $fitness_function = shift @_;

 foreach my $individual (@$population)
 {

     # Getting the motif given by the position
     my $motif = find_motif($individual->{position});

     # set the fitness to the result of invoking the fitness function
     # on the individual's position
     # If the subsequences have not been reset, the fitness will be accumulative
     if ($generation % $seq_number == $seq_number - 1)
     {
	 $individual->{fitness} = $fitness_function->($motif, $individual);
     } 
     else
     {
	 $individual->{fitness} += $fitness_function->($motif, $individual);
     }

 }
}

sub survive
{
 my $population = shift @_;
 my $min_fitness = shift @_;

 my @fitness_list = map { $_->{fitness} } @$population;
 my $f_mean = mean @fitness_list;
 my $f_stddev = stddev @fitness_list;
 if ($f_stddev == 0)
 {
     $f_stddev = 1;
 }

 foreach my $individual (@$population)
 {
  # set the fitness to the result of invoking the fitness function
  # on the individual's DNA
  $individual->{survived} = ($individual->{fitness} - $f_mean / $f_stddev >= $min_fitness);

  # Performing Mann-Whitney test
  my @scores = @{$individual->{mw_original}};
  if (scalar(@scores) > 5)
  {
      my @shuf_scores = @{$individual->{mw_shuffled}};
      my $prob = 1;
      if (sum(@scores) == 0)
      {
	  $prob = 1;
      }
      elsif (sum(@shuf_scores) == 0)
      {
	  $prob = 0.05;
      }
      else
      {
	  my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
	  
	  $wilcox_test->load_data(\@scores, \@shuf_scores);
	  $prob = $wilcox_test->probability();
      }

      if ($prob > $max_mw)
      {
	  $individual->{survived} = 0;
	  $individual->{fitness} = 0;
      }

      if ((scalar(@scores) > 9) && ($prob <= ($max_mw/2)) && (thinness($individual, $pop_ref) < $max_thinness))
      {
          # Getting the output file into an string an splitting it solution by solution
	  my @solutions;
	  if (-e $output) {
	      open (OUTPUTFILE, $output);
	      my @output = <OUTPUTFILE>;
	      my $output_string = join("", @output);
	      $output_string = substr($output_string, 1);
	      my @solutions_raw = split("\n", $output_string);
	  
	      foreach my $solution (@solutions_raw)
	      {
		  my @info = split(" ", $solution);
		  # Getting position, consensus sequence and fitness
		  my $pos = $info[0];
		  my $cons = $info[1];
		  my $fitn = $info[2];
	      
		  push @solutions, { position => $pos, consensus => $cons, fitness => $fitn };
	      }
	  
	      # Close the input file
	      close (OUTPUTFILE);
	  }

	  # For consensus sequence
	  my $aln = Bio::SimpleAlign->new();

	  my $position = $individual->{position};

	  my $motif = find_motif($position);
	  open (OUTPUT, ">>$output") or die "Could not open $args{output_file}. Please, use a valid file name";
	  # Get all the motifs for the matrix
	  my $id_counter = 0;
	  foreach my $sequence (@original_sequences)
	  {
	      foreach my $pos (0..(length($sequence) - $motif_min_width))
	      {
		  my $new_motif = substr ($sequence, $pos, $motif_min_width);
		  if (String::Approx::amatch($motif, ["i35%I0D0S$motif_max_mismatch"], $new_motif))
		  {
		      my $loc_seq = new Bio::LocatableSeq(-seq =>$new_motif, -id => $id_counter);
		      $aln->add_seq($loc_seq);
		      $id_counter++;
		  }
		  elsif (String::Approx::amatch($motif, ["i35%I0D0S$motif_max_mismatch"], rev_complement($new_motif)))
		  {
		      my $loc_seq = new Bio::LocatableSeq(-seq =>rev_complement($new_motif), -id => $id_counter);
		      $aln->add_seq($loc_seq);
		      $id_counter++;
		  }
		  if ($id_counter > 100)
		  {
		      last;
		  }
	      }   
	      if ($id_counter > 100)
	      {
		  last;
	      }
	  }

	  # Print consensus sequences for filtering
	  printf OUTPUT "%d %d %s %f\n", $position, $motif_min_width, $aln->consensus_string,  $prob / ($motif_max_width/10);

	  $individual->{survived} = 0;
	  $individual->{fitness} = 0;

	  $solutions++;

	  close(OUTPUT);
      }

  }

  # set the fitness to 0 for unfit individuals (so they won't procreate)
  $individual->{fitness} = 0 if ($individual->{fitness} - $f_mean / $f_stddev < $min_fitness);

 }
}

sub select_parents
{
 my $population = shift @_;
 my $pop_size = scalar @$population;	# population size

 # create the weights array: select only survivors from the population,
 # then use map to have only the fitness come through
 my @weights = map { $_->{fitness} } grep { $_->{survived} } @$population;

 # if we have less than 2 survivors, we're in trouble
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
    
	my $child = { survived => 1, parent => 0, fitness => 0 };

	# One-point crossover
	my $point = int(rand($motif_min_width - 1)) + 1;
	
	# this is breeding!
	my $parent1_pos = $parent1->{position};
	my $parent2_pos = $parent2->{position};

	my $parent1_motif = find_motif($parent1_pos);
	my $parent2_motif = find_motif($parent1_pos);

	# Child mixing the end of the parent1 with the beginning of the parent2
	my $child_motif = substr($parent1_motif, $point, $motif_min_width - $point) . substr($parent2_motif, 0, $point);

	# First motif found with 50% mismatch
	$child->{position} = make_unique_child($child_motif, \@new_population);

	if ($child->{position} == -1)
	{
	    last;
	}

	push @new_population, $child;		# the child is now a part of the new generation
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
  # mutate the position by adding one
  $individual->{position} = make_unique($individual->{position} + int(rand($motif_min_width - 1)) + 1, $population);

  if ($individual->{position} == -1)
  {
      last;
  }
 }
}

sub fitness
{
    my $motif = shift @_;
    my $individual = shift @_;

    # Generate shuffled motif and calculate the scores for it
    my @motif_array = split("", $motif);
    my @shuf_motif_array = shuffle(@motif_array);
    foreach my $i (1..50) {
	@shuf_motif_array = shuffle(@shuf_motif_array);
    }
    my $shuf_motif = join("", @shuf_motif_array);

    # Calculate the scores for the original motif and the shuffle motif
    my $seq = $sequence_list[$generation % $seq_number];

    my $score = similar_word_number($motif, $seq);
    my $mw_orig_ref = $individual->{mw_original};
    push @$mw_orig_ref, $score;
    $individual->{mw_original} = $mw_orig_ref;

    my $shuf_score = similar_word_number($shuf_motif, $seq);
    my $mw_shuf_ref = $individual->{mw_shuffled};
    push @$mw_shuf_ref, $shuf_score;
    $individual->{mw_shuffled} = $mw_shuf_ref;

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

# Function to check if the individual already exists in the population and, if so,
# assign him a new position
sub make_unique
{
    my $position = shift @_;
    my $population = shift @_;
    my @positions = map { $_->{position} } @$population;

    if ($position >= ($total_sequence_length - $motif_min_width))
    {
	$position = 0;
    } elsif (is_inside_sequence($position) == 0)
    {
	$position += $motif_min_width;
    }
    my $count = 0;
    while ((($position ~~ @positions) || (is_valid_motif(find_motif($position))  == 0)) && ($count < $total_sequence_length - ($motif_min_width * scalar(@original_sequences))))
    {
	if ($position >= $total_sequence_length - $motif_min_width)
	{
	    $position = 0;
	} elsif (is_inside_sequence($position))
	{
	    $position += $motif_min_width;
	} else
	{
	    $position++;
	}
	$count++;
    }
 
    return $position;
}

# Function to check if the newborn child individual already exists in the population and, if so,
# assign him a new position
sub make_unique_child
{
    my $motif = shift @_;
    my $population = shift @_;
    my @positions_used = map { $_->{position} } @$population;
    
    my @positions;

    my $mismatches = 0;
    my $position = int(rand($total_sequence_length - $motif_min_width - 1));
    my $initial_position = $position;
    my @previous_positions;
    push @previous_positions, $position;
    my $pos_found = -1;
    my $offset = $position;
    while ($pos_found == -1)
    {
	$pos_found = String::Approx::aindex($motif, ["i50%I0D0S$mismatches"], substr($super_seq, $position + 1, $total_sequence_length - $motif_min_width - $position));
	if ($pos_found == -1) {
	    if (($position < $initial_position) || ($position == 0))
	    {
		$position = int(rand($total_sequence_length - $motif_min_width - 1));
		while ($position >= $total_sequence_length - $motif_min_width - 1)
		{
		    $position = int(rand($total_sequence_length - $motif_min_width - 1));
		}
		@previous_positions = [];
		push @previous_positions, $position;
		$initial_position = $position;
		$mismatches++;
	    } 
            else  {		
		$position = 0;
		push @previous_positions, $position;
	    }
	    $offset = $position;
	} else {
	    $offset += $pos_found + 1;
	    $pos_found += $position;
	} 

	if (($pos_found > -1) && (($pos_found ~~ @positions_used) || (is_inside_sequence($pos_found)  == 0) || (is_valid_motif(find_motif($pos_found))  == 0)))
	{
	    $pos_found = -1;
	    $position = $offset;
	    push @previous_positions, $position;
	    if ($position ~~ @previous_positions)
	    {
	        $position = int(rand($total_sequence_length - $motif_min_width - 1));
		while ($position >= $total_sequence_length - $motif_min_width - 1)
		{
		    $position = int(rand($total_sequence_length - $motif_min_width - 1));
		}
		@previous_positions = [];
		push @previous_positions, $position;
		$initial_position = $position;
		$mismatches++;
	    }
	    elsif ($position >= $total_sequence_length - $motif_min_width - 1)
	    {
		$position = 0;
		$offset = $position;
		push @previous_positions, $position;
	    }
	} 		
    }
 
    return make_unique($pos_found);
}

# Get all the positions with similar motifs
sub get_all_positions
{
    my $motif = shift @_;
    my $mismatches = shift @_;
    my $population = shift @_;

    my $pos_found = 0;
    my $position = 0;
    my $offset = $position;
    my @positions;
    while ($pos_found > -1)
    {
	$pos_found = String::Approx::aindex($motif, ["i50%I0D0S$mismatches"], substr($super_seq, $position + 1, $total_sequence_length - $motif_min_width - $position));
        $offset += $pos_found + 1;
	if ($pos_found > -1)
	{
	    push @positions, $pos_found;
	    $position = $offset;
	}
    }
 
    return @positions;
}

# Function to round up a number
sub roundup {
    my $n = shift;
    return(($n == int($n)) ? $n : int($n + 1));
}

# Function that returns the motif given the absolute position
sub find_motif {
    my $pos = shift;
    my $motif;

    $motif = substr ($super_seq, $pos, $motif_min_width);

    return $motif;
}

# Function to get the number of similar words
sub similar_word_number {
    my $motif = shift @_;
    my $sequence = shift @_;
    my @motif_list;
    my $similar_word_number = 0;
    
    foreach my $pos (0..(length($sequence) - $motif_min_width))
    {

	if ((String::Approx::amatch($motif, ["i35%I0D0S$motif_max_mismatch"], substr ($sequence, $pos, $motif_min_width))) ||
	    (String::Approx::amatch(rev_complement($motif), ["i35%I0D0S$motif_max_mismatch"], substr ($sequence, $pos, $motif_min_width))))
	{
	    $similar_word_number++;
	}
    }

    return $similar_word_number;
}

# Function to get the kurtosis of the population
sub dist_kurtosis {
    my $individual = shift @_;
    my $population = shift @_;
    my @z_scores = map { $_->{fitness} } @$population;
    my $z_score = $individual->{fitness};

    # Get the distribution f(Z)
    my @distribution;
    foreach my $score (@z_scores)
    {
	my $occurrences = 0;
	foreach my $element (@z_scores) 
	{ 
	    if ($element == $score) 
	    { 
		$occurrences++;
	    }
	}
	push @distribution, $occurrences;
    }
    my $dist_mean = mean(@distribution);
    my $dist_stddev = stddev(@distribution);

    # Calculate kurtosis for the individual
    # Get number of seed words with the same Z-score
    my $occurrences = 0;
    foreach my $element (@z_scores) 
    { 
	if ($element == $z_score) 
	{ 
	    $occurrences++;
	}
    }

    my $kurtosis = 0;
    if ((scalar(@distribution) - 1) * pow($dist_stddev, 4) != 0)
    { 
	$kurtosis = pow(($occurrences - $dist_mean), 4) * scalar(@distribution) / ((scalar(@distribution) - 1) * pow($dist_stddev, 4)) - 3
    } 

    return $kurtosis;
}

# Function to get the thinness coefficient
sub thinness{
    my $individual = shift @_;
    my $population = shift @_;

    my $std_error = 2 * sqrt(6 / scalar (@$population));

    return (dist_kurtosis($individual, $population) + (2 * $std_error)) / (4 * $std_error);
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

# Function to check if the motif id formed only by valid nucleotides
sub is_valid_motif {
    my $motif = shift @_;

    if ($motif =~ /\A[acgt]+\z/i)  {
	return 1;
    }
    else  {
	return 0;
    }
}

# Function to check if the position is not at the end of the sequence so that there would be no space for the motif
sub is_inside_sequence {
    my $position = shift @_;

    my $shift = 0;
    foreach my $seq (@original_sequences)  {
	$shift += length $seq;
	if (($position < $shift) && ($position + $motif_min_width > $shift)) {
	    return 0;
	}
    }
    return 1;
}
