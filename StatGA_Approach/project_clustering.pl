
#!/usr/bin/perl -w
use Bio::Perl;

use strict;
use warnings;
use Data::Dumper;
use List::Util qw(shuffle max min sum);
use Statistics::Basic qw(:all);
use String::Approx;
use PDL::Basic;
use Getopt::Long;
use List::Compare;
use Switch;
use Benchmark qw(:all);

# Measuring computation time
my $t0 = Benchmark->new;

# Check the validity of the parameters
my %args;
GetOptions(\%args,
	   "output_instances=s",
	   "stats_file=s",
	   "input_sequences=s",
           "max_similarity=s",
) or die "Invalid arguments!";
die "Missing -output_instances!" unless $args{output_instances};
die "Missing -stats_file!" unless $args{stats_file};
die "Missing -input_sequences!" unless $args{input_sequences};
die "Missing -max_similarity!" unless $args{max_similarity};

my $max_similarity = $args{max_similarity};
my $unpaired_weight = 0.6;

# Open the input file
my $file = "inner_output.txt";
open (INPUTFILE, $file) or die "Could not open temporary results file";

# Getting the input file into an string an splitting it sequence by sequence
my @input = <INPUTFILE>;
my $input_string = join("", @input);
$input_string = substr($input_string, 1);
my $output_instances = $args{output_instances};
my @solutions_raw = split("\n", $input_string);

my @solutions;
foreach my $solution (@solutions_raw)
{
    my @info = split(" ", $solution);
    # Getting position, consensus sequence and fitness
    my $position = $info[0];
    my $width = $info[1];
    my $consensus = $info[2];
    my $fitness = $info[3];

    if ($fitness < 0.01)  {
	push @solutions, { position => $position, width => $width, consensus => $consensus, fitness => $fitness };
    }
}

# Close the input file
close (INPUTFILE);

# Open the input sequences (FASTA FORMAT)
my $input_seq = $args{input_sequences};
open (INPUTSEQ, $input_seq) or die "Could not open $args{input_sequences}. Please, use a valid file name";

# Getting the input sequences into an string an splitting it sequence by sequence
my @inputseqs = <INPUTSEQ>;
my $input_seq_string = join("", @inputseqs);
$input_seq_string = substr($input_seq_string, 1);
my @fasta_sequences = split("\n>", $input_seq_string);
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
    $total_sequence_length += length($seq);
    $super_seq .= $seq;
    push @original_sequences, $seq;
}

# Close the input file
close (INPUTSEQ);

my @sorted_solutions = sort { $a->{fitness} <=> $b->{fitness} || $b->{width} <=> $a->{width}} @solutions;

my @final_solutions;

foreach my $sol (@sorted_solutions)
{
    my $unique = 1;
    foreach my $final_sol (@final_solutions) {
	my ($distance, $rel_pos) = motif_distance($sol->{consensus}, $final_sol->{consensus});
	if ($distance <= $max_similarity)  {
	    $unique = 0;
	    last;
	}
    }
    if ($unique == 1) {
	push @final_solutions, $sol;
    }
}
# Printing the results in output_file (matrix) and in output_instances_file (instances)
open (OUTPUT_INSTANCES, ">>$output_instances") or die "Could not open $args{output_instances}. Please, use a valid file name";

foreach my $sol (@final_solutions)
{
    my @motif_list;
    my @instance_list;
    my $position = $sol->{position};
    my $fitness = $sol->{fitness};
    my $motif_width = length($sol->{consensus});
    my $adjusted_similarity= $max_similarity;
    my $first_motif = find_motif($sol->{position}, $motif_width);
    push @motif_list, $first_motif;
    my ($first_seq_number, $first_seq_pos) = seq_position($sol->{position});
    push @instance_list, { seq_number => $first_seq_number, seq_pos => $first_seq_pos, abs_pos => $sol->{position}, instance => $first_motif };
    foreach my $new_sol (@sorted_solutions) {
	my ($distance, $rel_pos) = motif_distance($sol->{consensus}, $new_sol->{consensus});
	my $overlap = 0;
	foreach my $inst (@instance_list) {
	    if (abs(($new_sol->{position} + $rel_pos) - $inst->{abs_pos}) < $motif_width)  {
		$overlap++;
	    }
	}
	if (($overlap == 0) && ($distance <= $adjusted_similarity) && (is_inside_sequence($new_sol->{position} + $rel_pos, $motif_width)))  {
	    my $motif = find_motif($new_sol->{position} + $rel_pos, $motif_width);
	    if (length($motif) == $motif_width)   {
		push @motif_list, $motif;
		my ($seq_number, $seq_pos) = seq_position($new_sol->{position} + $rel_pos);
		push @instance_list, { seq_number => $seq_number, seq_pos => $seq_pos, abs_pos => $new_sol->{position} + $rel_pos, instance => $motif };
	    } 
	} 
    }


    # Calculate the fitness of a CTM for the given solution

    # Vocabulary file for CTM
    my @words;
    foreach my $instance (@instance_list)
    {
        push @words, $instance->{instance};
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

    # Documents file for CTM
    open (OUTPUT_CTM_SEQS, ">./ctm_seqs.txt") or die "Could not open ctm_seqs.txt.";

	
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
    my $perplexity = `Rscript ./getPerplexity.R`;

    my $perp_str = substr($perplexity, 4, length($perplexity) - 4);

    my $perp = sprintf("%.6f", $perp_str);

    # Print only if fitness < 100
    if ($perp < 100) {

	# Print instances
	my @sorted_instances = sort { $a->{seq_number} <=> $b->{seq_number} || $a->{seq_pos} <=> $b->{seq_pos}} @instance_list;
	printf OUTPUT_INSTANCES "\n\n\n>Position %d, Fitness %f, Perplexity %f\n\n", $position, $fitness, $perp;
	printf OUTPUT_INSTANCES ">instances\n";
	foreach my $instance (@sorted_instances) {
	    printf OUTPUT_INSTANCES "%d,%d,%s\n", $instance->{seq_number}, $instance->{seq_pos}, $instance->{instance};
	}
    }
}

close(OUTPUT_INSTANCES);

my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);
my $stats = $args{stats_file};
open (STATS, ">>$stats") or die "Could not open $args{stats_file}. Please, use a valid file name";
print STATS "Elapsed time:",timestr($td),"\n";
close(STATS);

# Delete temporary files
unlink "inner_output.txt", "ctm_seqs.txt", "ctm_vocab.txt";


# Function to measure the distance between two motifs
sub motif_distance {
    my $motif1 = shift @_;
    my $motif2 = shift @_;

    my $offset;

    my $length1 = length $motif1;
    my $length2 = length $motif2;
    my $min_length = min ($length1, $length2);

    my $submotif1;
    my $submotif2;
    my @distances;
    foreach my $pos (1..$length1)  {
	if ($pos <= $min_length)  {
	    $submotif1 = substr($motif1, 0, $pos);
	    $submotif2 = substr($motif2, $length2 - $pos, $pos);
	} else  {
	    $submotif1 = substr($motif1, $pos - $min_length, $min_length);
	    $submotif2 = substr($motif2, 0, $min_length);
	}
	my $total = 0;
	foreach my $new_pos (0..min($pos, $min_length) - 1) {
	    my $indiv_distance =  symbol_distance(substr($submotif1, $new_pos, 1), substr($submotif2, $new_pos, 1));
	    $total += $indiv_distance;
	}
	my $unpaired = max($length1, $length2) - $pos;
	if ($length1 < $length2) {
	    if ($pos <= $min_length)  {
		$offset = $length2 - $pos;
	    } else  {
		$offset = $pos - $length2;
	    }
	} else  {
	    $offset = $length2 - $pos;
	}
	push  @distances, { position => $offset, distance => ($total + ($unpaired * $unpaired_weight)) };
    }

    foreach my $pos (1..$length2)  {
	if ($pos <= $min_length)  {
	    $submotif1 = substr($motif2, 0, $pos);
	    $submotif2 = substr($motif1, $length1 - $pos, $pos);
	} else  {
	    $submotif1 = substr($motif2, $pos - $min_length, $min_length);
	    $submotif2 = substr($motif1, 0, $min_length);
	}
	my $total = 0;
	foreach my $new_pos (0..min($pos, $min_length) - 1) {
	    my $indiv_distance =  symbol_distance(substr($submotif1, $new_pos, 1), substr($submotif2, $new_pos, 1));
	    $total += $indiv_distance;
	}
	my $unpaired = max($length1, $length2) - min($pos, $min_length);
	if ($length1 < $length2) {
		$offset = $pos - $length1;
	} else  {
	    if ($pos <= $min_length)  {
		$offset = $pos - $length1;
	    } else  {
		$offset = $length1 - $pos;
	    }
	}   
	push  @distances, { position => $offset, distance => ($total + ($unpaired * $unpaired_weight))};
    }

    my @sorted_distances = sort { $a->{distance} <=> $b->{distance} || $a->{position} <=> $b->{position}} @distances;
    return ($sorted_distances[0]->{distance} / mean(length($motif1), length($motif2))), $sorted_distances[0]->{position};
}

# Function to measure the distance between two IUPAC symbols
sub symbol_distance {
    my $symbol1 = shift @_;
    my $symbol2 = shift @_;

    my @list1 = iupac_to_list($symbol1);
    my @list2 = iupac_to_list($symbol2);

    my $lc = List::Compare->new(\@list1, \@list2);

    my $distance = 1 - 2* ($lc->get_intersection / (scalar(@list1) + scalar(@list2)));

    return $distance;
}

# Function to transform a IUPAC symbol into a list of possible values
sub iupac_to_list {
    my $symbol = shift @_;

    my @list;
    switch ($symbol) {
	case /A/i { @list = ("A"); }
	case /C/i { @list = ("C"); }
	case /G/i { @list = ("G"); }
	case /T/i { @list = ("T"); }
	case /R/i { @list = ("A", "G"); }
	case /Y/i { @list = ("C", "T"); }
	case /S/i { @list = ("G", "C"); }
	case /W/i { @list = ("A", "T"); }
	case /K/i { @list = ("G", "T"); }
        case /M/i { @list = ("A", "C"); }
	case /B/i { @list = ("C", "G", "T"); }
	case /D/i { @list = ("A", "G", "T"); }		
	case /H/i { @list = ("A", "C", "T"); }
	case /V/i { @list = ("A", "C", "G"); }
	case /N/i { @list = ("A", "C", "G", "T"); }
	else { @list = ("A", "C", "G", "T"); }	
    }

    return @list;
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

# Function to check if the motif is formed only by valid nucleotides
sub is_valid_motif {
    my $motif = shift @_;

    if ($motif =~ /\A[acgt]+\z/i)  {
	return 1;
    }
    else  {
	return 0;
    }
}

# Function that returns the motif given the absolute position
sub find_motif {
    my $pos = shift @_;
    my $width = shift @_;
    my $motif;

    $motif = substr ($super_seq, $pos, $width);

    return $motif;
}


# Function to find the number of sequence and the position within it
sub seq_position {
    my $position = shift @_;
    my $seq_number = 0;
    my $shift = 0;
    foreach my $seq (@original_sequences)  {
	$shift += length $seq;
	if ($position < $shift) {
	    return $seq_number, $position - $shift;
	}
	$seq_number++;
    }
}


# Function to check if the position is not at the end of the sequence so that there would be no space for the motif
sub is_inside_sequence {
    my $position = shift @_;
    my $width = shift @_;

    my $shift = 0;
    foreach my $seq (@original_sequences)  {
	$shift += length $seq;
	if (($position < $shift) && ($position + $width > $shift)) {
	    return 0;
	}
    }
    return 1;
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

# Function to round up a number
sub roundup {
    my $n = shift;
    return(($n == int($n)) ? $n : int($n + 1));
}

# Function to get the hamming distance between two motifs
sub hammingDistance
{ 

    return length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] );

}
