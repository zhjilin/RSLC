#!/usr/bin/perl

# count RSLs using only first given number of bases i.e sum together those that identical bases
# in the rest of the positions  
# T. Kivioja, 2016

use strict;
use warnings;


if (scalar(@ARGV) < 3) {
    print STDERR "Usage: perl countTruncatedRSLs.pl trunclen mincount input_countfile.csv\n";
    exit(1);
}

my $tlen = $ARGV[0];       # truncated length
my $min_count = $ARGV[1];  # min count  
my $input_file = $ARGV[2]; # input file containing all counts


open(IFH, "<$input_file") or die "Error: could not open input file $input_file: $!"; 

# put summary counts to memory - should not be two many if only 16 or 32 per guide
my %trunc_counts = ();

my $header_line = <IFH>;
my $lines = 0;
while (my $line = <IFH>) {
    chomp $line;

    my ($rsl_guide, $guide_set, $control_count, $treatment_count) = split /,/, $line;
    # remove the guide set from the rsl guide, leave only the label
    $rsl_guide =~ s/^[\w-]+_\d+_//;
    
    # take the first bases
    my $trsl = substr($rsl_guide, 0, $tlen);  
    if ($control_count >= $min_count || $treatment_count >= $min_count) {
	if (!defined($trunc_counts{$guide_set}->{$trsl})) {
	    $trunc_counts{$guide_set}->{$trsl}->[0] = $control_count; 
	    $trunc_counts{$guide_set}->{$trsl}->[1] = $treatment_count;
	} else {
	    $trunc_counts{$guide_set}->{$trsl}->[0] += $control_count;
	    $trunc_counts{$guide_set}->{$trsl}->[1] += $treatment_count;
	}
    }
    $lines++;
    print STDERR "Processed $lines lines\n" if ($lines % 1e6 == 0);
}


# Print the counts of the truncated labels
print $header_line;
for my $guide_set (sort keys %trunc_counts) {

    for my $trsl (sort keys %{ $trunc_counts{$guide_set} }) {
	my $control_count = $trunc_counts{$guide_set}->{$trsl}->[0];
	my $treatment_count = $trunc_counts{$guide_set}->{$trsl}->[1];
	
	print $guide_set."_".$trsl, ",", $guide_set, ",", $control_count, ",", $treatment_count, "\n";
	
    }

}






