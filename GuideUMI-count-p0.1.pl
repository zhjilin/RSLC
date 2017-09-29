#!/usr/bin/perl -w
=head1 USAGE


    perl GuideUMI-count-p0.1.pl --library TF_guides_file fastq
    
    options:
    --library   str, comma separated table containing the information of TF guides, shuold be "identifier,guide,gene_name"
    --n         int, cutoff to remove low reads count umi, default (no filtering)
    --nmerge    int, counts number cutoff for merging within 1 hamming distance
    --mismatch  int, allow mismatches in the guide. This function will use bowtie and slow down the process, not recommended. 
    --dual      str, dual index mode, XR1+6index1+6index2+0R2
	--group		specify to use the 1-hamming distance grouping function
    --help      to print the help information on the screen


=cut 


use strict;
use Data::Dumper;
use File::Basename qw/dirname basename/;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Getopt::Long;

my ($Num,$Alien,$Mismatch,$lnum,$Climit,$xLen,$Group,$Help);

GetOptions(
    'n:i'=>\$Num,
    'library:s'=>\$Alien,
    'mismatch:n' => \$Mismatch,
    'nmerge:i'=>\$lnum,
    'group!' => \$Group,
    'help!'=>\$Help
);
$Mismatch ||=0;
$Num ||= -1;  
$lnum ||=3;   
die `pod2text $0` if( not defined $Alien); 

my $fastq =shift;
my $fbase=basename $fastq;
$fbase=~s/.fq(.gz)?$|.fastq(.gz)?$//;
my $prefix_tag= $1 if ($fbase=~/^([^_]+)_.*/);



##Check if bowtie has been setup in PATH

my $status=`bowtie --version |grep Compiler`;
die "bowite was not available in \$PATH $!\n" if(not defined $status);

my (%bigtab,%alien,%cpool);


if(defined $Alien){
    if( $Mismatch >=1 ){
    `cut -d , -f2 $Alien |awk '{print ">"\$1"\\n"\$1}' > TF_guide_seq.fa` if(not -e "TF_guide_seq.fa");
    system("bowtie-build TF_guide_seq.fa TF_guide_seq_ref") if(not -e "TF_guide_seq_ref.1.ebwt");
    }
    &ReadTab($Alien,\%alien);
}
open OUT2,">$prefix_tag.count";
open OUTU1,">$prefix_tag.UMI";
if(defined $Group){
open OUTUP,">$prefix_tag.grp.UMI";
open OUT3,">$prefix_tag.count.group" ;
}
#my($ocount,$gocount,$utable,$gutable);

my $c=0;
## The function modified as the following description
# condition 1) the constant part starts from the exact position 26(0 based)
# condition 2) the constant sequence shifted to the right, it means stating from 26+
my $instream;
if($fastq=~/\.gz$/){
    $instream=IO::Uncompress::Gunzip->new($fastq,MultiStream =>1);
}
else{
    open $instream,'<', $fastq;
}

while(<$instream>){
    last if ($instream->eof());
        chomp;
        my $line1=$_;
        my $line2=<$instream>;
        my $line3=<$instream>;
        my $line4=<$instream>;
            chomp $line1;
            chomp $line2;
            my($umi,$ligand);
            $ligand=substr($line2,0,20);
            my $processed_ligand;
            if( $Mismatch >=1){
                $processed_ligand=&AllowMismatch($ligand, $Mismatch);
                next if(not defined $processed_ligand);
            }
            else{
                $processed_ligand=$ligand;
            }
            if(exists $alien{$processed_ligand} ){
                $umi=substr($line1,length($line1)-13,6);
                if($umi =~/N/g){next;}
                $bigtab{$processed_ligand}{$umi}+=1;
            }
}
    close $instream;


my (@stat,@rstat, %etest);
for my $k (keys %bigtab){
    next if(not exists $alien{$k});
    $etest{$k}=$alien{$k}[0];
    my @punch;
    my ($count,$cf)=(0,0);
    for my $uk( sort {$bigtab{$k}{$b} <=> $bigtab{$k}{$a}} keys %{$bigtab{$k}}){
        $count+=$bigtab{$k}{$uk};
        next if( $bigtab{$k}{$uk} <= $Num);
        print OUT2 "$alien{$k}[0]_$uk,$alien{$k}[0],$bigtab{$k}{$uk}\n";
        push @punch,$uk;
        $stat[$bigtab{$k}{$uk}]+=1;
        $cf++;
    }
    print OUTU1 "$k,$alien{$k}[0],$cf\n";
    if(defined $Group){
    my %new=&GroupUMI(\@punch,\%{$bigtab{$k}});
    my ($count2,$cf2 )=(0,0);
    for my $uk (sort {$new{$b} <=> $new{$a}} keys %new){
        $count2 += $new{$uk};
        next if( $new{$uk} <=$Num);
        print OUT3 "$alien{$k}[0]_$uk,$alien{$k}[0],$new{$uk}\n";
        $rstat[$new{$uk}]+=1;
        $cf2++;
    }
    print OUTUP "$k,$alien{$k}[0],$cf2\n";
    }
}

close OUT2;
close OUTU1;

if(defined $Group){
close OUT3;
close OUTUP;
}





sub ReadTab{
    my $f=shift;
    my $h=shift;
    open RT,$f;
    while(<RT>){
        chomp;
        s/\r//g;
        my @t=split ",";
        $h->{$t[1]}=[$t[0],$t[2]];
    }
    close RT;
}



sub GroupUMI{
    my $a=shift;
    my $h=shift;

    my %tem=%$h;
    for(my $i=0; $i <= $#$a ; $i++){
        for(my $j= $#$a; $j> $i; $j--){
            my $hd=&HammingDist($a->[$i], $a->[$j]);
            if($hd == 1 && (exists $tem{$a->[$j]})){
                next if $tem{$a->[$j]} > $lnum;
                $tem{$a->[$i]}+=$tem{$a->[$j]};
                delete $tem{$a->[$j]};
            }
        }
    }
    return %tem;
}

sub HammingDist{
    my $s1=shift;
    my $s2=shift;
    my $td=0;
    for(my $i=0; $i <length($s1); $i++){
        my $s1e=substr($s1,$i,1);
        my $s2e=substr($s2,$i,1);
        if($s1e ne $s2e){
            $td+=1;
        }
    }
    return $td;
}

sub AllowMismatch{
    my $query=shift;
    my $mis=shift;
    my $bowtieout=`bowtie -n $mis -l 20 -c TF_guide_seq_ref $query 2>&1 |grep -v ^# `;
    if($bowtieout !~/^No/){
    my @field=split /\s+/, $bowtieout;
    return $field[2] if(defined $field[2]);
    }
    else {
        return();
        }
}
