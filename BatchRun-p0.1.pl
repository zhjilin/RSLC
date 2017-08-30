#!/usr/bin/perl -w

=head1 Usage

    This is a wrapper based on the UMI counting script.Thus, the GuideUMI-p0.1.pl must be contained in the same folder. 
	It contains two functions:
    1) Parallel counting for a group of fastq files
    2) Merging the UMI tables

    usage:

    --fastq     str, input fastq files separate by comma, gzip supported
    --library   str, the table that contains guides
    --n         int, throw UMI with reads count lower than this number (default: no filtering)
    --nmerge    int, UMI with read counts <= this number(default 3) will be allowed to group into the UMI 1-hamming dist
    --out       str, output prefix for the combined table, default "output"
    --mismatch  int
    --step      str, processing steps, default 12;
	--group		specify to use group function
   	--help		see the help

=cut 


use strict;
use Data::Dumper;
use File::Basename qw/dirname basename/;
use Getopt::Long;
use FindBin qw($Bin);
my($Prefix,$List, $Alien, $Mismatch,$Shift, $Const, $Num,$Xlen,$Lnum,$step,$Group, $Help);

GetOptions(
    'fastq:s' => \$List,
    'library:s' => \$Alien,
    'n:i' => \$Num,
    'mismatch:n' =>\$Mismatch,
    'read_length:i' => \$Xlen,
    'nmerge:i' => \$Lnum,
    'out:s' =>\$Prefix,
    'group!'    =>\$Group,
    'step:s' => \$step,
    'help!' => \$Help
);
#$Xlen ||= 55;
$Num ||= -1;
$Lnum ||=3;
$Mismatch ||=0;
$step ||="12";
$Prefix||="output";


die `pod2text $0` if(not defined $Alien);
die `pod2text $0` if(defined $Help);
my @file;
if($List=~/\,/){
    @file=split(/,/ ,$List);
}
else{
    push @file,$List;
}

my $program="GuideUMI-p0.1.pl";
my $para="--library $Alien --nmerge $Lnum ";
$para.="--n $Num " if($Num >0);
$para.="--mismatch $Mismatch " if(defined $Mismatch && $Mismatch >=1);
$para.="--group" if(defined $Group );
#partition the task    
if($step =~/1/){
my @pids;
print STDERR "@file are being processed, cross your fingers...\n";
for my $subf(@file){
    my $pid=fork();
    if($pid== -1){
        die;
    } elsif( $pid ==0){
        print STDERR "start counting file $subf \n";
        exec "perl $Bin/$program $para $subf" or die"failed to run \n$!";
        
    }
    push @pids, $pid;
}
while (wait() != -1) {}
print STDERR "Done for counting UMIs\n";
}
#Merging tables resulted from the first step;
my @tf1=@file;
my (@r1, @r2, @outtag);

@r1= map {$_=~s/.*\///g;$_=~s/_.*fastq(.gz)?$|_.*fq(.gz)?$/\.count/;$_ } @file;
@r2= map {$_=~s/.*\///g;$_=~s/_.*fastq(.gz)?$|_.*fq(.gz)?$/\.count\.group/;$_} @tf1;

print STDERR "@r1\n";
print STDERR "@r2\n"  if(defined $Group );
if($step=~/2/){
    my %alien;
    &ReadTab($Alien,\%alien);
    print STDERR "Starting merging tables...\n";
    &Combine(\@r1,"$Prefix.raw", \%alien);
    &Combine(\@r2,"$Prefix.group", \%alien) if(defined $Group );
}
print STDERR "Done\n";

sub ReadTab{
    my $f=shift;
    my $h=shift;
    open RT,$f;
    while(<RT>){
        chomp;
        s/\r//g;
        my @t=split ",";
        $h->{$t[0]}=[$t[1],$t[2]];
  }
  close RT;
}


sub Combine{
    my $h=shift;
    my $f=shift;
    my $al=shift;
    open OI, ">summary_count.$f";
    my (%uniq,%q);
## first get a non-redundant hash table from all the tables (key,value)=(ID_UMI)    
    foreach my $fn (@$h){
        open FN, $fn;
        while(<FN>){
            chomp;
            my($s,$n)=($1,$2) if(/^([^,]+),([^,]+)/);
            my ($s1,$s2)=($1,$2) if($s =~ /^(\S+)_([^_]+)$/);
            $uniq{$s}=$n if(not exists $uniq{$s}); #to store the list of guide with the unique UMI found in the counting files
            $q{$s1}=$s1 if(not exists $q{$s1});  #to store the list of unique guide found in the counting files
        }
        close FN;
    }
#    print Dumper(%uniq);

     foreach my $fn (@$h){
         open FN, $fn;
         my %tq; #a temporary hash to store the current founded guide with a set of UMI
         while(<FN>){
             chomp;
             my @t=split /,/;
             if(exists $uniq{$t[0]}){
                 $uniq{$t[0]}.=",$t[2]";
                 $tq{$t[0]}=$t[1];
             }
        }
         close FN;
#the following is to make 0 for missing umi for each time course
         for my $k1 (keys  %uniq){
             next if(exists $tq{$k1});
             $uniq{$k1}.=",0";
        }
     }

#final step to compensate the missing guide
     my $sfnum=@$h;
     for my $ak (keys %$al){
         next if(exists $q{$ak});
         $uniq{"$ak\_NNNNNN"}= $al->{$ak}[1].(",0" x $sfnum);
     }
     my @newname=map {$_ =~s/\.count.*$// ; $_ } @$h;
     my @cname=@newname;
     shift @cname; 
     print OI "RSL.guide,guide.set,Control,".join(",",@cname)."\n";
     for my $k (sort {$a cmp $b }keys %uniq){
         my ($k1,$k2)=($1,$2) if($k=~/^(\S+)_([^_]+)$/);
         print OI "$k1\_$k2,$uniq{$k}\n";
     }
     close OI;
}
