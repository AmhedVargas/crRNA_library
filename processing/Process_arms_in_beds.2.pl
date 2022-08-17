#!/usr/bin/perl

use List::Util qw/sum/;

#####Documentation######
$num_args = $#ARGV + 1;
if ($num_args != 2) {
  print "\nUsage: Process_arms.pl Bed_intersections_with_arms Pattern_files STDOUT\n\n";
  exit;
}
############
##Define hash tables
my %aacode = (
  TTT => "F", TTC => "F", TTA => "L", TTG => "L",
  TCT => "S", TCC => "S", TCA => "S", TCG => "S",
  TAT => "Y", TAC => "Y", TAA => "STOP", TAG => "STOP",
  TGT => "C", TGC => "C", TGA => "STOP", TGG => "W",
  CTT => "L", CTC => "L", CTA => "L", CTG => "L",
  CCT => "P", CCC => "P", CCA => "P", CCG => "P",
  CAT => "H", CAC => "H", CAA => "Q", CAG => "Q",
  CGT => "R", CGC => "R", CGA => "R", CGG => "R",
  ATT => "I", ATC => "I", ATA => "I", ATG => "M",
  ACT => "T", ACC => "T", ACA => "T", ACG => "T",
  AAT => "N", AAC => "N", AAA => "K", AAG => "K",
  AGT => "S", AGC => "S", AGA => "R", AGG => "R",
  GTT => "V", GTC => "V", GTA => "V", GTG => "V",
  GCT => "A", GCC => "A", GCA => "A", GCG => "A",
  GAT => "D", GAC => "D", GAA => "E", GAG => "E",
  GGT => "G", GGC => "G", GGA => "G", GGG => "G",
); # this is the hash table for the amino acids

my %HoA = (
    F => [ "TTT", "TTC" ],
    L => [ "TTA", "TTG", "CTT" , "CTC", "CTA", "CTG"],
    S => [ "TCT", "TCC" , "TCA", "TCG", "AGT", "AGC"],
    Y => [ "TAT", "TAC"],
    STOP => [ "TGA", "TAA", "TAG"],
    C => [ "TGT", "TGC"],
    W => [ "TGG"],
    P => [ "CCT", "CCC" , "CCA","CCG"],
    H => [ "CAT", "CAC"],
    Q => [ "CAA", "CAG"],
    R => [ "CGT", "CGC" , "CGA", "CGG","AGA", "AGG"],
    I => [ "ATT", "ATC" , "ATA"],
    M => [ "ATG"],
    T => [ "ACT", "ACC" , "ACA", "ACG"],
    N => [ "AAT", "AAC"],
    K => [ "AAA", "AAG"],
    V => [ "GTT", "GTC" , "GTA", "GTG"],
    A => [ "GCT", "GCC" , "GCA", "GCG"],
    D => [ "GAT", "GAC"],
    E => [ "GAA", "GAG"],
    G => [ "GGT", "GGC" , "GGA", "GGG"],
);


##Reading and verifying input commands
my $DB=$ARGV[1];
my $file=$ARGV[0];

###If the following files are not present -> stop the program
if (!(-e $file."")){die $file." file not found!\n";}
if (!(-e $DB."")){die $DB." file not found!\n";}

##Define hashes
my @pat;

##Reading and indexing db
open(op,$DB) or die "cannot open DB file\n";
while($line=<op>){
chomp($line);
push(@pat,uc($line));
}
close(op);

##Reading and indexing fragment
open(op,$file) or die "cannot open intersect file\n";
while($line=<op>){
chomp($line);
#printf $line."\n";
@info = split(/\t/,$line);
$left=uc($info[6]);
$right=uc($info[7]);

$force=1;
while($force ==1){
    $force=0;
for my $check(@pat){
    #printf $check."\n";
    $loc=index($left,$check);
    while($loc >= 0){
    @leftori=split(//,$left);
    $dif = ($loc % 3);
    $dif = $dif - 3;
    if(($loc+$dif) < 0){$dif = $dif +3;}
    ##First AA
    $a=substr($left,($loc+$dif),3);
    ##Second AA
    $b=substr($left,($loc+$dif)+3,3);
    #$c=substr($left,($loc+$dif)+6,3);
    #print "Fourth codon\n".$d."\n";
    $A=alternate_aa($a);
    #print "First mod codon\n".$A."\n";
    $B=alternate_aa($b);
    #$C=alternate_aa($c);
    #print "Second mod codon\n".$B."\n";
    ($leftori[($loc+$dif)],$leftori[($loc+$dif)+1],$leftori[($loc+$dif)+2])=split(//,$A);
    ($leftori[($loc+$dif)+3],$leftori[($loc+$dif+3)+1],$leftori[($loc+$dif+3)+2])=split(//,$B);
    #($leftori[($loc+$dif)+6],$leftori[($loc+$dif+6)+1],$leftori[($loc+$dif+6)+2])=split(//,$C);
    $left=join('',@leftori);
    $loc=index($left,$check);
    $force=1;
    }

    $loc=index($right,$check);
    while($loc >= 0){
    @rightori=split(//,$right);
    $dif = ($loc % 3);
    $dif= $dif - 3;
    if(($loc+$dif) < 0){$dif = $dif +3;}
    ##First AA
    $a=substr($right,($loc+$dif),3);
    ##Second AA
    $b=substr($right,($loc+$dif)+3,3);
    #$c=substr($right,($loc+$dif)+6,3);
    #print "Fourth codon\n".$d."\n";
    $A=alternate_aa($a);
    #print "First mod codon\n".$A."\n";
    $B=alternate_aa($b);
    #$C=alternate_aa($c);
    #print "Second mod codon\n".$B."\n";
    ($rightori[($loc+$dif)],$rightori[($loc+$dif)+1],$rightori[($loc+$dif)+2])=split(//,$A);
    ($rightori[($loc+$dif)+3],$rightori[($loc+$dif+3)+1],$rightori[($loc+$dif+3)+2])=split(//,$B);
    #($rightori[($loc+$dif)+6],$rightori[($loc+$dif+6)+1],$rightori[($loc+$dif+6)+2])=split(//,$C);
    $right=join('',@rightori);
    $loc=index($right,$check);
    $force=1;
    }
}
}
printf $line."\t".$left."\t".$right."\n";

}

close(op);


exit;

sub revcomp{
$Revcomp = reverse $_[0];
$Revcomp =~ tr/ACGTacgt/TGCAtgca/;
return($Revcomp)
}

sub match_positions {
    my ($regex, $string) = @_;
    return if not $string =~ /$regex/;
    return ($-[0], $+[0]);
}
sub match_all_positions {
    my ($regex, $string) = @_;
    my @ret;
    while ($string =~ /$regex/g) {
        push @ret, [ $-[0], $+[0] ];
    }
    return @ret
}

sub alternate_aa{
$pat=$_[0];
$aa=$aacode{$pat};
#print "Size:".$#{ $HoA{$aa} }."\n";
if($#{ $HoA{$aa} } < 1){return($pat)}
else{
    $reto=rand($#{ $HoA{$aa} });
if($HoA{$aa}[$reto]==$pat){
    $reto=rand($#{ $HoA{$aa} });
    return($HoA{$aa}[$reto])
    }else{return($HoA{$aa}[$reto])}
}
}
