#!/usr/bin/perl

use List::Util qw/sum/;

#####Documentation######
$num_args = $#ARGV + 1;
if ($num_args != 4) {
  print "\nUsage: Process_arms4ATG.pl Genebodies_data.hom.tsv Bed_intersect Column_for_crRNA Type\n\n";
  print "Type: being CDS, ATG, or STOP\n";
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
my $DB=$ARGV[0];
my $file=$ARGV[1];
my $column=$ARGV[2];
my $Tipin=$ARGV[3];

###If the following files are not present -> stop the program
if (!(-e $file."")){die $file." file not found!\n";}
if (!(-e $DB."")){die $DB." file not found!\n";}

##Define hashes
my %hash;
my $key;

##Reading and indexing db
open(op,$DB) or die "cannot open DB file\n";
while($line=<op>){
chomp($line);
@info = split(/\t/,$line);
@datata=split(/;/,$info[3]);
$key=$datata[0].";".$datata[1];
$hash{$key."strand"}=$info[5];
$hash{$key."ori"}=$info[8];
$hash{$key."frame"}=$info[9];
$hash{$key."exon"}=$info[10];
$hash{$key."seq"}=$info[11];
}
close(op);

##Reading and indexing fragment
open(op,$file) or die "cannot open intersect file\n";
while($line=<op>){
chomp($line);
@info = split(/\t/,$line);
@datata=split(/;/,$info[3]);
$key=$datata[0].";".$datata[1];
$flag=0;
if($hash{$key."strand"} ne $info[11]){$flag=1;}

@data=split(/;/,$info[9]);

$crRNA=$data[$column-1];
if($flag==1){$crRNA=revcomp($crRNA);}

($start,$end)=match_positions($crRNA,$hash{$key."seq"});

@seq=split(//,$hash{$key."seq"});
@ori=split(//,$hash{$key."ori"});
@frame=split(//,$hash{$key."frame"});
@exon=split(//,$hash{$key."exon"});

#Make sure first and last 3bp are not considered as exons
$ats=index($hash{$key."exon"},"1");
##Now index inverse
$RevEx=reverse($hash{$key."exon"});
$tas=(scalar(@exon)-index($RevEx,"1")-1);

#print join('',@ori[($ats)..($ats+2)])."\t".join('',@ori[($tas-2)..($tas)])."\n";

if($Tipin eq 'ATG'){
$cutsite=$ats +3;
}
if($Tipin eq 'STOP'){
$cutsite=$tas-2;
}

if($Tipin eq 'CDS'){
$exon[$ats]=2;
$exon[$ats+1]=2;
$exon[$ats+2]=2;
$exon[$tas-2]=2;
$exon[$tas-1]=2;
$exon[$tas]=2;
##Cutting is 3bp before the end
$cutsite=($end-3);
##Unless sequence was in the reverse complement
if($flag==1){$cutsite=($start+3);}
##If in intron and 1 before exon, go before
if(($exon[$cutsite]==0)&($exon[$cutsite-1]==1)){$cutsite=$cutsite-1;}
##If in intron and 1 after exon, go after
if(($exon[$cutsite]==0)&($exon[$cutsite+1]==1)){$cutsite=$cutsite+1;}
##If in intron and 2 before exon, go before
if(($exon[$cutsite]==0)&($exon[$cutsite-2]==1)){$cutsite=$cutsite-2;}
##If in intron and 2 after exon, go after
if(($exon[$cutsite]==0)&($exon[$cutsite+2]==1)){$cutsite=$cutsite+2;}
##If in intron and 3 before exon, go before
if(($exon[$cutsite]==0)&($exon[$cutsite-3]==1)){$cutsite=$cutsite-3;}
##If in intron and 3 after exon, go after
if(($exon[$cutsite]==0)&($exon[$cutsite+3]==1)){$cutsite=$cutsite+3;}
##When looking for ATG or STOP
##If in intron and 6 before exon, go before
if(($exon[$cutsite]==0)&($exon[$cutsite-6]==1)){$cutsite=$cutsite-6;}
##If in intron and 6 after exon, go after
if(($exon[$cutsite]==0)&($exon[$cutsite+6]==1)){$cutsite=$cutsite+6;}
##If in intron and 9 before exon, go before
if(($exon[$cutsite]==0)&($exon[$cutsite-9]==1)){$cutsite=$cutsite-9;}
##If in intron and 9 after exon, go after
if(($exon[$cutsite]==0)&($exon[$cutsite+9]==1)){$cutsite=$cutsite+9;}
##If in intron and 9 before exon, go before
if(($exon[$cutsite]==0)&($exon[$cutsite-10]==1)){$cutsite=$cutsite-10;}
##If in intron and 9 after exon, go after
if(($exon[$cutsite]==0)&($exon[$cutsite+10]==1)){$cutsite=$cutsite+10;}
#If any, stay there
##If frame 0 remain the same
#If frame 2 add 1
if($frame[$cutsite] == 2){$cutsite++;}
#If frame 1 remove 1
if($frame[$cutsite] == 1){$cutsite--;}
}
##Now check if spacer sequence is on homology arm
#return original sequence

#printf $crRNA."\t".$hash{$key."strand"}."\t".$info[11]."\t".$start."\t".$end."\t".join('',@seq[$start..($end-1)])."\t".$cutsite."\t".$seq[$cutsite]."\t".$frame[$cutsite]."\t".join('',@ori[$start..($end-1)])."\t".join('',@frame[$start..($end-1)])."\t".join('',@exon[$start..($end-1)])."\t";
#printf join('',@ori[($cutsite-45)..($cutsite-1)])."\t".join('',@ori[($cutsite)..($cutsite+44)])."\n";
$LeftOri=join('',@ori[($cutsite-45)..($cutsite-1)]);
$RightOri=join('',@ori[($cutsite)..($cutsite+44)]);

$LeftSeq=join('',@seq[($cutsite-45)..($cutsite-1)]);
$RightSeq=join('',@seq[($cutsite)..($cutsite+44)]);

$LeftFrame=join('',@frame[($cutsite-45)..($cutsite-1)]);
$RightFrame=join('',@frame[($cutsite)..($cutsite+44)]);

$LeftExon=join('',@exon[($cutsite-45)..($cutsite-1)]);
$RightExon=join('',@exon[($cutsite)..($cutsite+44)]);

##Check if spacer in left Hom arm
$loc=index($LeftSeq,$crRNA);
if($loc > 0){
    @leftframe=split(//,$LeftFrame);
    @leftori=split(//,$LeftOri);
    @leftseq=split(//,$LeftSeq);
    @leftexon=split(//,$LeftExon);
    #if inside frame
    $inframe= sum(@leftframe[($loc)..($loc+19)]);
    $dif=0;
    if($inframe > 1){
        $dif=index(join('',@leftframe[($loc)..($loc+19)]),"2");
        $dif = ($dif % 3);
    }
    ##First AA
    #$a=substr($LeftSeq,($loc+$dif),3);
    ##Second AA
    #$b=substr($LeftSeq,($loc+$dif)+3,3);
    ##Third AA
    #$c=substr($LeftSeq,($loc+$dif)+6,3);
    ##Fourth AA
    $d=substr($LeftSeq,($loc+$dif)+9,3);
    ##Fifth AA
    #$e=substr($LeftSeq,($loc+$dif)+18,3);
    #print "Fourth codon\n".$d."\n";
    #$A=alternate_aa($a);
    #print "First mod codon\n".$A."\n";
    #$B=alternate_aa($b);
    #print "Second mod codon\n".$B."\n";
    #$C=alternate_aa($c);
    #print "Third mod codon\n".$C."\n";
    $D=alternate_aa($d);
    #print "Third mod codon\n".$C."\n";
    #$E=alternate_aa($e);
    #($leftori[($loc+$dif)],$leftori[($loc+$dif)+1],$leftori[($loc+$dif)+2])=split(//,$A);
    #($leftori[($loc+$dif)+3],$leftori[($loc+$dif+3)+1],$leftori[($loc+$dif+3)+2])=split(//,$B);
    #($leftori[($loc+$dif)+6],$leftori[($loc+$dif+6)+1],$leftori[($loc+$dif+6)+2])=split(//,$C);
    ($leftori[($loc+$dif)+9],$leftori[($loc+$dif+9)+1],$leftori[($loc+$dif+9)+2])=split(//,$D);
    #($leftori[($loc+$dif)+18],$leftori[($loc+$dif+18)+1],$leftori[($loc+$dif+18)+2])=split(//,$E);
    $LeftOri=join('',@leftori);
}

##Check if spacer in left Hom arm
$loc=index($RightSeq,$crRNA);
if($loc > 0){
    @rightframe=split(//,$RightFrame);
    @rightori=split(//,$RightOri);
    @rightseq=split(//,$RightSeq);
    @rightexon=split(//,$RightExon);
    #if inside frame
    $inframe=sum(@rightframe[($loc)..($loc+19)]);
    $dif=0;
    if($inframe > 1){
        $dif=index(join('',@rightframe[($loc)..($loc+19)]),"2");
        $dif = ($dif % 3);
    }
    ##First AA
    #$a=substr($RightSeq,($loc+$dif),3);
    ##Second AA
    #$b=substr($RightSeq,($loc+$dif)+3,3);
    ##Third AA
    #$c=substr($RightSeq,($loc+$dif)+6,3);
    ##Fourth AA
    $d=substr($RightSeq,($loc+$dif)+9,3);
    ##Fifth AA
    #$e=substr($RightSeq,($loc+$dif)+18,3);
    #print "Fourth codon\n".$d."\n";
    #$A=alternate_aa($a);
    #print "First mod codon\n".$A."\n";
    #$B=alternate_aa($b);
    #print "Second mod codon\n".$B."\n";
    #$C=alternate_aa($c);
    #print "Third mod codon\n".$C."\n";
    $D=alternate_aa($d);
    #print "Third mod codon\n".$C."\n";
    #$E=alternate_aa($e);
    #($rightori[($loc+$dif)],$rightori[($loc+$dif)+1],$rightori[($loc+$dif)+2])=split(//,$A);
    #($rightori[($loc+$dif)+3],$rightori[($loc+$dif+3)+1],$rightori[($loc+$dif+3)+2])=split(//,$B);
    #($rightori[($loc+$dif)+6],$rightori[($loc+$dif+6)+1],$rightori[($loc+$dif+6)+2])=split(//,$C);
    ($rightori[($loc+$dif)+9],$rightori[($loc+$dif+9)+1],$rightori[($loc+$dif+9)+2])=split(//,$D);
    #($rightori[($loc+$dif)+18],$rightori[($loc+$dif+18)+1],$rightori[($loc+$dif+18)+2])=split(//,$E);
    $RightOri=join('',@rightori);
}

#printf $data[$column-1]."\t".$crRNA."\t".$hash{$key."strand"}."\t".$info[11]."\t".$start."\t".$end."\t".join('',@seq[$start..($end-1)])."\t".$cutsite."\t".$seq[$cutsite]."\t".$frame[$cutsite]."\t".join('',@ori[$start..($end-1)])."\t".join('',@frame[$start..($end-1)])."\t".join('',@exon[$start..($end-1)])."\n";
printf $data[$column-1]."\t".$LeftOri."\t".$RightOri."\n";
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
if($HoA{$aa}[0]==$pat){return($HoA{$aa}[1])}else{return($HoA{$aa}[0])}
}
}
