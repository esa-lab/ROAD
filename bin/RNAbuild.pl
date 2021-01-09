#!/usr/bin/env perl   
#  -*- perl -*-
#
#	RNAbuild1.5pl
#
#	Written by Cody Geary, cody@dna.caltech.edu May 2016
#
#	Requires RNA_nts.pdb datafile
#
#	Usage:  perl RNAbuild.pl input.file
#	
#	Input file is a traceable text map with >Name as the first line
#	Outputs Name.pdb PDB file
#
#   Uses rigid-body rotations to build RNA structures from a dot-paren and a sequence.
#	Recognizes secondary structural patterns in the dot-paren as 3D motifs:
#			(....) maps to GAAA tetraloop	
#			(.........) (..[[[[[[.) and (..]]]]]].) map to 180KL with the sequence AANNNNNNA
#			^	maps to a crossover point with AE symmetry
#
#	Support for Spinach, Mango, PP7, MS2, Tar/Tat added
#	Support for L7Ae and bKL
	

#Copyright 2020  Cody Geary and Ebbe S. Andersen

#Permission is hereby granted, free of charge, to any person obtaining
#a copy of this software and associated documentation files (the
#"Software"), \ to deal in the Software without restriction, including
#without limitation the rights to use, copy, modify, merge, publish,
#distribute, sublicense, \ and/or sell copies of the Software, and to
#permit persons to whom the Software is furnished to do so, subject to
#the following conditions:

#The above copyright notice and this permission notice shall be
#included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY,\ FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\ LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALIN\ GS IN THE
# SOFTWARE.

use strict;
use warnings FATAL => qw ( all );
#use Data::Dumper;
use Math::Trig ':radial';
use Math::Trig;

my @nt_step = ();	my @angles = (); 	my @point = ();	my @moved = ();
my @ref_frame = ();	my @ref_frame_angle = ();	my @ref_trans = ();
my @motif = (); #is each nt_pos a motif, 1 or 0
my $roundedx =0; my $roundedy=0; my $roundedz=0;
my $ref_frame_position = 1;
my $end_of_chain = 0;

####################
###Begin Program
####################

#####regular version
my $file1 = "";	( $file1 ) = @ARGV;

#####server version
#my ( $file1, $file2, $file3 ) = @ARGV;
#$file2 //= "RNA_lib.pdb";
#$file3 //= "output.pdb";


sub read_file
{
    my ( $file ) = @_;

    my ( @lines, $content );

    if ( not open FILE, $file ) {
        die "Not readable: $file";
    }

    {
        # Reads the whole file by temporarily undefining the
        # record separator,
        
        local $/ = undef;
        $content = <FILE>;

        # Make Unix line ends from either Mac or DOS,

        $content =~ s/\r\n/\n/g;
        $content =~ s/\r/\n/g;
    }

    close FILE;

    @lines = split "\n", $content;

    return wantarray ? @lines : \@lines;
}


#################################
#Trace the Strandpath Input
#################################


my $input_structure = "";
my $input_sequence = "";
my $name = 'mytestoutput';

&trace;

my $find = " ";
my $replace = "";
$name =~ s/$find/$replace/g;
$name=$name.".pdb";

#################################
#Setup Output Spool
#################################


##### regular version
open(my $output_spool, '>', "$name"); ##Output buffer for the screen/file outputs.  Spool will hold all screen outputs and send them to a file at the end.

##### server version
#open(my $output_spool, '>', "$file3"); ##Output buffer for the screen/file outputs.  Spool will hold all screen outputs and send them to a file at the end.



my ( $file, $line, @lines, @cols, $seq, @seq, $ATOM, @ATOM );
my @map = (); #stores the dot_paren map
( $file ) = @ARGV;

my @pat = split(//, $input_structure);

#####regular version
if ( not open FILE, "< RNA_lib.pdb" ) { die "RNA_lib.pdb file not found!";}	#some error checking stuff


#####server version
#if ( not open FILE, "< $file2" ) { die "RNA library file not found!";}	#some error checking stuff


if ( defined $input_sequence ) { @seq = split(//, $input_sequence);}
else { die "no sequence provided";}

my @PDB = ();   my $PDB = 0;  #this will store the final PDB file that we build


###Scan in the PDB library file
$ATOM = 0;
@ATOM = ( );  #this hash stores the PDB coords of the nt and motif library
while ( $line = <FILE> ) {   #This parses the information from the FILE into $lines.  Each line that begins with ATOM is further parsed into substrings that are stored into vars of ATOM hash-file.
    @cols = ( );
    @cols = split(/\s+/, $line);
    if ( defined $cols[0] ) {
        if ( $cols[0] eq "ATOM") { 
            push @ATOM, {
                "a" => substr("$line", 13, 3),  #atom type
                "n" => substr("$line", 19, 1),  #residue type
                "p" => substr("$line", 17, 3),  #residue type (for proteins)
                "i" => int(substr("$line", 22, 4)),  #residue number
                "h" => substr("$line", 21, 1), #chain name
                "x" => substr("$line", 30, 8),
                "y" => substr("$line", 38, 8),
                "z" => substr("$line", 46, 8),
            }
        }
    }
}
close FILE || die "couldn't close file.";

###Additional Variable Declarations
my $i = 0;	my $x = 1;	my $y = 1;	my $z = 1;
my $c = $ATOM[0]->{i}; 


my $strand_length = scalar(@seq);
my $pattern_length = scalar(@pat);

&printer("REMARK -  File Generated by RNA_build.pl - written by Cody Geary \n");
&printer("REMARK - Building a total of $strand_length nts. \n");


#&printer("REMARK - Sequence 5- $input_sequence -3 \n");  #remarks for long sequences will crash PDB viewers

my $dot_paren = $input_structure;
$find = "^";
$replace = "";
$find = quotemeta $find; # escape regex metachars if present
$dot_paren =~ s/$find/$replace/g;

&map($dot_paren); #generate the nt pairing map with no ^ symbols for RNA building


for ($i=0; $i<$strand_length; $i++){$motif[$i]=1;} #pre-block all sites, and uncheck them as clean helix is produced.

my $parsed_struc = $input_structure;

$find = qr/.[\^]/;  #temporarily replace "?^" with "E" so the strand length is the correct one  This will be updated again with the correct pattern
$replace = "E";
$parsed_struc =~ s/$find/$replace/g;

my @temp_structure = split(//, $parsed_struc);

$find =    "))\[\[\[\[\[\[.))";     #bKL A
$replace = "V----------"; #$replace = quotemeta $replace;
my $offset = 0;
my $search_pos = index($parsed_struc, $find, $offset);
	while ($search_pos != -1) {  #find the first position of each case, and each subsequent case
		$temp_structure[$map[$search_pos]]="-";   #bKLB
		$temp_structure[$map[$search_pos]-1]="-";	
		$temp_structure[$map[$search_pos]-2]="-";
		$temp_structure[$map[$search_pos]-3]="W";	
			
		$offset = $search_pos + 1;
		$search_pos = index($parsed_struc, $find, $offset);

	}

	$parsed_struc = join "",(@temp_structure);
	$find = quotemeta $find;
	$parsed_struc =~ s/$find/$replace/g; # now go and swap $find with the equal-length fragment $replace
	@temp_structure = split(//, $parsed_struc);


$find =    "))\]\]\]\]\]\].))";     #bKL A - closing KL form
$replace = "V----------"; #$replace = quotemeta $replace;
$offset = 0;
$search_pos = index($parsed_struc, $find, $offset);
	while ($search_pos != -1) {  #find the first position of each case, and each subsequent case
		$temp_structure[$map[$search_pos]]="-";   #bKLB
		$temp_structure[$map[$search_pos]-1]="-";	
		$temp_structure[$map[$search_pos]-2]="-";
		$temp_structure[$map[$search_pos]-3]="W";	
			
		$offset = $search_pos + 1;
		$search_pos = index($parsed_struc, $find, $offset);

	}

	$parsed_struc = join "",(@temp_structure);
	$find = quotemeta $find;
	$parsed_struc =~ s/$find/$replace/g; # now go and swap $find with the equal-length fragment $replace
	@temp_structure = split(//, $parsed_struc);


$find =    ")).......))";     #bKL A unpaired
$replace = "V----------"; #$replace = quotemeta $replace;
$offset = 0;
$search_pos = index($parsed_struc, $find, $offset);
	while ($search_pos != -1) {  #find the first position of each case, and each subsequent case
		$temp_structure[$map[$search_pos]]="-";   #bKLB
		$temp_structure[$map[$search_pos]-1]="-";	
		$temp_structure[$map[$search_pos]-2]="-";
		$temp_structure[$map[$search_pos]-3]="W";			
		$offset = $search_pos + 1;
		$search_pos = index($parsed_struc, $find, $offset);
	}
	$parsed_struc = join "",(@temp_structure);
	$find = quotemeta $find;
	$parsed_struc =~ s/$find/$replace/g; # now go and swap $find with the equal-length fragment $replace
	@temp_structure = split(//, $parsed_struc);



$find =    "(.....(";     #90deg motif AACUA
$replace = "F------"; #$replace = quotemeta $replace;
	$offset = 0;
	$search_pos = index($parsed_struc, $find, $offset);
	while ($search_pos != -1) {  #find the first position of each case, and each subsequent case
		$temp_structure[$map[$search_pos]]="-";   #90deg motif
		$temp_structure[$map[$search_pos]-1]="G";		
		$offset = $search_pos + 1;
		$search_pos = index($parsed_struc, $find, $offset);
	}
	$parsed_struc = join "",(@temp_structure);
	$find = quotemeta $find;
	$parsed_struc =~ s/$find/$replace/g; # now go and swap $find with the equal-length fragment $replace
	@temp_structure = split(//, $parsed_struc);



$find =    ").....)";     #90deg motif AACUA
$replace = "F------"; #$replace = quotemeta $replace;
	$offset = 0;
	$search_pos = index($parsed_struc, $find, $offset);
	while ($search_pos != -1) {  #find the first position of each case, and each subsequent case
		$temp_structure[$map[$search_pos]]="-";   #90deg motif
		$temp_structure[$map[$search_pos]-1]="G";		
		$offset = $search_pos + 1;
		$search_pos = index($parsed_struc, $find, $offset);
	}
	$parsed_struc = join "",(@temp_structure);
	$find = quotemeta $find;
	$parsed_struc =~ s/$find/$replace/g; # now go and swap $find with the equal-length fragment $replace
	@temp_structure = split(//, $parsed_struc);


#character "-" will be ignored, used as fill-spacer in the diagram


$find =    "((((((.(.((((....)))))))))))";  #Tar/Tat RNA length28
$replace = "R---------------------------";
$find = quotemeta $find; # escape regex metachars if present
$parsed_struc =~ s/$find/$replace/g;

$find =    "(((((.((((......)))))))))";  #PP7  length25
$replace = "N------------------------";
$find = quotemeta $find;
$parsed_struc =~ s/$find/$replace/g;

$find =    "((.((....))))";  #MS2 length13
$replace = "M------------";
$find = quotemeta $find;
$parsed_struc =~ s/$find/$replace/g;

$find =    "(((..((............(.(";  #iSpinach A
$replace = "S---------------------";
$find = quotemeta $find;
$parsed_struc =~ s/$find/$replace/g;

$find =    ")))..))............).)";  #iSpinach A
$replace = "S---------------------";
$find = quotemeta $find;
$parsed_struc =~ s/$find/$replace/g;

$find =    "(((..((..............(";  #iSpinach A
$replace = "S---------------------";
$find = quotemeta $find;
$parsed_struc =~ s/$find/$replace/g;

$find =    ")))..))..............)";  #iSpinach A
$replace = "S---------------------";
$find = quotemeta $find;
$parsed_struc =~ s/$find/$replace/g;


$find =    ").).........)))))";  #iSpinach B
$replace = "U----------------";
$find = quotemeta $find;
$parsed_struc =~ s/$find/$replace/g;

$find =    "(.(.........(((((";  #iSpinach B
$replace = "U----------------";
$find = quotemeta $find;
$parsed_struc =~ s/$find/$replace/g;


$find =    "(((.......................)))";  #Mango  length29
$replace = "Q----------------------------";
$find = quotemeta $find;
$parsed_struc =~ s/$find/$replace/g;




$find =    "(....)";
$replace = "T-----";
$find = quotemeta $find;
$parsed_struc =~ s/$find/$replace/g;

$find =    "(.........)";
$replace = "K----------";
$find = quotemeta $find;
$parsed_struc =~ s/$find/$replace/g;

$find =    "(..[[[[[[.)";
$replace = "K----------";
$find = quotemeta $find;
$parsed_struc =~ s/$find/$replace/g;

$find =    "(..]]]]]].)";
$replace = "K----------";
$find = quotemeta $find;
$parsed_struc =~ s/$find/$replace/g;

$find =    "(.......)";
$replace = "L--------";
$find = quotemeta $find;
$parsed_struc =~ s/$find/$replace/g;

$find =    "([[[[[[[)";
$replace = "L--------";
$find = quotemeta $find;
$parsed_struc =~ s/$find/$replace/g;

$find =    "(]]]]]]])";
$replace = "L--------";
$find = quotemeta $find;
$parsed_struc =~ s/$find/$replace/g;




$find =    "(...((";   #search for K-turn, note that this actually reads any 3nt or 6nt bulge as a Kturn pattern, but no other motifs use 3nt or 6nt at this stage.
$replace = "B-----";
$find = quotemeta $find;
$parsed_struc =~ s/$find/$replace/g;

$find =    ")...))";
$replace = "B-----";
$find = quotemeta $find;
$parsed_struc =~ s/$find/$replace/g;

$find =    "((......(";  # B is the 2nd side of the Kturn
$replace = "C--------";
$find = quotemeta $find;
$parsed_struc =~ s/$find/$replace/g;

$find =    "))......)";
$replace = "C--------";
$find = quotemeta $find;
$parsed_struc =~ s/$find/$replace/g;

#&printer("$parsed_struc \n");  ##Structures larger than 1000nts can crash SwissPDB by overflowing the remark buffer 

$find = qr[E..];
$replace = "E"; $replace = quotemeta $replace;
$parsed_struc =~ s/$find/$replace/g;


$find = "-"; #remove all of the "-" spacers now that we are at the end
$replace = "";
$parsed_struc =~ s/$find/$replace/g;


$find = qr/\W/;
$replace = "H";
$parsed_struc =~ s/$find/$replace/g;

$find = qr[\.];
$replace = "H";
$parsed_struc =~ s/$find/$replace/g;

my $nt_pos = 1; 



#&printer("$parsed_struc \n");  ##Structures larger than 1000nts can crash SwissPDB by overflowing the remark buffer 
#&printer("$input_structure \n");  ##Strutures larger than 1000nts can crash SwissPDB by overflowing the remark buffer 



&map($dot_paren); #generate the nt pairing map without the ^ symbols for RNA building

my @assembly_instructions = split(//, $parsed_struc);
my $k =0;
my $nt_disp = $nt_pos+1;
my $nt_disp2 = $nt_pos;
for ( $k=0; $k< scalar(@assembly_instructions) ; $k++){

	if ($assembly_instructions[$k] eq 'H') { &add_nts;}
	if ($assembly_instructions[$k] eq 'T') {&printer("REMARK Add GNRA $nt_pos-"); &add_tetraloop; &printer("$nt_pos\n");}
	if ($assembly_instructions[$k] eq 'K') {&printer("REMARK Add 180KL $nt_pos-"); &add_KL; $nt_disp=$nt_pos-1;   &printer("$nt_disp\n");}
	if ($assembly_instructions[$k] eq 'L') {&printer("REMARK Add 120KL $nt_pos-"); &add_KLbend;&printer("$nt_pos\n");}
	if ($assembly_instructions[$k] eq 'E') {
		$nt_disp =$map[$nt_pos-1]+1; &printer("REMARK Add Crossover $nt_pos,$nt_disp");
		&crossover;
		$nt_disp=$nt_pos-2; $nt_disp2 =$map[$nt_disp-1]+1; &printer(",$nt_disp,$nt_disp2\n");
	}
	
	if ($assembly_instructions[$k] eq 'R') {&printer("REMARK Add TAR $nt_pos-"); &add_tar;&printer("$nt_pos\n");}
	if ($assembly_instructions[$k] eq 'M') {&printer("REMARK Add MS2 $nt_pos-"); &add_ms2;&printer("$nt_pos\n");}
	if ($assembly_instructions[$k] eq 'N') {&printer("REMARK Add PP7 $nt_pos-"); &add_pp7;&printer("$nt_pos\n");}
	
	
	if ($assembly_instructions[$k] eq 'S') {&printer("REMARK Add Spinach-A $nt_pos-"); &add_spinacha;&printer("$nt_pos\n");}
	if ($assembly_instructions[$k] eq 'U') {&printer("REMARK Add Spinach-B $nt_pos-"); &add_spinachb;&printer("$nt_pos\n");}

	if ($assembly_instructions[$k] eq 'B') {&printer("REMARK Add Kturn-A $nt_pos-"); &add_kturna;&printer("$nt_pos\n");}
	if ($assembly_instructions[$k] eq 'C') {&printer("REMARK Add Kturn-B $nt_pos-"); &add_kturnb;&printer("$nt_pos\n");}

	if ($assembly_instructions[$k] eq 'Q') {&printer("REMARK Add MANGO $nt_pos-"); &add_mango;&printer("$nt_pos\n");}	

	if ($assembly_instructions[$k] eq 'V') {&printer("REMARK Add BKL-A $nt_pos-"); &add_bkla; &printer("$nt_pos\n");}
	if ($assembly_instructions[$k] eq 'W') {&printer("REMARK Add BKL-B $nt_pos-"); &add_bklb; &printer("$nt_pos\n");}

	if ($assembly_instructions[$k] eq 'F') {&printer("REMARK Add 90deg AACUA bend\n"); &add_ninetya;}
	if ($assembly_instructions[$k] eq 'G') {&printer("REMARK Add 90deg AACUA hinge\n"); &add_ninetyb;}

		
}
&printer("REMARK Assembly completed\n");
####Print PDB File


my $t = 1;
my $seam = 0;
my $filecount=1;
foreach $PDB ( @PDB ) {  #print out the PDB file
	if ($t<99000){$seam = $PDB->{i};}
	if ($PDB->{i}<=$seam){
	
		if (defined $PDB->{n}){  # PDB->{n} or {p} tells if it is RNA or Protein, Here we only put RNAs 	
			printf $output_spool "%-6s", "ATOM";
			printf $output_spool "%5s", $t;
			print $output_spool "  ";
			print $output_spool $PDB->{a};

			print $output_spool "   ";
			print $output_spool $PDB->{n};
		
			print $output_spool " ";
			print $output_spool $PDB->{h};  #Retain chain names
			printf $output_spool "%4s", $PDB->{i};
			print $output_spool "    ";					
			printf $output_spool ("%8s",$PDB->{x});
			printf $output_spool ("%8s",$PDB->{y});
			printf $output_spool ("%8s",$PDB->{z});
			print $output_spool "  1.00100.00   ";
			print $output_spool "\n";
			$t++;
		}
		
	} else {
		close $output_spool;
	
		substr $name, length($name)-4 , 4, "";
		$name=$name."_$filecount.pdb";
		
		open( $output_spool, '>', "$name"); #open a new file buffer
		$t=1;
		
		if (defined $PDB->{n}){  # PDB->{n} or {p} tells if it is RNA or Protein, Here we only put RNAs
			printf $output_spool "%-6s", "ATOM";
			printf $output_spool "%5s", $t;
			print $output_spool "  ";
			print $output_spool $PDB->{a};
			print $output_spool "   ";
			print $output_spool $PDB->{n};
			print $output_spool " ";
			print $output_spool $PDB->{h};  #Retain chain names
			printf $output_spool "%4s", $PDB->{i};
			print $output_spool "    ";					
			printf $output_spool ("%8s",$PDB->{x});
			printf $output_spool ("%8s",$PDB->{y});
			printf $output_spool ("%8s",$PDB->{z});
			print $output_spool "  1.00100.00   ";
			print $output_spool "\n";
			$t++;
		}
	}
}
printf $output_spool "TER\n";

my $chain_trigger=0;

foreach $PDB ( @PDB ) {  #print out the PDB file non-RNA part
	
	if (!defined $PDB->{n}){	# PDB->{n} or {p} tells if it is RNA or Protein, Here we only put Protein

		if ( $PDB->{h} eq "P" ) {
			
			if ( $PDB->{t} eq "0" && $chain_trigger==1){
				$chain_trigger=0;
				printf $output_spool "TER\n";
			}
			
			printf $output_spool "%-6s", "ATOM";
			printf $output_spool "%5s", $t;
			print $output_spool "  ";
			print $output_spool $PDB->{a};

			print $output_spool " ";
			print $output_spool $PDB->{p};
	
			print $output_spool " ";
			print $output_spool $PDB->{h};  #Retain chain names
			printf $output_spool "%4s", $PDB->{i};
			print $output_spool "    ";					
			printf $output_spool ("%8s",$PDB->{x});
			printf $output_spool ("%8s",$PDB->{y});
			printf $output_spool ("%8s",$PDB->{z});
			print $output_spool "  1.00100.00   ";
			print $output_spool "\n";
			
			if ( $PDB->{t} eq "1" ){ $chain_trigger=1;}
			
			$t++;
		} 
	}
		
}



close $output_spool;

#################################################
## >>>>>>>>>>>>>>>>>> SUB ROUTINES <<<<<<<<<<<<<<<<<<<<
#################################################

sub translate_matrix {
	my @r_tr1 = split(' ',$_[0]);
	my @r_tr2 = split(' ',$_[1]);
	my $tx=$r_tr2[0];
	my $ty=$r_tr2[1];
	my $tz=$r_tr2[2];
	my @Translate=(
	    1, 0, 0, 0,
    	0, 1, 0, 0, 
    	0, 0, 1, 0,
    	$tx, $ty, $tz, 1);
	transform ("@r_tr1","@Translate");
}

sub rotate_x {  #point followed by theta in rad
	my @r_tr1 = split(' ',$_[0]);
	my $theta = $_[1];
	my @Translate=(
	    1, 0, 0, 0,
    	0, cos($theta), -sin($theta), 0, 
    	0, sin($theta), cos($theta), 0,
    	0, 0, 0, 1);
	transform ("@r_tr1","@Translate");
}

sub rotate_y {  #point followed by theta in rad
	my @r_tr1 = split(' ',$_[0]);
	my $theta = $_[1];
	my @Translate=(
	    cos($theta), 0, sin($theta), 0,
    	0, 1, 0, 0, 
    	-sin($theta), 0, cos($theta), 0,
    	0, 0, 0, 1);
	transform ("@r_tr1","@Translate");
}

sub rotate_z {  #point followed by theta in rad
	my @r_tr1 = split(' ',$_[0]);
	my $theta = $_[1];
	my @Translate=(
	    cos($theta), -sin($theta), 0, 0,
    	sin($theta), cos($theta), 0, 0, 
    	0, 0, 1, 0,
    	0, 0, 0, 1);
	transform ("@r_tr1","@Translate");
}

sub transform {   #point (x,y,z,1) multiplied by the transform_matrix (16 elements long)
	my @m1 = split(' ',$_[0]);
	my @m2 = split(' ',$_[1]);
    my $result = [];

	my $xout=$m1[0]*$m2[0] + $m1[1]*$m2[4] + $m1[2]*$m2[8] + $m1[3]*$m2[12];
	my $yout=$m1[0]*$m2[1] + $m1[1]*$m2[5] + $m1[2]*$m2[9] + $m1[3]*$m2[13];
	my $zout=$m1[0]*$m2[2] + $m1[1]*$m2[6] + $m1[2]*$m2[10] + $m1[3]*$m2[14];
	my $vout=$m1[0]*$m2[3] + $m1[1]*$m2[7] + $m1[2]*$m2[11] + $m1[3]*$m2[15];  
	($xout, $yout, $zout, $vout);
}

sub add_nts {
	my $no = 1;
	if ($map[$nt_pos-1]<($nt_pos-1) && $motif[$map[$nt_pos-1]]==0) {  #if map<nt_pos then we are a closing nt, and we calculate from the partner's ref-frame (later we may average the two options)

		$ref_frame_position = $map[$nt_pos-1]+2;
		@ref_frame_angle = &update_ref_frame($ref_frame_position);
							
		foreach $ATOM ( @ATOM ) {			#paste in the nt as paired to the map[nt_pos]
			if ( ($ATOM->{n} eq $seq[$nt_pos-1]) && ($ATOM->{h} eq 'A')) {   #$c is the current resi, atom(i) is the resi number element in the PDB file.							

				$nt_step[0]=2.2260;	$nt_step[1]=2.0800;	$nt_step[2]=-10.4200;	$nt_step[3]=1;			#flip strand
				$angles[0]=2.00510013;	$angles[1]=2.95283932;	$angles[2]=1.1652;		

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"n" => $ATOM->{n},  #residue type
					"i" => $nt_pos,  #residue number
					"h" => "A", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),																	
				}			
			}
		}
		$nt_pos += 1;	#increment nt_counter		
	
	} else {

		$ref_frame_position = $nt_pos;
		@ref_frame_angle = &update_ref_frame($ref_frame_position);
							
		foreach $ATOM ( @ATOM ) {			#we scan the datafile for the nt type and paste it into the active site
			if ( ($ATOM->{n} eq $seq[$nt_pos-1]) && ($ATOM->{h} eq 'A')) {   #$c is the current resi, atom(i) is the resi number element in the PDB file.							

				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"n" => $ATOM->{n},  #residue type
					"i" => $nt_pos,  #residue number
					"h" => "A", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),																	
				}			
			}
		}
		$nt_pos += 1;	#increment nt_counter
	}
	
	$motif[$nt_pos-2] = 0;
}

sub update_ref_frame {  #takes ($nt_pos) and aligns to that sugar
	my $c_one_x=0; my $c_one_y=0; my $c_one_z=0;  #coords of all the C1-primes
	my $c_four_x=0; my $c_four_y=0; my $c_four_z=0;  #coords of all the C4-primes
	my $n_x=0; my $n_y=0; my $n_z=0;  #coords of all the N1 or N9

	foreach $PDB ( @PDB ) {
	
		if ( $PDB->{i} == ($_[0]-1) ) {   #$nt_pos is the current position that we grab the ref frame of here
			if ( $PDB->{a} eq "O3*" || $PDB->{a} eq "O3'") {   #look for two different spellings of O3'
				$x = $PDB->{x};
				$y = $PDB->{y};
				$z = $PDB->{z};
			}
			if ( $PDB->{a} eq "C1*" || $PDB->{a} eq "C1'") {   #look for two different spellings of C1'
				$c_one_x = $PDB->{x};
				$c_one_y = $PDB->{y};
				$c_one_z = $PDB->{z};
			}
			if ( $PDB->{a} eq "C4*" || $PDB->{a} eq "C4'") {   #look for two different spellings of C4'
				$c_four_x = $PDB->{x};
				$c_four_y = $PDB->{y};
				$c_four_z = $PDB->{z};
			}
			if ( $PDB->{a} eq "C2*" || $PDB->{a} eq "C2'") {   #look for atoms
				$n_x = $PDB->{x};
				$n_y = $PDB->{y};
				$n_z = $PDB->{z};
			}			
		}
	}	
	&calc_ref_frame ($c_one_x, $c_one_y, $c_one_z, $c_four_x, $c_four_y, $c_four_z, $n_x, $n_y, $n_z);  #calculate the ref frame of the current nt		
} 

sub calc_ref_frame {  #needs (c1'x, c1'y, c1'z, c4'x, c4'y, c4'z, c3'x, c3'y, c3'z)
	my @point2 = ();   #this is C1'
	$point2[0]=$_[0];	$point2[1]=$_[1];	$point2[2]=$_[2];	$point2[3]=1;

	my @point3 = ();   #this is C4'
	$point3[0]=$_[3];	$point3[1]=$_[4];	$point3[2]=$_[5];	$point3[3]=1;
	
	my @point4 = ();   #this is C3'
	$point4[0]=$_[6];	$point4[1]=$_[7];	$point4[2]=$_[8];	$point4[3]=1;

	#Calculate the translation vector we need	
	my @vector1 = ();   #translates C1' to the origin
	$vector1[0]=-$_[0];	$vector1[1]=-$_[1];	$vector1[2]=-$_[2];	$vector1[3]=1;
	
	my @c4vect = &translate_matrix ("@point3","@vector1");	#look at C4' relative to C1'	
	my($rho, $theta, $phi)  = cartesian_to_spherical($c4vect[0], $c4vect[1], $c4vect[2]); #calculate the rotation angles we need to put C4 on the x-axis

	my @nvect = &translate_matrix ("@point4","@vector1"); #look at N(now c3') relative to C1, so we can follow and calculate the y-axis rotation	
	@nvect = &rotate_z ("@nvect",$theta);  #first rotation
	@nvect = &rotate_y ("@nvect",$phi); #second rotation

	my $last_angle = 0;		
	if ($nvect[0]==0 && $nvect[1]==0){  #avoid divide by zero error for nts on the axis.
		$last_angle = 0;  #at the asymptote we get 0
	} elsif ($nvect[0]==0 && ($nvect[1]>0 || $nvect[1]<0) ) {
		$last_angle = 3.141592/2;  #at the asymptote we get pi/2
	} else {
		$last_angle = atan($nvect[1]/$nvect[0]);
	}
	@nvect = &rotate_z ("@nvect",$last_angle); #check the position of the c3' after the last move
	if ($nvect[0]<0){$last_angle += 3.141592654;}  #there are 2 possible orientations in the plane, so we flip all the negative ones

	$ref_trans[0]=$_[0];	$ref_trans[1]=$_[1];	$ref_trans[2]=$_[2];

	($theta, $phi, $last_angle);  #these angles correspond to rotate_z, rotate_y, and rotate_z in that order
}

####################################
### MOTIFS
####################################

sub add_tetraloop {

	if ($seq[$nt_pos] eq "G"){	
		#&printer("REMARK     Loop is a GNRA shown as GAAA\n");
		
		
		$ref_frame_position = $nt_pos;
		@ref_frame_angle = &update_ref_frame($ref_frame_position);
		foreach $ATOM ( @ATOM ) {			#Copy the backbone of resi 1 of this motif,  o2' has to come at the end, so we pop it off at the end.
			if ( ($ATOM->{h} eq 'T')  &&  ($ATOM->{i} == 1) )	{							
				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"n" => $seq[$nt_pos-1],  #residue type
					"i" => $nt_pos,  #residue number
					"h" => "A", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),															
				}			
			}		
		}	
		my $nt_temp = pop @PDB;		#remove the O2' atom

		foreach $ATOM ( @ATOM ) {						#copy only the base atoms from the lib
			if ( ( ($ATOM->{h} eq 'A')  &&  ($ATOM->{n} eq $seq[$nt_pos-1]) &&
				( ($ATOM->{a} eq 'N9 ') || ($ATOM->{a} eq 'C8 ') ||
				  ($ATOM->{a} eq 'N7 ') || ($ATOM->{a} eq 'C5 ') ||
				  ($ATOM->{a} eq 'C6 ') || ($ATOM->{a} eq 'N6 ') ||
				  ($ATOM->{a} eq 'N1 ') || ($ATOM->{a} eq 'C2 ') ||
				  ($ATOM->{a} eq 'N3 ') || ($ATOM->{a} eq 'C4 ') ||
				  ($ATOM->{a} eq 'O6 ') || ($ATOM->{a} eq 'N2 ') ||
				  ($ATOM->{a} eq 'O2 ') || ($ATOM->{a} eq 'N4 ') || ($ATOM->{a} eq 'O4 ')) ) ) {   	
			  
				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"n" => $seq[$nt_pos-1],  #residue type
					"i" => $nt_pos,  #residue number
					"h" => "A", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),																
				}			
			}
		}		
		push @PDB, $nt_temp;	#add the O2' atom back at the end
		for ($i=2;$i<6;$i++) {  
			foreach $ATOM ( @ATOM ) {			#we scan the datafile for the nt type and paste it in
				if ( ($ATOM->{h} eq 'T') &&  ($ATOM->{i} == $i)) {   #$c is the current resi, atom(i) is the resi number element in the PDB file.							
					$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
					$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

					$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

					@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
					@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
					@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
					@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

					@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
					@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
					@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

					@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

					$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

					push @PDB, {
						"a" => $ATOM->{a},  #atom type
						"n" => $ATOM->{n},  #residue type
						"i" => $nt_pos-1+$i,  #residue number
						"h" => "A", #chain name		
						"x" => sprintf("%8s",$roundedx),
						"y" => sprintf("%8s",$roundedy),
						"z" => sprintf("%8s",$roundedz),																	
					}			
				}
			}	
		}
	
		$nt_pos += 5;
		
		foreach $ATOM ( @ATOM ) {			#Copy the backbone of resi 1 of this motif,  o2' has to come at the end, so we pop it off at the end.
			if ( ($ATOM->{h} eq 'T')  &&  ($ATOM->{i} == 6) )	{							

				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"n" => $seq[$nt_pos-1],  #residue type
					"i" => $nt_pos,  #residue number
					"h" => "A", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),																	
				}			
			}		
		}	
		$nt_temp = pop @PDB;		#remove the O2' atom
	
	
		$ref_frame_position = $nt_pos+1;
		@ref_frame_angle = &update_ref_frame($ref_frame_position);
	
	
		foreach $ATOM ( @ATOM ) {						#copy only the base atoms from the lib
			if ( ( ($ATOM->{h} eq 'A')  &&  ($ATOM->{n} eq $seq[$nt_pos-1]) &&
				( ($ATOM->{a} eq 'N9 ') || ($ATOM->{a} eq 'C8 ') ||
				  ($ATOM->{a} eq 'N7 ') || ($ATOM->{a} eq 'C5 ') ||
				  ($ATOM->{a} eq 'C6 ') || ($ATOM->{a} eq 'N6 ') ||
				  ($ATOM->{a} eq 'N1 ') || ($ATOM->{a} eq 'C2 ') ||
				  ($ATOM->{a} eq 'N3 ') || ($ATOM->{a} eq 'C4 ') ||
				  ($ATOM->{a} eq 'O6 ') || ($ATOM->{a} eq 'N2 ') ||
				  ($ATOM->{a} eq 'O2 ') || ($ATOM->{a} eq 'N4 ') || ($ATOM->{a} eq 'O4 ')) ) ) {   						
					$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

					@moved = &rotate_z ("@point",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
					@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
					@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

					@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

					$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

					push @PDB, {
						"a" => $ATOM->{a},  #atom type
						"n" => $seq[$nt_pos-1],  #residue type
						"i" => $nt_pos,  #residue number
						"h" => "A", #chain name		
						"x" => sprintf("%8s",$roundedx),
						"y" => sprintf("%8s",$roundedy),
						"z" => sprintf("%8s",$roundedz),																	
					}			
			}
		}		
		push @PDB, $nt_temp;	#add the O2' atom back at the end
	
		$nt_pos += 1;

		
	} else {
		#&printer("REMARK     Loop is a UNCG shown as UUCG\n");
		
		
		$ref_frame_position = $nt_pos;
		@ref_frame_angle = &update_ref_frame($ref_frame_position);
		foreach $ATOM ( @ATOM ) {			#Copy the backbone of resi 1 of this motif,  o2' has to come at the end, so we pop it off at the end.
			if ( ($ATOM->{h} eq 'U')  &&  ($ATOM->{i} == 1) )	{							
				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"n" => $seq[$nt_pos-1],  #residue type
					"i" => $nt_pos,  #residue number
					"h" => "A", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),															
				}			
			}		
		}	
		my $nt_temp = pop @PDB;		#remove the O2' atom

		foreach $ATOM ( @ATOM ) {						#copy only the base atoms from the lib
			if ( ( ($ATOM->{h} eq 'A')  &&  ($ATOM->{n} eq $seq[$nt_pos-1]) &&
				( ($ATOM->{a} eq 'N9 ') || ($ATOM->{a} eq 'C8 ') ||
				  ($ATOM->{a} eq 'N7 ') || ($ATOM->{a} eq 'C5 ') ||
				  ($ATOM->{a} eq 'C6 ') || ($ATOM->{a} eq 'N6 ') ||
				  ($ATOM->{a} eq 'N1 ') || ($ATOM->{a} eq 'C2 ') ||
				  ($ATOM->{a} eq 'N3 ') || ($ATOM->{a} eq 'C4 ') ||
				  ($ATOM->{a} eq 'O6 ') || ($ATOM->{a} eq 'N2 ') ||
				  ($ATOM->{a} eq 'O2 ') || ($ATOM->{a} eq 'N4 ') || ($ATOM->{a} eq 'O4 ')) ) ) {   	
			  
				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"n" => $seq[$nt_pos-1],  #residue type
					"i" => $nt_pos,  #residue number
					"h" => "A", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),																
				}			
			}
		}		
		push @PDB, $nt_temp;	#add the O2' atom back at the end
		for ($i=2;$i<6;$i++) {  
			foreach $ATOM ( @ATOM ) {			#we scan the datafile for the nt type and paste it in
				if ( ($ATOM->{h} eq 'U') &&  ($ATOM->{i} == $i)) {   #$c is the current resi, atom(i) is the resi number element in the PDB file.							
					$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
					$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

					$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

					@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
					@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
					@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
					@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

					@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
					@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
					@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

					@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

					$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

					push @PDB, {
						"a" => $ATOM->{a},  #atom type
						"n" => $ATOM->{n},  #residue type
						"i" => $nt_pos-1+$i,  #residue number
						"h" => "A", #chain name		
						"x" => sprintf("%8s",$roundedx),
						"y" => sprintf("%8s",$roundedy),
						"z" => sprintf("%8s",$roundedz),																	
					}			
				}
			}	
		}
	
		$nt_pos += 5;
		
		foreach $ATOM ( @ATOM ) {			#Copy the backbone of resi 1 of this motif,  o2' has to come at the end, so we pop it off at the end.
			if ( ($ATOM->{h} eq 'U')  &&  ($ATOM->{i} == 6) )	{							

				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"n" => $seq[$nt_pos-1],  #residue type
					"i" => $nt_pos,  #residue number
					"h" => "A", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),																	
				}			
			}		
		}	
		$nt_temp = pop @PDB;		#remove the O2' atom
	
	
		$ref_frame_position = $nt_pos+1;
		@ref_frame_angle = &update_ref_frame($ref_frame_position);
	
	
		foreach $ATOM ( @ATOM ) {						#copy only the base atoms from the lib
			if ( ( ($ATOM->{h} eq 'A')  &&  ($ATOM->{n} eq $seq[$nt_pos-1]) &&
				( ($ATOM->{a} eq 'N9 ') || ($ATOM->{a} eq 'C8 ') ||
				  ($ATOM->{a} eq 'N7 ') || ($ATOM->{a} eq 'C5 ') ||
				  ($ATOM->{a} eq 'C6 ') || ($ATOM->{a} eq 'N6 ') ||
				  ($ATOM->{a} eq 'N1 ') || ($ATOM->{a} eq 'C2 ') ||
				  ($ATOM->{a} eq 'N3 ') || ($ATOM->{a} eq 'C4 ') ||
				  ($ATOM->{a} eq 'O6 ') || ($ATOM->{a} eq 'N2 ') ||
				  ($ATOM->{a} eq 'O2 ') || ($ATOM->{a} eq 'N4 ') || ($ATOM->{a} eq 'O4 ')) ) ) {   						
					$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

					@moved = &rotate_z ("@point",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
					@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
					@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

					@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

					$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

					push @PDB, {
						"a" => $ATOM->{a},  #atom type
						"n" => $seq[$nt_pos-1],  #residue type
						"i" => $nt_pos,  #residue number
						"h" => "A", #chain name		
						"x" => sprintf("%8s",$roundedx),
						"y" => sprintf("%8s",$roundedy),
						"z" => sprintf("%8s",$roundedz),																	
					}			
			}
		}		
		push @PDB, $nt_temp;	#add the O2' atom back at the end
	
		$nt_pos += 1;
		
	}
}

sub crossover {	
	$ref_frame_position = $nt_pos-1;
	@ref_frame_angle = &update_ref_frame($ref_frame_position);
	
	foreach $ATOM ( @ATOM ) {			#Copy the backbone of resi 1 of this motif,  o2' has to come at the end, so we pop it off at the end.
		if ( ($ATOM->{h} eq 'E')  &&  ($ATOM->{i} == 2) )	{	#we start with nt 2,and moved back 1 frame						
			$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

			$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
			$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

			@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
			@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
			@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

			@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
			@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	
			@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

			$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

			push @PDB, {
				"a" => $ATOM->{a},  #atom type
				"n" => $seq[$nt_pos-1],  #residue type
				"i" => $nt_pos,  #residue number
				"h" => "A", #chain name		
				"x" => sprintf("%8s",$roundedx),
				"y" => sprintf("%8s",$roundedy),
				"z" => sprintf("%8s",$roundedz),															
			}			
		}		
	}	
	my $nt_temp = pop @PDB;		#remove the O2' atom  (before rotation since the crossover ref frame is offset by 1)
	
	foreach $ATOM ( @ATOM ) {						#copy only the base atoms from the lib
		if ( ( ($ATOM->{h} eq 'A')  &&  ($ATOM->{n} eq $seq[$nt_pos-1]) &&
			( ($ATOM->{a} eq 'N9 ') || ($ATOM->{a} eq 'C8 ') ||
			  ($ATOM->{a} eq 'N7 ') || ($ATOM->{a} eq 'C5 ') ||
			  ($ATOM->{a} eq 'C6 ') || ($ATOM->{a} eq 'N6 ') ||
			  ($ATOM->{a} eq 'N1 ') || ($ATOM->{a} eq 'C2 ') ||
			  ($ATOM->{a} eq 'N3 ') || ($ATOM->{a} eq 'C4 ') ||
			  ($ATOM->{a} eq 'O6 ') || ($ATOM->{a} eq 'N2 ') ||
			  ($ATOM->{a} eq 'O2 ') || ($ATOM->{a} eq 'N4 ') || ($ATOM->{a} eq 'O4 ') ) ) ) {   						
			$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

			$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
			$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

			@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
			@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
			@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

			@moved = &rotate_z ("@moved",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
			@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
			@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

			@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
			@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

			@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

			$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

			push @PDB, {
				"a" => $ATOM->{a},  #atom type
				"n" => $seq[$nt_pos-1],  #residue type
				"i" => $nt_pos,  #residue number
				"h" => "A", #chain name		
				"x" => sprintf("%8s",$roundedx),
				"y" => sprintf("%8s",$roundedy),
				"z" => sprintf("%8s",$roundedz),																
			}			
		}
	}		
	$nt_pos += 1;	
	push @PDB, $nt_temp;	#add the O2' atom back at the end	
	
	foreach $ATOM ( @ATOM ) {			#Copy the backbone of resi 2 of this motif,  o2' has to come at the end, so we pop it off at the end.
		if ( ($ATOM->{h} eq 'E')  &&  ($ATOM->{i} == 3) )	{							
			$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

			$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
			$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

			@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
			@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
			@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

			@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
			@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

			@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

			$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

			push @PDB, {
				"a" => $ATOM->{a},  #atom type
				"n" => $seq[$nt_pos-1],  #residue type
				"i" => $nt_pos,  #residue number
				"h" => "A", #chain name		
				"x" => sprintf("%8s",$roundedx),
				"y" => sprintf("%8s",$roundedy),
				"z" => sprintf("%8s",$roundedz),					
			}			
		}		
	}
	$nt_temp = pop @PDB;		#remove the O2' atom
	#now, recalculate a reference frame from residue 4, and step back 1 nt_step to place the crossover base.
	my $c_one_x=0; my $c_one_y=0; my $c_one_z=0;  #coords of all the C1-primes
	my $c_four_x=0; my $c_four_y=0; my $c_four_z=0;  #coords of all the C4-primes
	my $n_x=0; my $n_y=0; my $n_z=0;  #coords of all the N1 or N9

	foreach $ATOM ( @ATOM ) {	
		if ( ($ATOM->{i} == 4) && ($ATOM->{h} eq 'E') ) {   #$nt_pos is the current position that we grab the ref frame of here

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector
				
				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	
				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 
			

		
			if ( $ATOM->{a} eq "C1*" || $ATOM->{a} eq "C1'") {   #look for two different spellings of C1'			
				$c_one_x = $moved[0];
				$c_one_y = $moved[1];
				$c_one_z = $moved[2];
			}
			if ( $ATOM->{a} eq "C4*" || $ATOM->{a} eq "C4'") {   #look for two different spellings of C4'
				$c_four_x = $moved[0];
				$c_four_y = $moved[1];
				$c_four_z = $moved[2];
			}
			if ( $ATOM->{a} eq "C2*" || $ATOM->{a} eq "C2'") {   #look for atoms
				$n_x = $moved[0];
				$n_y = $moved[1];
				$n_z = $moved[2];
			}			
		}
	}	
	my @base_ref_frame_angle = &calc_ref_frame ($c_one_x, $c_one_y, $c_one_z, $c_four_x, $c_four_y, $c_four_z, $n_x, $n_y, $n_z);  #calculate the ref frame of the current nt		
	my @base_ref_trans = @ref_trans;
	
	foreach $ATOM ( @ATOM ) {						#copy only the base atoms from the lib
		if ( ( ($ATOM->{h} eq 'A')  &&  ($ATOM->{n} eq $seq[$nt_pos-1]) &&
			( ($ATOM->{a} eq 'N9 ') || ($ATOM->{a} eq 'C8 ') ||
			  ($ATOM->{a} eq 'N7 ') || ($ATOM->{a} eq 'C5 ') ||
			  ($ATOM->{a} eq 'C6 ') || ($ATOM->{a} eq 'N6 ') ||
			  ($ATOM->{a} eq 'N1 ') || ($ATOM->{a} eq 'C2 ') ||
			  ($ATOM->{a} eq 'N3 ') || ($ATOM->{a} eq 'C4 ') ||
			  ($ATOM->{a} eq 'O6 ') || ($ATOM->{a} eq 'N2 ') ||
			  ($ATOM->{a} eq 'O2 ') || ($ATOM->{a} eq 'N4 ') || ($ATOM->{a} eq 'O4 ')) ) ) {  
			   						
			$moved[0]=$ATOM->{x};	$moved[1]=$ATOM->{y};	$moved[2]=$ATOM->{z};	$moved[3]=1;	

			$nt_step[0]=-5.05026531;	$nt_step[1]=-0.63351020;	$nt_step[2]=2.27143878;	$nt_step[3]=1;	#take a forward step
			$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

			@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			
			@moved = &rotate_z ("@moved",$angles[0]);  #first rotation  -we move the nt to the untranslocated pos
			@moved = &rotate_y ("@moved",$angles[1]); #second rotation	
			@moved = &rotate_z ("@moved",$angles[2]); #final rotation

			@moved = &rotate_z ("@moved",-$base_ref_frame_angle[2]);  #first rotation    -first move the new nt to be in the frame of the last nt of motif
			@moved = &rotate_y ("@moved",-$base_ref_frame_angle[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$base_ref_frame_angle[0]); #final rotation	
			@moved = &translate_matrix ("@moved","@base_ref_trans");  #      Now moving to add to the 3prime end 
						
			$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

			push @PDB, {
				"a" => $ATOM->{a},  #atom type
				"n" => $seq[$nt_pos-1],  #residue type
				"i" => $nt_pos,  #residue number
				"h" => "A", #chain name		
				"x" => sprintf("%8s",$roundedx),
				"y" => sprintf("%8s",$roundedy),
				"z" => sprintf("%8s",$roundedz),											
			}			
		}
	}		
	push @PDB, $nt_temp;	#add the O2' atom back at the end
	$nt_pos += 1;
	
	$ref_frame_position = $nt_pos-3;
	@ref_frame_angle = &update_ref_frame($ref_frame_position);
	
	
	foreach $ATOM ( @ATOM ) {			#Copy the backbone of resi 4 of this motif,  o2' has to come at the end, so we pop it off at the end.
		if ( ($ATOM->{h} eq 'E')  &&  ($ATOM->{i} == 4) )	{							
			$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

			$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
			$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

			@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
			@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
			@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			


			@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
			@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

			@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

			$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

			push @PDB, {
				"a" => $ATOM->{a},  #atom type
				"n" => $seq[$nt_pos-1],  #residue type
				"i" => $nt_pos,  #residue number
				"h" => "A", #chain name		
				"x" => sprintf("%8s",$roundedx),
				"y" => sprintf("%8s",$roundedy),
				"z" => sprintf("%8s",$roundedz),					
			}			
		}		
	}
	$nt_temp = pop @PDB;		#remove the O2' atom
	#now, the reference frame from residue 4 is already calculated from before and step back 1 nt_step to place the crossover base.
	
	foreach $ATOM ( @ATOM ) {						#copy only the base atoms from the lib
		if ( ( ($ATOM->{h} eq 'A')  &&  ($ATOM->{n} eq $seq[$nt_pos-1]) &&
			( ($ATOM->{a} eq 'N9 ') || ($ATOM->{a} eq 'C8 ') ||
			  ($ATOM->{a} eq 'N7 ') || ($ATOM->{a} eq 'C5 ') ||
			  ($ATOM->{a} eq 'C6 ') || ($ATOM->{a} eq 'N6 ') ||
			  ($ATOM->{a} eq 'N1 ') || ($ATOM->{a} eq 'C2 ') ||
			  ($ATOM->{a} eq 'N3 ') || ($ATOM->{a} eq 'C4 ') ||
			  ($ATOM->{a} eq 'O6 ') || ($ATOM->{a} eq 'N2 ') ||
			  ($ATOM->{a} eq 'O2 ') || ($ATOM->{a} eq 'N4 ') || ($ATOM->{a} eq 'O4 ')) ) ) {  
			   						
			$moved[0]=$ATOM->{x};	$moved[1]=$ATOM->{y};	$moved[2]=$ATOM->{z};	$moved[3]=1;	

			@moved = &rotate_z ("@moved",-$base_ref_frame_angle[2]);  #first rotation    -first move the new nt to be in the frame of the last nt of motif
			@moved = &rotate_y ("@moved",-$base_ref_frame_angle[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$base_ref_frame_angle[0]); #final rotation	
			@moved = &translate_matrix ("@moved","@base_ref_trans");  #      Now moving to add to the 3prime end 
						
			$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

			push @PDB, {
				"a" => $ATOM->{a},  #atom type
				"n" => $seq[$nt_pos-1],  #residue type
				"i" => $nt_pos,  #residue number
				"h" => "A", #chain name		
				"x" => sprintf("%8s",$roundedx),
				"y" => sprintf("%8s",$roundedy),
				"z" => sprintf("%8s",$roundedz),											
			}			
		}
	}		
	push @PDB, $nt_temp;	#add the O2' atom back at the end
	$nt_pos += 1;
	
}

sub add_KL {
	$ref_frame_position = $nt_pos;
	@ref_frame_angle = &update_ref_frame($ref_frame_position);
	foreach $ATOM ( @ATOM ) {			#Copy the backbone of resi 1 of this motif,  o2' has to come at the end, so we pop it off at the end.
		if ( ($ATOM->{h} eq 'K')  &&  ($ATOM->{i} == 1) )	{							
			$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

			$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
			$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

			@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
			@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
			@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

			@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    
			@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

			@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

			$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

			push @PDB, {
				"a" => $ATOM->{a},  #atom type
				"n" => $seq[$nt_pos-1],  #residue type
				"i" => $nt_pos,  #residue number
				"h" => "A", #chain name		
				"x" => sprintf("%8s",$roundedx),
				"y" => sprintf("%8s",$roundedy),
				"z" => sprintf("%8s",$roundedz),															
			}			
		}		
	}	
	my $nt_temp = pop @PDB;		#remove the O2' atom

	foreach $ATOM ( @ATOM ) {						#copy only the base atoms from the lib
		if ( ( ($ATOM->{h} eq 'A')  &&  ($ATOM->{n} eq $seq[$nt_pos-1]) &&
			( ($ATOM->{a} eq 'N9 ') || ($ATOM->{a} eq 'C8 ') ||
			  ($ATOM->{a} eq 'N7 ') || ($ATOM->{a} eq 'C5 ') ||
			  ($ATOM->{a} eq 'C6 ') || ($ATOM->{a} eq 'N6 ') ||
			  ($ATOM->{a} eq 'N1 ') || ($ATOM->{a} eq 'C2 ') ||
			  ($ATOM->{a} eq 'N3 ') || ($ATOM->{a} eq 'C4 ') ||
			  ($ATOM->{a} eq 'O6 ') || ($ATOM->{a} eq 'N2 ') ||
			  ($ATOM->{a} eq 'O2 ') || ($ATOM->{a} eq 'N4 ') || ($ATOM->{a} eq 'O4 ') )  ) ) {   	
			  
			$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
			$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

			$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

			@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
			@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
			@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

			@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
			@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

			@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

			$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

			push @PDB, {
				"a" => $ATOM->{a},  #atom type
				"n" => $seq[$nt_pos-1],  #residue type
				"i" => $nt_pos,  #residue number
				"h" => "A", #chain name		
				"x" => sprintf("%8s",$roundedx),
				"y" => sprintf("%8s",$roundedy),
				"z" => sprintf("%8s",$roundedz),																
			}			
		}
	}		
	push @PDB, $nt_temp;	#add the O2' atom back at the end
	
	for ($i=2;$i<4;$i++) {    #this puts the -AA- bulge in
		foreach $ATOM ( @ATOM ) {			#we scan the datafile for the nt type and paste it in
			if ( ($ATOM->{h} eq 'K') &&  ($ATOM->{i} == $i)) {   #$c is the current resi, atom(i) is the resi number element in the PDB file.							
				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"n" => $ATOM->{n},  #residue type
					"i" => $nt_pos+($i-1),  #residue number
					"h" => "A", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),																	
				}			
			}
		}	
	}
	
	for ($i=4; $i<10; $i++){	#now it puts the programmable part in
		$ref_frame_position = $nt_pos;
		@ref_frame_angle = &update_ref_frame($ref_frame_position);
		
		foreach $ATOM ( @ATOM ) {			#Copy the backbone of resi 4 of this motif,  o2' has to come at the end, so we pop it off at the end.
			if ( ($ATOM->{h} eq 'K')  &&  ($ATOM->{i} == $i) )	{							

				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"n" => $seq[$nt_pos+$i-1-1],  #residue type
					"i" => $nt_pos+$i-1,  #residue number
					"h" => "A", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),																	
				}			
			}		
		}	
		$nt_temp = pop @PDB;		#remove the O2' atom
	
		$ref_frame_position = $nt_pos+$i;
		@ref_frame_angle = &update_ref_frame($ref_frame_position);
	
		foreach $ATOM ( @ATOM ) {						#copy only the base atoms from the lib
			if (  ($ATOM->{h} eq 'A')  &&  ($ATOM->{n} eq $seq[$nt_pos+$i-1-1]) &&
				( ($ATOM->{a} eq 'N9 ') || ($ATOM->{a} eq 'C8 ') ||
				  ($ATOM->{a} eq 'N7 ') || ($ATOM->{a} eq 'C5 ') ||
				  ($ATOM->{a} eq 'C6 ') || ($ATOM->{a} eq 'N6 ') ||
				  ($ATOM->{a} eq 'N1 ') || ($ATOM->{a} eq 'C2 ') ||
				  ($ATOM->{a} eq 'N3 ') || ($ATOM->{a} eq 'C4 ') ||
				  ($ATOM->{a} eq 'O6 ') || ($ATOM->{a} eq 'N2 ') ||
				  ($ATOM->{a} eq 'O2 ') || ($ATOM->{a} eq 'N4 ') || ($ATOM->{a} eq 'O4 ') )  ) {   						
					$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

					@moved = &rotate_z ("@point",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
					@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
					@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

					@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

					$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

					push @PDB, {
						"a" => $ATOM->{a},  #atom type
						"n" => $seq[$nt_pos+$i-1-1],  #residue type
						"i" => $nt_pos+$i-1,  #residue number
						"h" => "A", #chain name		
						"x" => sprintf("%8s",$roundedx),
						"y" => sprintf("%8s",$roundedy),
						"z" => sprintf("%8s",$roundedz),																	
					}			
			}
		}		
		push @PDB, $nt_temp;	#add the O2' atom back at the end
	}

	$ref_frame_position = $nt_pos;
	@ref_frame_angle = &update_ref_frame($ref_frame_position);
	
	for ($i=10;$i<11;$i++) {  
		foreach $ATOM ( @ATOM ) {			#we scan the datafile for the nt type and paste it in
			if ( ($ATOM->{h} eq 'K') &&  ($ATOM->{i} == $i)) {   #$c is the current resi, atom(i) is the resi number element in the PDB file.							
				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"n" => $ATOM->{n},  #residue type
					"i" => $nt_pos-1+$i,  #residue number
					"h" => "A", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),																	
				}			
			}
		}	
	}	
	$nt_pos += 10;

	foreach $ATOM ( @ATOM ) {			#Copy the backbone of the last resi of this motif,  o2' has to come at the end, so we pop it off at the end.
		if ( ($ATOM->{h} eq 'K')  &&  ($ATOM->{i} == 11) )	{							

			$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
			$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

			$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

			@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
			@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
			@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

			@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
			@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

			@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

			$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

			push @PDB, {
				"a" => $ATOM->{a},  #atom type
				"n" => $seq[$nt_pos-1],  #residue type
				"i" => $nt_pos,  #residue number
				"h" => "A", #chain name		
				"x" => sprintf("%8s",$roundedx),
				"y" => sprintf("%8s",$roundedy),
				"z" => sprintf("%8s",$roundedz),																	
			}			
		}		
	}	
	$nt_temp = pop @PDB;		#remove the O2' atom	
	
	$ref_frame_position = $nt_pos+1;
	@ref_frame_angle = &update_ref_frame($ref_frame_position);
		
	foreach $ATOM ( @ATOM ) {						#copy only the base atoms from the lib
		if ( ( ($ATOM->{h} eq 'A')  &&  ($ATOM->{n} eq $seq[$nt_pos-1]) &&
			( ($ATOM->{a} eq 'N9 ') || ($ATOM->{a} eq 'C8 ') ||
			  ($ATOM->{a} eq 'N7 ') || ($ATOM->{a} eq 'C5 ') ||
			  ($ATOM->{a} eq 'C6 ') || ($ATOM->{a} eq 'N6 ') ||
			  ($ATOM->{a} eq 'N1 ') || ($ATOM->{a} eq 'C2 ') ||
			  ($ATOM->{a} eq 'N3 ') || ($ATOM->{a} eq 'C4 ') ||
			  ($ATOM->{a} eq 'O6 ') || ($ATOM->{a} eq 'N2 ') ||
			  ($ATOM->{a} eq 'O2 ') || ($ATOM->{a} eq 'N4 ') || ($ATOM->{a} eq 'O4 ')) ) ) {   						
				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"n" => $seq[$nt_pos-1],  #residue type
					"i" => $nt_pos,  #residue number
					"h" => "A", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),																	
				}			
		}
	}		
	push @PDB, $nt_temp;	#add the O2' atom back at the end
	$nt_pos += 1;	
}



sub add_KLbend {
	$ref_frame_position = $nt_pos;
	@ref_frame_angle = &update_ref_frame($ref_frame_position);
	foreach $ATOM ( @ATOM ) {			#Copy the backbone of resi 1 of this motif,  o2' has to come at the end, so we pop it off at the end.
		if ( ($ATOM->{h} eq 'L')  &&  ($ATOM->{i} == 1) )	{							
			$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

			$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
			$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

			@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
			@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
			@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

			@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
			@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

			@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

			$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

			push @PDB, {
				"a" => $ATOM->{a},  #atom type
				"n" => $seq[$nt_pos-1],  #residue type
				"i" => $nt_pos,  #residue number
				"h" => "A", #chain name		
				"x" => sprintf("%8s",$roundedx),
				"y" => sprintf("%8s",$roundedy),
				"z" => sprintf("%8s",$roundedz),															
			}			
		}		
	}	
	my $nt_temp = pop @PDB;		#remove the O2' atom

	foreach $ATOM ( @ATOM ) {						#copy only the base atoms from the lib
		if ( ($ATOM->{h} eq 'A')  &&  ($ATOM->{n} eq $seq[$nt_pos-1]) &&
			( ($ATOM->{a} eq 'N9 ') || ($ATOM->{a} eq 'C8 ') ||
			  ($ATOM->{a} eq 'N7 ') || ($ATOM->{a} eq 'C5 ') ||
			  ($ATOM->{a} eq 'C6 ') || ($ATOM->{a} eq 'N6 ') ||
			  ($ATOM->{a} eq 'N1 ') || ($ATOM->{a} eq 'C2 ') ||
			  ($ATOM->{a} eq 'N3 ') || ($ATOM->{a} eq 'C4 ') ||
			  ($ATOM->{a} eq 'O6 ') || ($ATOM->{a} eq 'N2 ') ||
			  ($ATOM->{a} eq 'O2 ') || ($ATOM->{a} eq 'N4 ') || ($ATOM->{a} eq 'O4 ') ) ) {   	
			  
			$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
			$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

			$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

			@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
			@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
			@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

			@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
			@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

			@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

			$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

			push @PDB, {
				"a" => $ATOM->{a},  #atom type
				"n" => $seq[$nt_pos-1],  #residue type
				"i" => $nt_pos,  #residue number
				"h" => "A", #chain name		
				"x" => sprintf("%8s",$roundedx),
				"y" => sprintf("%8s",$roundedy),
				"z" => sprintf("%8s",$roundedz),																
			}			
		}
	}		
	push @PDB, $nt_temp;	#add the O2' atom back at the end
	for ($i=2;$i<9;$i++) {  
		foreach $ATOM ( @ATOM ) {			#we scan the datafile for the nt type and paste it in
			if ( ($ATOM->{h} eq 'L') &&  ($ATOM->{i} == $i)) {   #$c is the current resi, atom(i) is the resi number element in the PDB file.							
				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"n" => $ATOM->{n},  #residue type
					"i" => $nt_pos-1+$i,  #residue number
					"h" => "A", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),																	
				}			
			}
		}	
	}
	
	$nt_pos += 8;
		
	foreach $ATOM ( @ATOM ) {			#Copy the backbone of resi 9 of this motif,  o2' has to come at the end, so we pop it off at the end.
		if ( ($ATOM->{h} eq 'L')  &&  ($ATOM->{i} == 9) )	{							

			$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
			$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

			$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

			@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
			@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
			@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

			@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
			@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

			@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

			$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

			push @PDB, {
				"a" => $ATOM->{a},  #atom type
				"n" => $seq[$nt_pos-1],  #residue type
				"i" => $nt_pos,  #residue number
				"h" => "A", #chain name		
				"x" => sprintf("%8s",$roundedx),
				"y" => sprintf("%8s",$roundedy),
				"z" => sprintf("%8s",$roundedz),																	
			}			
		}		
	}	
	$nt_temp = pop @PDB;		#remove the O2' atom
	
	
	$ref_frame_position = $nt_pos+1;
	@ref_frame_angle = &update_ref_frame($ref_frame_position);
	
	
	foreach $ATOM ( @ATOM ) {						#copy only the base atoms from the lib
		if ( ( ($ATOM->{h} eq 'A')  &&  ($ATOM->{n} eq $seq[$nt_pos-1]) &&
			( ($ATOM->{a} eq 'N9 ') || ($ATOM->{a} eq 'C8 ') ||
			  ($ATOM->{a} eq 'N7 ') || ($ATOM->{a} eq 'C5 ') ||
			  ($ATOM->{a} eq 'C6 ') || ($ATOM->{a} eq 'N6 ') ||
			  ($ATOM->{a} eq 'N1 ') || ($ATOM->{a} eq 'C2 ') ||
			  ($ATOM->{a} eq 'N3 ') || ($ATOM->{a} eq 'C4 ') ||
			  ($ATOM->{a} eq 'O6 ') || ($ATOM->{a} eq 'N2 ') ||
			  ($ATOM->{a} eq 'O2 ') || ($ATOM->{a} eq 'N4 ') || ($ATOM->{a} eq 'O4 ')) ) ) {   						
				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"n" => $seq[$nt_pos-1],  #residue type
					"i" => $nt_pos,  #residue number
					"h" => "A", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),																	
				}			
		}
	}		
	push @PDB, $nt_temp;	#add the O2' atom back at the end
	
	$nt_pos += 1;
}



##################

sub add_tar {
	$ref_frame_position = $nt_pos;
	@ref_frame_angle = &update_ref_frame($ref_frame_position);
	foreach $ATOM ( @ATOM ) {			#Copy the backbone of resi 1 of this motif,  o2' has to come at the end, so we pop it off at the end.
		if ( ($ATOM->{h} eq 'R')  &&  ($ATOM->{i} == 1) )	{							
			$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

			$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
			$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

			@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
			@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
			@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

			@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
			@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

			@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

			$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

			push @PDB, {
				"a" => $ATOM->{a},  #atom type
				"n" => $seq[$nt_pos-1],  #residue type
				"i" => $nt_pos,  #residue number
				"h" => "A", #chain name		
				"x" => sprintf("%8s",$roundedx),
				"y" => sprintf("%8s",$roundedy),
				"z" => sprintf("%8s",$roundedz),															
			}			
		}		
	}	
	my $nt_temp = pop @PDB;		#remove the O2' atom

	foreach $ATOM ( @ATOM ) {						#copy only the base atoms from the lib
		if ( ( ($ATOM->{h} eq 'A')  &&  ($ATOM->{n} eq $seq[$nt_pos-1]) &&
			( ($ATOM->{a} eq 'N9 ') || ($ATOM->{a} eq 'C8 ') ||
			  ($ATOM->{a} eq 'N7 ') || ($ATOM->{a} eq 'C5 ') ||
			  ($ATOM->{a} eq 'C6 ') || ($ATOM->{a} eq 'N6 ') ||
			  ($ATOM->{a} eq 'N1 ') || ($ATOM->{a} eq 'C2 ') ||
			  ($ATOM->{a} eq 'N3 ') || ($ATOM->{a} eq 'C4 ') ||
			  ($ATOM->{a} eq 'O6 ') || ($ATOM->{a} eq 'N2 ') ||
			  ($ATOM->{a} eq 'O2 ') || ($ATOM->{a} eq 'N4 ') || ($ATOM->{a} eq 'O4 ')) ) ) {   	
			  
			$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
			$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

			$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

			@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
			@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
			@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

			@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
			@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

			@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

			$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

			push @PDB, {
				"a" => $ATOM->{a},  #atom type
				"n" => $seq[$nt_pos-1],  #residue type
				"i" => $nt_pos,  #residue number
				"h" => "A", #chain name		
				"x" => sprintf("%8s",$roundedx),
				"y" => sprintf("%8s",$roundedy),
				"z" => sprintf("%8s",$roundedz),																
			}			
		}
	}		
	push @PDB, $nt_temp;	#add the O2' atom back at the end
	for ($i=2;$i<28;$i++) {  
		foreach $ATOM ( @ATOM ) {			#we scan the datafile for the nt type and paste it in
			if ( ($ATOM->{h} eq 'R') &&  ($ATOM->{i} == $i)) {   #$c is the current resi, atom(i) is the resi number element in the PDB file.							
				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"n" => $ATOM->{n},  #residue type
					"i" => $nt_pos-1+$i,  #residue number
					"h" => "A", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),																	
				}			
			}
		}	
	}

	for ($i=100;$i<115;$i++) {  
		foreach $ATOM ( @ATOM ) {			#we scan the datafile for the nt type and paste it in
			if ( ($ATOM->{h} eq 'R') &&  ($ATOM->{i} == $i)) {   #$c is the current resi, atom(i) is the resi number element in the PDB file.							
				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				if ($i==113) {	$end_of_chain=1;} else{$end_of_chain=0;} # this gets set to 1 at the end of the chain

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"p" => $ATOM->{p},  #residue type
					"i" => $nt_pos-1+$i,  #residue number
					"h" => "P", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),
					"t" => $end_of_chain,																	
				} #end_push			
			}
		}	
	}


	
	$nt_pos += 27;
		
	foreach $ATOM ( @ATOM ) {			#Copy the backbone of resi 1 of this motif,  o2' has to come at the end, so we pop it off at the end.
		if ( ($ATOM->{h} eq 'R')  &&  ($ATOM->{i} == 28) )	{							

			$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
			$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

			$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

			@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
			@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
			@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

			@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
			@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
			@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

			@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

			$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

			push @PDB, {
				"a" => $ATOM->{a},  #atom type
				"n" => $seq[$nt_pos-1],  #residue type
				"i" => $nt_pos,  #residue number
				"h" => "A", #chain name		
				"x" => sprintf("%8s",$roundedx),
				"y" => sprintf("%8s",$roundedy),
				"z" => sprintf("%8s",$roundedz),																	
			}			
		}		
	}	
	$nt_temp = pop @PDB;		#remove the O2' atom
	
	
	$ref_frame_position = $nt_pos+1;
	@ref_frame_angle = &update_ref_frame($ref_frame_position);
	
	
	foreach $ATOM ( @ATOM ) {						#copy only the base atoms from the lib
		if ( ( ($ATOM->{h} eq 'A')  &&  ($ATOM->{n} eq $seq[$nt_pos-1]) &&
			( ($ATOM->{a} eq 'N9 ') || ($ATOM->{a} eq 'C8 ') ||
			  ($ATOM->{a} eq 'N7 ') || ($ATOM->{a} eq 'C5 ') ||
			  ($ATOM->{a} eq 'C6 ') || ($ATOM->{a} eq 'N6 ') ||
			  ($ATOM->{a} eq 'N1 ') || ($ATOM->{a} eq 'C2 ') ||
			  ($ATOM->{a} eq 'N3 ') || ($ATOM->{a} eq 'C4 ') ||
			  ($ATOM->{a} eq 'O6 ') || ($ATOM->{a} eq 'N2 ') ||
			  ($ATOM->{a} eq 'O2 ') || ($ATOM->{a} eq 'N4 ') || ($ATOM->{a} eq 'O4 ')) ) ) {   						
				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"n" => $seq[$nt_pos-1],  #residue type
					"i" => $nt_pos,  #residue number
					"h" => "A", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),																	
				}			
		}
	}		
	push @PDB, $nt_temp;	#add the O2' atom back at the end
	
	$nt_pos += 1;
}

sub add_pp7 {
	$ref_frame_position = $nt_pos;
	@ref_frame_angle = &update_ref_frame($ref_frame_position);

	for ($i=1;$i<26;$i++) {  
		foreach $ATOM ( @ATOM ) {			#we scan the datafile for the nt type and paste it in
			if ( ($ATOM->{h} eq 'N') &&  ($ATOM->{i} == $i)) {   #$c is the current resi, atom(i) is the resi number element in the PDB file.							
				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"n" => $ATOM->{n},  #residue type
					"i" => $nt_pos-1+$i,  #residue number
					"h" => "A", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),																	
				}			
			}
		}	
	}

	for ($i=100;$i<344;$i++) {  
		foreach $ATOM ( @ATOM ) {			#we scan the datafile for the nt type and paste it in
			if ( ($ATOM->{h} eq 'N') &&  ($ATOM->{i} == $i)) {   #$c is the current resi, atom(i) is the resi number element in the PDB file.							
				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				if ($i==343) {	$end_of_chain=1;} else{$end_of_chain=0;} # this gets set to 1 at the end of the chain

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"p" => $ATOM->{p},  #residue type
					"i" => $nt_pos-1+$i,  #residue number
					"h" => "P", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),
					"t" => $end_of_chain,
																						
				}			
			}
		}	
	}
	

	
	$nt_pos += 25;
		
}


sub add_ms2 {
	$ref_frame_position = $nt_pos;
	@ref_frame_angle = &update_ref_frame($ref_frame_position);

	for ($i=1;$i<14;$i++) {  
		foreach $ATOM ( @ATOM ) {			#we scan the datafile for the nt type and paste it in
			if ( ($ATOM->{h} eq 'M') &&  ($ATOM->{i} == $i)) {   #$c is the current resi, atom(i) is the resi number element in the PDB file.							
				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"n" => $ATOM->{n},  #residue type
					"i" => $nt_pos-1+$i,  #residue number
					"h" => "A", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),																	
				}			
			}
		}	
	}

	for ($i=100;$i<358;$i++) {  
		foreach $ATOM ( @ATOM ) {			#we scan the datafile for the nt type and paste it in
			if ( ($ATOM->{h} eq 'M') &&  ($ATOM->{i} == $i)) {   #$c is the current resi, atom(i) is the resi number element in the PDB file.							
				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				if ($i==356) {	$end_of_chain=1;} else{$end_of_chain=0;} # this gets set to 1 at the end of the chain

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"p" => $ATOM->{p},  #residue type
					"i" => $nt_pos-1+$i,  #residue number
					"h" => "P", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),
					"t" => $end_of_chain,																	
				}			
			}
		}	
	}
	
	$nt_pos += 13;
}

sub add_mango {
	$ref_frame_position = $nt_pos;
	@ref_frame_angle = &update_ref_frame($ref_frame_position);

	for ($i=1;$i<30;$i++) {  
		foreach $ATOM ( @ATOM ) {			#we scan the datafile for the nt type and paste it in
			if ( ($ATOM->{h} eq 'Q') &&  ($ATOM->{i} == $i)) {   #$c is the current resi, atom(i) is the resi number element in the PDB file.							
				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"n" => $ATOM->{n},  #residue type
					"i" => $nt_pos-1+$i,  #residue number
					"h" => "A", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),																	
				}			
			}
		}	
	}

	for ($i=100;$i<101;$i++) {  
		foreach $ATOM ( @ATOM ) {			#we scan the datafile for the nt type and paste it in
			if ( ($ATOM->{h} eq 'Q') &&  ($ATOM->{i} == $i)) {   #$c is the current resi, atom(i) is the resi number element in the PDB file.							
				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"p" => $ATOM->{p},  #residue type
					"i" => $nt_pos-1+$i,  #residue number
					"h" => "P", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),	
					"t" => "1",																
				}			
			}
		}	
	}
	
	$nt_pos += 29;		

}

sub add_spinacha {
	$ref_frame_position = $nt_pos;
	@ref_frame_angle = &update_ref_frame($ref_frame_position);
	
	for ($i=100;$i<122;$i++) {  
		foreach $ATOM ( @ATOM ) {			#we scan the datafile for the nt type and paste it in
			if ( ($ATOM->{h} eq 'S') &&  ($ATOM->{i} == $i)) {   #$c is the current resi, atom(i) is the resi number element in the PDB file.							
				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"n" => $ATOM->{n},  #residue type
					"i" => $nt_pos-1+$i-99,  #residue number
					"h" => "A", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),																	
				}			
			}
		}	
	}
	
	$nt_pos += 22;		
}

sub add_spinachb {
	$ref_frame_position = $nt_pos;
	@ref_frame_angle = &update_ref_frame($ref_frame_position);

	for ($i=1;$i<18;$i++) {  
		foreach $ATOM ( @ATOM ) {			#we scan the datafile for the nt type and paste it in
			if ( ($ATOM->{h} eq 'S') &&  ($ATOM->{i} == $i)) {   #$c is the current resi, atom(i) is the resi number element in the PDB file.							
				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"n" => $ATOM->{n},  #residue type
					"i" => $nt_pos-1+$i,  #residue number
					"h" => "A", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),																	
				}			
			}
		}	
	}
	
	for ($i=200;$i<201;$i++) {  
		foreach $ATOM ( @ATOM ) {			#we scan the datafile for the nt type and paste it in
			if ( ($ATOM->{h} eq 'S') &&  ($ATOM->{i} == $i)) {   #$c is the current resi, atom(i) is the resi number element in the PDB file.							
				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"p" => $ATOM->{p},  #residue type
					"i" => $nt_pos-1+$i,  #residue number
					"h" => "P", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),
					"t" => "1",																		
				}			
			}
		}	
	}

	
	$nt_pos += 17;
		
}


sub add_kturna {
	$ref_frame_position = $nt_pos;
	@ref_frame_angle = &update_ref_frame($ref_frame_position);
	
	for ($i=1;$i<7;$i++) {  
		foreach $ATOM ( @ATOM ) {			#we scan the datafile for the nt type and paste it in
			if ( ($ATOM->{h} eq 'B') &&  ($ATOM->{i} == $i)) {   #$c is the current resi, atom(i) is the resi number element in the PDB file.							
				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"n" => $ATOM->{n},  #residue type
					"i" => $nt_pos-1+$i,  #residue number
					"h" => "A", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),																	
				}			
			}
		}	
	}

	for ($i=200;$i<313;$i++) {  
		foreach $ATOM ( @ATOM ) {			#we scan the datafile for the nt type and paste it in
			if ( ($ATOM->{h} eq 'B') &&  ($ATOM->{i} == $i)) {   #$c is the current resi, atom(i) is the resi number element in the PDB file.							
				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				if ($i==113) {	$end_of_chain=1;} else{$end_of_chain=0;} # this gets set to 1 at the end of the chain

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"p" => $ATOM->{p},  #residue type
					"i" => $nt_pos-1+$i,  #residue number
					"h" => "P", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),	
					"t" => $end_of_chain,													
				}			
			}
		}	
	}
	
	$nt_pos += 6;		
}

sub add_kturnb {
	$ref_frame_position = $nt_pos;
	@ref_frame_angle = &update_ref_frame($ref_frame_position);

	for ($i=101;$i<110;$i++) {  
		foreach $ATOM ( @ATOM ) {			#we scan the datafile for the nt type and paste it in
			if ( ($ATOM->{h} eq 'B') &&  ($ATOM->{i} == $i)) {   #$c is the current resi, atom(i) is the resi number element in the PDB file.							
				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"n" => $ATOM->{n},  #residue type
					"i" => $nt_pos-1+$i-100,  #residue number
					"h" => "A", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),																	
				}			
			}
		}	
	}
	
	$nt_pos += 9;
		
}

####################################################


sub add_bkla {
	$ref_frame_position = $nt_pos;
	@ref_frame_angle = &update_ref_frame($ref_frame_position);
	
	for ($i=1;$i<12;$i++) {  
		foreach $ATOM ( @ATOM ) {			#we scan the datafile for the nt type and paste it in
			if ( ($ATOM->{h} eq 'V') &&  ($ATOM->{i} == $i)) {   #$c is the current resi, atom(i) is the resi number element in the PDB file.							
				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"n" => $ATOM->{n},  #residue type
					"i" => $nt_pos-1+$i,  #residue number
					"h" => "A", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),																	
				}			
			}
		}	
	}
	
	$nt_pos += 11;		
}

sub add_bklb {
	$ref_frame_position = $nt_pos;
	@ref_frame_angle = &update_ref_frame($ref_frame_position);

	for ($i=1;$i<5;$i++) {  
		foreach $ATOM ( @ATOM ) {			#we scan the datafile for the nt type and paste it in
			if ( ($ATOM->{h} eq 'W') &&  ($ATOM->{i} == $i)) {   #$c is the current resi, atom(i) is the resi number element in the PDB file.							
				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"n" => $ATOM->{n},  #residue type
					"i" => $nt_pos-1+$i,  #residue number
					"h" => "A", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),																	
				}			
			}
		}	
	}
	
	$nt_pos += 4;
		
}


sub add_ninetya {
	$ref_frame_position = $nt_pos;
	@ref_frame_angle = &update_ref_frame($ref_frame_position);
	
	for ($i=1;$i<10;$i++) {  
		foreach $ATOM ( @ATOM ) {			#we scan the datafile for the nt type and paste it in
			if ( ($ATOM->{h} eq 'B') &&  ($ATOM->{i} == $i)) {   #$c is the current resi, atom(i) is the resi number element in the PDB file.							
				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"n" => $ATOM->{n},  #residue type
					"i" => $nt_pos-1+$i,  #residue number
					"h" => "A", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),																	
				}			
			}
		}	
	}


	
	$nt_pos += 9;		
}

sub add_ninetyb {
	$ref_frame_position = $nt_pos;
	@ref_frame_angle = &update_ref_frame($ref_frame_position);

	for ($i=1;$i<5;$i++) {  
		foreach $ATOM ( @ATOM ) {			#we scan the datafile for the nt type and paste it in
			if ( ($ATOM->{h} eq 'B') &&  ($ATOM->{i} == $i)) {   #$c is the current resi, atom(i) is the resi number element in the PDB file.							
				$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
				$angles[0]=-1.05152505;	$angles[1]=0.46918242;	$angles[2]=1.37954160;			

				$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

				@moved = &rotate_z ("@point",-$angles[2]);  #first rotation  -we move the nt to the untranslocated pos
				@moved = &rotate_y ("@moved",-$angles[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$angles[0]); #final rotation
				@moved = &translate_matrix ("@moved","@nt_step");  #POINT followed by translation vector			

				@moved = &rotate_z ("@moved",-$ref_frame_angle[2]);  #first rotation    -now we rotate aout 
				@moved = &rotate_y ("@moved",-$ref_frame_angle[1]); #second rotation	
				@moved = &rotate_z ("@moved",-$ref_frame_angle[0]); #final rotation	

				@moved = &translate_matrix ("@moved","@ref_trans");  #      Now moving to add to the 3prime end 

				$roundedx = sprintf("%4.3f", $moved[0]);	$roundedy = sprintf("%4.3f", $moved[1]);	$roundedz = sprintf("%4.3f", $moved[2]);					

				push @PDB, {
					"a" => $ATOM->{a},  #atom type
					"n" => $ATOM->{n},  #residue type
					"i" => $nt_pos-1+$i,  #residue number
					"h" => "A", #chain name		
					"x" => sprintf("%8s",$roundedx),
					"y" => sprintf("%8s",$roundedy),
					"z" => sprintf("%8s",$roundedz),																	
				}			
			}
		}	
	}
	
	$nt_pos += 4;
		
}









####################################

sub map {
	@map = ();
	my @secondary = split('', $_[0]);
	my @last_bracket;	my @last_brace;	my @last_curlybrace;  #markers to remember the 	
	my @last_a_brace;	my @last_b_brace;	my @last_c_brace;
	
	my $counter = 0;
	foreach(@secondary){
		if ($_ eq "."){
			$map[$counter] = $counter;	 ##Single strands and map to themselves.
		}
		elsif ($_ eq "("){
			push @last_bracket, $counter;		   
		}
		elsif ($_ eq ")"){
			my $partner = pop @last_bracket;
			$map[$counter]=$partner;
			$map[$partner]=$counter;
		}
		elsif ($_ eq "["){
			push @last_brace, $counter;	
			#$map[$counter] = $counter;	 ##Single strands and map to themselves.
	 
		}
		elsif ($_ eq "]"){
			my $buddy = pop @last_brace;
			$map[$counter]=$buddy;
			$map[$buddy]=$counter;
			#$map[$counter] = $counter;	 ##Single strands and map to themselves.
		}
		$counter ++;
	}
	
	######## Scrub the sequence and fill all non-AUCG space with a random N
	
	$counter=0;
	foreach(@seq){
		if($seq[$counter] eq 'A' || $seq[$counter] eq 'U' || $seq[$counter] eq 'C' || $seq[$counter] eq 'G'){
			
		} else {
			my $pick = int rand(4);			
			if ($pick == 0){$seq[$counter]= 'A'; $seq[$map[$counter]]='U'}
			if ($pick == 1){$seq[$counter]= 'U'; $seq[$map[$counter]]='A'}
			if ($pick == 2){$seq[$counter]= 'G'; $seq[$map[$counter]]='C'}
			if ($pick == 3){$seq[$counter]= 'C'; $seq[$map[$counter]]='G'}
		}
		$counter ++;
	}
	
}

#####################################

sub trace{
	my ( $line, @cols, $seq );
	my @pri = ( );

	@lines = &read_file( $file1 ); #use the read_file sub to read in the input file


	# translate unicode to symbols:
	# ----------------------------------------
	# symbol unicode       unique description
	# ----------------------------------------
	# _      \302\240      240    space       
	# -      \342\224\200  200    straight    
	# i      \342\224\202  202    up-down
	# p      \342\224\212  212    pair        
	# x      \342\224\274  274    cross       
	# L      \342\225\255  255    down-right
	# J      \342\225\256  256    down-left
	# 7      \342\225\257  257    up-left
	# r      \342\225\260  260    up-right
	# b                           base-pair vertical
	# ----------------------------------------

	my @m = ( );	my @n = ( );	my @t = ( );	@cols = ( );
	my $i = 0;	my $j = 0;	my $cols = 0;	my $l = 0;
	my $k = 0;	my $maxlength = 0;	my $endpadding = 0;	my $numtees = 0;
	my $wide =0;
	my $widest = 0;
	my @cross = ();
	my $KL_pattern = '';



	foreach $line ( @lines ) {
		$l = length $line;
		
		$k = 0;
		$wide=0;
		
		if ($j<1){
			if($l<500){	$maxlength = 500;}else{$maxlength=$l;}  #make the width at least 500, unless the user specifies something longer in the first line
		}


		for ($i=0;$i<$l;$i++) {
			$cols = substr("$line", $i, 1);
			if ( $cols eq "\t" ) { die "Error: Pattern file contains tabs, please remove them and replace with spaces.\n";} 			
			if ( $cols eq "#" ) { $i=$i+$maxlength; }
			if ( $cols eq ">" ) { $name = substr("$line", $i+1, 32); $i=$i+$maxlength; } #Extract the File Name Here with a max 32 chars, skip to the next line
			if ( $cols eq "@" ) { $KL_pattern = substr("$line", $i+1, 32);  $i=$i+$maxlength; } #Extract the KL matching pattern with a max 32 chars, skip to the next line

			if ( $k==0 ) { $m[$j][$k] = " "; $k++; } #add a space at the left edge
			
			if ( $cols eq " " ) { $m[$j][$k] = " "; $k++; }     
			if ( $cols eq "*" ) { $m[$j][$k] = "*"; $k++; $wide=$k;}        
			if ( $cols eq "\257" ) { $m[$j][$k] = "J"; $k++; $wide=$k;}
			if ( $cols eq "\/" &&  substr("$line", $i+1, 1)!~ /[+xNXGACURYKMSWVHBDT35-]/ ) { $m[$j][$k] = "J"; $k++; $wide=$k;}     
			if ( $cols eq "\200" ) { $m[$j][$k] = "-"; $k++; $wide=$k;}
			if ( $cols eq "-" ) { $m[$j][$k] = "-"; $k++; $wide=$k;}     
			if ( $cols eq "\260" ) { $m[$j][$k] = "L"; $k++; $wide=$k;}
			if ( $cols eq "\\" &&  substr("$line", $i+1, 1)=~ /[+xNXGACURYKMSWVHBDT35-]/ ) { $m[$j][$k] = "L"; $k++; $wide=$k;}      
			if ( $cols eq "\255" ) { $m[$j][$k] = "r"; $k++; $wide=$k;}
			if ( $cols eq "\/" &&  substr("$line", $i+1, 1)=~ /[+xNXGACURYKMSWVHBDT35-]/ ) { $m[$j][$k] = "r"; $k++; $wide=$k;}     
			if ( $cols eq "\256" ) { $m[$j][$k] = "7"; $k++; $wide=$k;} 
			if ( $cols eq "\\" &&  substr("$line", $i+1, 1)!~ /[+xNXGACURYKMSWVHBDT35-]/ ) { $m[$j][$k] = "7"; $k++; $wide=$k;} 
			
			if ( $cols eq "\202" ) { $m[$j][$k] = "i"; $k++; $wide=$k;}
			if ( $cols eq "|" ) { $m[$j][$k] = "i"; $k++; $wide=$k;}   
			if ( $cols eq "^" ) { $m[$j][$k] = "^"; $cross[$j][$k] = "^"; $k++;$wide=$k;}  
			   
			if ( $cols eq "\212" ) { $m[$j][$k] = "p"; $k++; $wide=$k;} 
			if ( $cols eq ":" ) { $m[$j][$k] = "p"; $k++; $wide=$k;} 
			if ( $cols eq "\;" ) { $m[$j][$k] = "p"; $k++; $wide=$k;}  
			  
			if ( $cols eq "!"  ) { $m[$j][$k] = "!"; $k++; $wide=$k;}   
			       
			if ( $cols eq "\274" ) { $m[$j][$k] = "x"; $k++; $wide=$k;}
			if ( $cols eq "+" ) { $m[$j][$k] = "x"; $k++; $wide=$k;}
			
			if ( $cols eq "="    ) { $m[$j][$k] = "b"; $k++; $wide=$k;}   
			  
			if ( $cols =~ /\w/ ) { $m[$j][$k] = "$cols"; $k++; $wide=$k;}
			
			if ($k>$maxlength){  #if we are at the new longest row, then go back and add space buffers to the ends of previous rows and update the maxlength.
				for (my $mi=0;$mi<($j+1); $mi++){
					for (my $pn = $maxlength; $pn < ($k); $pn++) { 
						if (defined $m[$mi][$pn])  { 
							if ($m[$mi][$pn] eq "\n"){$m[$mi][$pn]= " "; $m[$mi][$pn+1]="\n";}  # the \n needs to be recognized so we can remove it and extend the row.
						} else { $m[$mi][$pn]= "U"; $m[$mi][$pn+1]="\n"; }
					}					
				}
				$maxlength = $k;
			}
			
			if ($wide>$widest){$widest =$wide;}
			
		}

		if ($k<$maxlength){$endpadding=$maxlength-$k;} #if this row is shorter than the max, then we need to fill out the buffers
		for ($i=0;$i<$endpadding;$i++) {
			$m[$j][$k] = " "; $k++;
		}
		$m[$j][$k] = "\n"; $k++;

		$j++;
	}
	my $tallest = $j;
	
	#trim off the excess whitespace on the right edge
	
	for ($i=0; $i<$l; $i++){
		for ($j=$widest+1; $j<$maxlength; $j++){
			$m[$i][$j] = "";
		}		
		$m[$widest][$j]= "\n";
	}	
	
	#add extra row of spaces at the end
	for ($i=0; $i<$widest+2; $i++) {
		$m[$tallest][$i] = " ";
		$m[$tallest+1][$i] = " ";
	}
	$m[$tallest][$widest+2]="\n";
	$m[$tallest+1][$widest+2]="\n";

	#scrub filename and KL-pattern inputs
	$find = " ";	#remove extra spaces at the end
	$replace = "";
	$name =~ s/$find/$replace/g;
	$find = "\n";	#remove carriage return if its in the file name
	$replace = "";
	$name =~ s/$find/$replace/g;
	$KL_pattern =~ s/$find/$replace/g;
	my @KL_pat = split(//,$KL_pattern);
	my $current_KL=0;




	# Find 5 prime end and assign directionality based on the sequence in adjacent elements 

	my $r = 0;
	my $c = 0;
	my $d = "left";
	for ($i=0;$i<1000;$i++) {
		for ($j=0;$j<1000;$j++) {
			if ( defined $m[$i][$j] ) { 
				if ($m[$i][$j] =~ /\d+/ ) {
					if ( scalar $m[$i][$j] && $m[$i][$j] == 5 ) { 
						$r = $i;
						$c = $j;
						if ( $m[$r][$c+1] =~ /[NXGACURYKMSWVHBDT-]/ ) { $d = "right"; }
						if ( $m[$r][$c-1] =~ /[NXGACURYKMSWVHBDT-]/ ) { $d = "left"; }
						if ( $m[$r+1][$c] =~ /[NXGACURYKMSWVHBDTi^]/ ) { $d = "down"; }
						if ( $m[$r-1][$c] =~ /[NXGACURYKMSWVHBDTi^]/ ) { $d = "up"; }  #note: in the case of multiple options it will pick one
																						#adding detection of invalid strand path could make this more foolproof
					}
				}
			}
		}
	}

	# Find total number of nucleotides in blueprint

	my $nt = 0;
	for ($i=0;$i<1000;$i++) {
		for ($j=0;$j<1000;$j++) {
			if ( defined $m[$i][$j] ) {
				if ($m[$i][$j] =~ /[NXGACURYKMSWVHBDT]/ ) {
					$nt++;
				}
			}
		} 
	}

	#Trace the strand

	my $r2 = scalar ($r);
	my $c2 = scalar ($c);
	my $d2 = $d;
	my $num = 0;
	my @seq = ( );
	my $test = 'test';

	for ($k=0;$k<$nt+10000;$k++) { 
		# TRACE HORIZONTAL STRAND

		if ( $d eq "right" && $m[$r][$c+1] =~ /[xNXGACURYKMSWVHBDT-]/ ) { 
			if ( $m[$r][$c+1] =~ /[NXGACURYKMSWVHBDT]/ ) { $num++; $n[$r][$c+1] = $num; push @seq, $m[$r][$c+1]; } 
			$c++;
		}
		if ( $d eq "left"  && $m[$r][$c-1] =~ /[xNXGACURYKMSWVHBDT-]/ ) { 
			if ( $m[$r][$c-1] =~ /[NXGACURYKMSWVHBDT]/ ) { $num++; $n[$r][$c-1] = $num; push @seq, $m[$r][$c-1]; }
			$c--;
		}

		if ( $d eq "up" && $m[$r-1][$c] =~ /[xNXGACURYKMSWVHBDTi^]/ ) { 
			if ( $m[$r-1][$c] =~ /[NXGACURYKMSWVHBDT^]/ ) { $num++; $n[$r-1][$c] = $num; push @seq, $m[$r-1][$c]; } 
			$r--; 
		}
		if ( $d eq "down" && $m[$r+1][$c] =~ /[xNXGACURYKMSWVHBDTi^]/ ) { 
			if ( $m[$r+1][$c] =~ /[NXGACURYKMSWVHBDT^]/ ) { $num++; $n[$r+1][$c] = $num; push @seq, $m[$r+1][$c]; } 
			$r++; 
		}
	
		# CROSS-OVER
		if ( $d eq "right" ) {
			if ( $m[$r][$c+1] eq "7" ) { $d = "down"; $c++;}
			if ( $m[$r][$c+1] eq "J" ) { $d = "up"; $c++;}
		}    
		if ( $d eq "left" ) {
			if ( $m[$r][$c-1] eq "L" ) { $d = "up"; $c--;}
			if ( $m[$r][$c-1] eq "r" ) { $d = "down";  $c--;}
		}
		if ( $d eq "down" ) {
			if ( $m[$r+1][$c] eq "L" ) { $d = "right"; $r++;}
			if ( $m[$r+1][$c] eq "J" ) { $d = "left"; $r++;}
		}    
		if ( $d eq "up" ) {
			if ( $m[$r-1][$c] eq "7" ) { $d = "left"; $r--;}
			if ( $m[$r-1][$c] eq "r" ) { $d = "right"; $r--;}
		}

		# FIND 3-PRIME END
		if ( $m[$r][$c+1] =~ /\d/ ) { if ( $m[$r][$c+1] == 3 ) { $test = "success"; last; } }
		if ( $m[$r][$c-1] =~ /\d/ ) { if ( $m[$r][$c-1] == 3 ) { $test = "success"; last; } }
		if ( $m[$r+1][$c] =~ /\d/ ) { if ( $m[$r+1][$c] == 3 ) { $test = "success"; last; } }
		if ( $m[$r-1][$c] =~ /\d/ ) { if ( $m[$r-1][$c] == 3 ) { $test = "success"; last; } }
	}

	# Find base pairs

	$r = scalar ($r2); # reset row
	$c = scalar ($c2); # reset column
	$d = $d2; # restore 5p direction
	my @a = ( );
	my @b = ( );
	my @p = ( );

	for ($k=0;$k<$nt+10000;$k++) {
		# TRACE HORIZONTAL STRAND
		if ( $d eq "right" && $m[$r][$c+1] =~ /[xNXGACURYKMSWVHBDT-]/ ) { 
			if ( $m[$r][$c+1] =~ /[NXGACURYKMSWVHBDT]/ ) { 
				if ( $m[$r+1][$c+1] =~ /[!p\*]/ ) {
					push @a, $n[$r][$c+1];
					push @b, $n[$r+2][$c+1]; 
					push @p, $m[$r+1][$c+1]; 
				} elsif ( $m[$r-1][$c+1] =~ /[!p\*]/ ) {
					push @a, $n[$r][$c+1];
					push @b, $n[$r-2][$c+1]; 
					push @p, $m[$r-1][$c+1]; 
				} else {
					push @a, $n[$r][$c+1];
					push @b, 0; 
					push @p, "-";
				}
			}
			$c++;
		}
		if ( $d eq "left"  && $m[$r][$c-1] =~ /[xNXGACURYKMSWVHBDT-]/ ) { 
			if ( $m[$r][$c-1] =~ /[NXGACURYKMSWVHBDT]/ ) { 
				if ( $m[$r+1][$c-1] =~ /[!p\*]/ ) {
					push @a, $n[$r][$c-1];
					push @b, $n[$r+2][$c-1]; 
					push @p, $m[$r+1][$c-1]; 
				} elsif ( $m[$r-1][$c-1] =~ /[!p\*]/ ) {
					push @a, $n[$r][$c-1];
					push @b, $n[$r-2][$c-1]; 
					push @p, $m[$r-1][$c-1]; 
				} else {
					push @a, $n[$r][$c-1];
					push @b, 0; 
					push @p, "-";
				}
			}
			$c--; 
		}

		if ( $d eq "down" && $m[$r+1][$c] =~ /[xNXGACURYKMSWVHBDTi^]/ ) { 
			if ( $m[$r+1][$c] =~ /[NXGACURYKMSWVHBDT^]/ ) { 
				if ( $m[$r+1][$c+1] =~ /[b\*]/ ) {
					push @a, $n[$r+1][$c];
					push @b, $n[$r+1][$c+2]; 
					push @p, $m[$r+1][$c+1];  
				} elsif ( $m[$r+1][$c-1] =~ /[b\*]/ ) {
					push @a, $n[$r+1][$c];
					push @b, $n[$r+1][$c-2]; 
					push @p, $m[$r+1][$c-1]; 
				} else {
					push @a, $n[$r+1][$c];
					push @b, 0; 
					push @p, "i";
				}
			}
			$r++;
		}

		if ( $d eq "up" && $m[$r-1][$c] =~ /[xNXGACURYKMSWVHBDTi^]/ ) { 
			if ( $m[$r-1][$c] =~ /[NXGACURYKMSWVHBDT^]/ ) { 
				if ( $m[$r-1][$c+1] =~ /[b\*]/ ) {
					push @a, $n[$r-1][$c];
					push @b, $n[$r-1][$c+2]; 
					push @p, $m[$r-1][$c+1]; 
				} elsif ( $m[$r-1][$c-1] =~ /[b\*]/ ) {
					push @a, $n[$r-1][$c];
					push @b, $n[$r-1][$c-2]; 
					push @p, $m[$r-1][$c-1]; 
				} else {
					push @a, $n[$r-1][$c];
					push @b, 0; 
					push @p, "i";
				}
			}
			$r--;
		}

		if ( $d eq "right" && $m[$r][$c+1] =~ /[x-]/ ) { $c++; }
		if ( $d eq "left"  && $m[$r][$c-1] =~ /[x-]/ ) { $c--; }

		if ( $d eq "down" && $m[$r+1][$c] =~ /[xi]/ ) { $r++; }
		if ( $d eq "up" && $m[$r-1][$c] =~ /[xi]/ ) { $r--; }

		# CROSS-OVER

		if ( $d eq "right" ) {
			if ( $m[$r][$c+1] eq "7" ) { $d = "down"; $c++;}
			if ( $m[$r][$c+1] eq "J" ) { $d = "up"; $c++;}
		}    
		if ( $d eq "left" ) {
			if ( $m[$r][$c-1] eq "L" ) { $d = "up"; $c--;}
			if ( $m[$r][$c-1] eq "r" ) { $d = "down";  $c--;}
		}

		if ( $d eq "down" ) {
			if ( $m[$r+1][$c] eq "L" ) { $d = "right"; $r++;}
			if ( $m[$r+1][$c] eq "J" ) { $d = "left"; $r++;}
		}    
		if ( $d eq "up" ) {
			if ( $m[$r-1][$c] eq "7" ) { $d = "left"; $r--;}
			if ( $m[$r-1][$c] eq "r" ) { $d = "right"; $r--;}
		}
	}


	# print structure
	$input_structure = "";
	$input_sequence = "";

	my $a = 0;
	my $b = 0;
	$i = 0;
	foreach $a ( @a ) {
		if ( defined $b[$i] && defined $a ) {
			if ( $a > $b[$i] and $b[$i] != 0 ) { 
				if ( $p[$i] eq "p" ) { $input_structure = $input_structure.")"; } 
				if ( $p[$i] eq "b" ) { $input_structure = $input_structure.")"; } 
				if ( $p[$i] eq "!" ) { $input_structure = $input_structure.")"; } 
				if ( $p[$i] eq "*" ) { $input_structure = $input_structure."]"; } 
			}
			if ( $a < $b[$i] and $b[$i] != 0 ) { 
				if ( $p[$i] eq "p" ) { $input_structure = $input_structure."("; } 
				if ( $p[$i] eq "b" ) { $input_structure = $input_structure."("; } 
				if ( $p[$i] eq "!" ) { $input_structure = $input_structure."("; } 
				if ( $p[$i] eq "*" ) { $input_structure = $input_structure."["; } 
			}
			if ( $b[$i] == 0 ) {
				if ($seq[$i] eq "^" ) {
					$input_structure = $input_structure."^";
				} else {
				 $input_structure = $input_structure."."; 
				}
			} 
		} 
		$i++;
	}
	print("\n");

	# print sequence

	my $pri = 0;
    foreach $seq ( @seq ) {
        if ( $seq =~ /[NGACURYKMSWVHBD]/ ) { $input_sequence=$input_sequence.$seq; }
        if ( $seq eq "X" ) { $input_sequence=$input_sequence."N"; }
        if ( $seq eq "T" ) { $input_sequence=$input_sequence."U"; }
    }
	print("\n");
}

#####################################
###  I/O
#####################################


sub printer {
	print $output_spool "$_[0]";
}
