

#!/usr/bin/env perl   
#  -*- perl -*-
#
#	RNApath.pl
#
#	Written by Cody Geary,  2019
#
#	Requires RNA_nts.pdb datafile
#
#	Usage:  perl RNApath.pl input.file
#	
#	Input file is a traceable text map with >Name as the first line
#	Outputs Name.pdb PDB file keyframes and Chimera .cmd script for making morph movies.
#			

#Copyright 2019  Cody Geary and Ebbe S. Andersen

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

use Cwd;  #to get the current working directory
my $dir = getcwd;


####################
###Begin Program
####################
my $file1 = "";	
my $synthesis_step_size;
my $KL_delay;
my $output_file_path;

($file1, $synthesis_step_size, $KL_delay, $output_file_path ) = @ARGV;

########################
## User Specified Variables
########################

if ( !defined $synthesis_step_size ) { $synthesis_step_size=15; }  # Default number of nts between keyframes
if ( !defined $KL_delay ) { $KL_delay=150; }  # Default number of nts delay before KLs snap closed


my $eccentricity=.1;  #Determines how much single-stranded regions wobble around (0.5 looks like very random motion)
my $frozen = 0; #this var finds the first base-pair, and freezes the nascent chain before this position so it doesn't wobble too much.

my $nascent_chain = 20;  
my $stacking_delay = 60;

my $hidden_offset=.01;  #this sets the offset for both x,y,z for the 'hidden' superimposed residues.
						# an offset of 0 leads to divide by zero errors in some ribbon-models in PDB viewers, so it's important
						# to have his set to a small, but non-zero amount, so the ribbon will have a clear directionality.

my $pmode = 1;   #pmode =0 is all-atom, pmode=1 is P-trace only.  pmode=1 saves a lot of rendering time and storage space.

my $max_design_size = 0;  #0 sets unlimited size, otherwise sets the max accepted length in nts (for the webserver)

#################################################################################

my $synthesized = 0;
my $timesteps= 0;

my @nt_step = ();	my @angles = (); 	my @point = ();	my @moved = ();
my @ref_frame = ();	my @ref_frame_angle = ();	my @ref_trans = ();
my @motif = (); #is each nt_pos a motif, 1 or 0
my $roundedx =0; my $roundedy=0; my $roundedz=0;
my $ref_frame_position = 1;

my $cross_rand_a = rand(1)*$eccentricity-$eccentricity/2;  #randomize variables related to the relaxed crossover (use same random number for the whole construct each pdb iteration)
my $cross_rand_b = rand(1)*$eccentricity-$eccentricity/2;
my $cross_rand_c = rand(1)*$eccentricity-$eccentricity/2;

#################################
#Trace the Strandpath Input
#################################

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


my $input_structure = "";
my $barriers_input_structure = "";
my $input_sequence = "";
my $name = 'mytestoutput';

my @m = ( );	my @n = ( );
my @strand_dir= ( );
&trace;  #trace the strandpath


my @input_seq = split(//,$input_sequence);


my $find = " ";
my $replace = "";
$name =~ s/$find/$replace/g;
#$name=$name.".pdb";

	my ( $file, $line, @lines, @cols, $seq, @seq, $ATOM, @ATOM );
	my @map = (); #stores the dot_paren map
	( $file ) = @ARGV;

	my @pat = split(//, $input_structure);

	if ( not open FILE, "< RNA_lib.pdb" ) { die "file not found!";}	#some error checking stuff
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
	my $i = 0;	my $p = 0; my $x = 1;	my $y = 1;	my $z = 1;
	my $c = $ATOM[0]->{i}; 

	my $strand_length = scalar(@seq);
	my $pattern_length = scalar(@pat);
	my $nt_pos = 1; 
	
	
	####Variables for mapping multiloops
	my $multiloops = 0;
	my @multi = ();
	my @multi_mask = ();
	my @multi_type = ();

#########
# Check strand-length for the server
if ($max_design_size > 0){  #0 sets unlimited design size
	if($strand_length>$max_design_size){
		die "Strand Length is greater than $max_design_size, try a smaller design.";
	}
}

my $keyframes = int(($strand_length+$KL_delay)/$synthesis_step_size)+1;
my $framecount=0;

#######
# Make the Working DIR and write a Command file for Chimera
######

if ( !defined $output_file_path ){
	qx(mkdir data\_$name);
	$output_file_path = "data\_$name";
} else {
	$output_file_path = $output_file_path."\/";
}

open(my $output_spool, '>', "$output_file_path\/$name\_commands.cmd"); ##Output buffer for the screen/file outputs.  Spool will hold all screen outputs and send them to a file at the end.

open(my $output_spool2, '>', "$output_file_path\/$name\_Substructures.txt"); ##Output buffer for the screen/file outputs.  Spool will hold all screen outputs and send them to a file at the end.

########
# Show the Structural Barriers
#########

print $output_spool2 "\n\nHighlighting Structural Barriers\n\n";
print $output_spool2 "  Plot of delay = $KL_delay nts before closing KL interactions. \n\n";

			
	#now parse out the ^ symbols
	my $dot_paren = $input_structure;
	$find = "^";
	$replace = "";
	$find = quotemeta $find; # escape regex metachars if present
	$dot_paren =~ s/$find/$replace/g;

	&map($dot_paren);  #generate the nt pairing map with no ^ symbols first, to initialize the random sequences right
	&map($input_structure);

my @struc_array = split(//, $input_structure);
my $length=scalar @struc_array;

my @barriers = ();

my $topo_count =0;
my $k=0; my $j=0;

my $orange_list = "";
my $red_list = "";

my $crossover_count=0;
for ($i=0; $i<$length; $i++){
	$barriers[$i]="\342\227\246";  #blank
	if ($struc_array[$i]eq"^"){$barriers[$i]="^";}
	if ($struc_array[$i]eq"."){$barriers[$i]="\342\227\246"; $topo_count=0;} #unpaired are not marked
	if ($struc_array[$i]eq"("){$barriers[$i]="x"; $topo_count=0;} #mark open pairs as potential barriers
	if ($struc_array[$i]eq")"){	
		if ($barriers[$map[$i]] eq "x"){  #closing pairs remove barriers
			$barriers[$i]="\342\227\246"; 
			$barriers[$map[$i]]="\342\227\246"; 
			$topo_count=0;
		} 
		if ($barriers[$map[$i]] eq "X"){ #if the closing pair is already blocked, then we mark the barrier reached and start counting
			if($topo_count>5){
				$barriers[$i]="X"; 
				$barriers[$map[$i]]="~";
				
			} else { $barriers[$i]="~"; $barriers[$map[$i]]="~";}
			$topo_count++;  #increment the topology counter
		}  
	}
	if ($i>$KL_delay){
		if ($struc_array[$i-$KL_delay]eq"]"){   #closing KLs marks all the intervening sequence is blocked
			$barriers[$i-$KL_delay]="\342\227\246";  
			$barriers[$map[$i-$KL_delay]]="\342\227\246";
			for($k=$map[$i-$KL_delay]; $k<$i-$KL_delay; $k++){
				if ($barriers[$k]eq"x"){$barriers[$k]="X"; } #mark all the blocked residues
			}

		}
	}
}

my $offsetnum = 0;
for ($i=0; $i<$length; $i++){   #scan through and find the positions to color orange and red
	if ($barriers[$i] eq "^"){$crossover_count++}
	$offsetnum = $i-$crossover_count;
	if ($barriers[$i] eq "~"){ $orange_list=$orange_list."$offsetnum,"}
	if ($barriers[$i] eq "X"){ $red_list=$red_list."$offsetnum,"}
}
chop($orange_list); #chop off the tailing comma
chop($red_list);

sub num2range {   #regex code to convert series of numbers to ranges of numbers
  local $_ = join ',' => @_;
  s/(?<!\d)(\d+)(?:,((??{$++1}))(?!\d))+/$1-$+/g;
  return $_;
}

$orange_list = num2range($orange_list);
$red_list = num2range($red_list);


    for ($i=0;$i<1000;$i++) {
        for ($j=0;$j<1000;$j++) {
            if ( defined $m[$i][$j] || defined $n[$i][$j] ) {  
                if ( $m[$i][$j] =~ /[NXACGURYKMSTWVHBD]/ ) { print $output_spool2 "$barriers[$n[$i][$j]-1]"; }                 
                if ( $m[$i][$j] =~ /[35]/ ) { print $output_spool2 "$m[$i][$j]"; } 
                if ( $m[$i][$j] eq "T" ) { print $output_spool2 "\342\224\200"; } 
                if ( $m[$i][$j] eq "\n" ) { print $output_spool2 "$m[$i][$j]"; }
                if ( $m[$i][$j] eq " " ) { print $output_spool2 "$m[$i][$j]"; }
                if ( $m[$i][$j] eq "^" ) { print $output_spool2 "^"; } 
                if ( $m[$i][$j] eq "7" ) { print $output_spool2 "\342\225\256"; }
                if ( $m[$i][$j] eq "-" ) { print $output_spool2 "\342\224\200"; }
                if ( $m[$i][$j] eq "r" ) { print $output_spool2 "\342\225\255"; }
                if ( $m[$i][$j] eq "L" ) { print $output_spool2 "\342\225\260"; }
                if ( $m[$i][$j] eq "J" ) { print $output_spool2 "\342\225\257"; }
                if ( $m[$i][$j] eq "x" ) { print $output_spool2 "\342\224\274"; }
				if ( $m[$i][$j] eq "i" ) { print $output_spool2 "\342\224\202"; }            
                if ( $m[$i][$j] eq "b" ) { print $output_spool2 "\342\224\200"; }
                
                if (defined $strand_dir[$i][$j] ) {
                	if ( $strand_dir[$i][$j] eq "down" || $strand_dir[$i][$j] eq "up" ){               	
						if ( $m[$i][$j] eq "p" ) { print $output_spool2 "\342\224\200"; }
						if ( $m[$i][$j] eq "!" ) { print $output_spool2 "\342\224\200"; } 
						if ( $m[$i][$j] eq "\*") { print $output_spool2 "\342\224\200"; }   
					} else {
						if ( $m[$i][$j] eq "p" ) { print $output_spool2 "\342\224\202"; }
						if ( $m[$i][$j] eq "!" ) { print $output_spool2 "\342\224\202"; } 
						if ( $m[$i][$j] eq "\*") { print $output_spool2 "\342\224\202"; }   					
					}        
				}
		 
            }
        }
    }



#####################
#  Make the Chimera Script
#####################

	my $numbercount =0;
	my $models_created = 0;

	for ($numbercount=$synthesis_step_size; $numbercount<($keyframes*$synthesis_step_size+1); $numbercount += $synthesis_step_size ){
		&printer ("open $name\_$numbercount.pdb\n");
		$models_created++;
	}
	
	if ($pmode==1){		###in pmode=1, the phosphate trace gives a bug and doesn't render unless we use the final-model as the first frame.
		$models_created+= -1;
		&printer ("morph start \#$models_created\n");  ##sloppy way to show $models-created -1 in the text (decrement, print, increment)
		&printer ("morph interpolate \#0\n");
		$models_created++;
	} else {			
		&printer ("morph start \#0\n");
	}
	
	for ($i=1; $i<$models_created; $i++){
		&printer ("morph interpolate \#$i\n");
	}
	&printer ("~modeldisp #0-$models_created\n");
	&printer ("morph movie\n"); 
	
	&printer ("color blue #$models_created\n");
	if($orange_list ne "" ){&printer ("color orange #$models_created:$orange_list \n");}
	if($red_list ne "" ){&printer ("color red #$models_created:$red_list \n");}
	
	
	#&printer ("roll y 1\n");
	#&printer ("perframe \"center \#$models_created:$strand_length\"");
	
##### Optionally, you can add (in Chimera) this per-frame script to add a polymerase model and center on the 3' end.
###
### close <morph model# +1>
### center <morph model>:<strand_length>
### cofr <model>:<strand_length>
### shape sphere radius 10 center <model>:<strand_length>
###

close $output_spool;

for ($framecount=0; $framecount<$keyframes; $framecount++){

	$timesteps+=$synthesis_step_size; #keep track of this separately since $synthesized maxes at the strand-length.

	$synthesized+=$synthesis_step_size;
	if ($synthesized>$strand_length){  
		$synthesized=$strand_length;
	}  

	#################################
	#Setup Output Spool
	#################################

	open($output_spool, '>', "$output_file_path\/$name\_$timesteps.pdb"); ##Output buffer for the screen/file outputs.  Spool will hold all screen outputs and send them to a file at the end.

	&printer("REMARK -  File Generated by RNA_vis.pl - written by Cody Geary \n");
	&printer("REMARK - Building a total of $strand_length nts. \n");
	##&printer("REMARK - Sequence 5- @seq -3 \n");
	##&printer("REMARK - $input_structure\n");

	for ($i=0; $i<$strand_length; $i++){$motif[$i]=1;} #pre-block all sites, and uncheck them as clean helix is produced.

	my $sub_structure = "";  #this will contain the current substructure
	@struc_array = split(//, $input_structure);
			
	#now parse out the ^ symbols
	$dot_paren = $input_structure;
	$find = "^";
	$replace = "";
	$find = quotemeta $find; # escape regex metachars if present
	$dot_paren =~ s/$find/$replace/g;

	&map($dot_paren);  #generate the nt pairing map with no ^ symbols first, to initialize the random sequences right
	&map($input_structure);  #next re-generate the nt pairing map with ^ symbols added back for making the structural array

	my @dot_paren_array = split(//,$dot_paren);

	my $crossover_symbol_count = scalar(@struc_array)-scalar(@dot_paren_array);
	my $first_pair = 0;
	for ($i=0; $i<scalar( @struc_array ); $i++){
		if ($struc_array[$i] eq "(" && $map[$i]<($synthesized+$crossover_symbol_count) ){$first_pair = 1;}
		if ($i>($synthesized+$crossover_symbol_count)) {
			$sub_structure = $sub_structure."-";
		} elsif ( ( $struc_array[$i]eq"(" || $struc_array[$i]eq"[" ) && $map[$i]>($synthesized+$crossover_symbol_count) ){
			if ($first_pair==0){
				$sub_structure = $sub_structure.",";  #make this "." if you want the start strand to be rigid helix
				$frozen=$i;
			}elsif($first_pair==1){
				$sub_structure = $sub_structure.",";
			}
		} elsif ( $struc_array[$i]eq"^" ) {
			$sub_structure = $sub_structure."^";
		} else {
			$sub_structure = $sub_structure.$struc_array[$i];
		}
	}
	
	&printer2($synthesized,"$sub_structure\n");

	$find = "-";
	$replace = ".";
	$find = quotemeta $find; # escape regex metachars if present
	$sub_structure =~ s/$find/$replace/g;
	
	&map($dot_paren);  #generate the nt pairing map without ^ symbols for RNA building
	@dot_paren_array = split(//,$dot_paren);	
	
	########
	### Find Multiloops
	########
	
	my $mpos = 0;
	$multiloops = 1;
	for ($i=0; $i<$strand_length;$i++){$multi_mask[$i]="-";}
	for ($i=0; $i<$strand_length-1; $i++){
		if ( ($map[$i] ne ($map[$i+1]+1) ) && ( $map[$i] ne $i ) && $dot_paren_array[$i+1]ne"[" && $dot_paren_array[$i+1]ne"]" && $dot_paren_array[$i+1]ne"." ) {		
			if (defined $multi[$i]){
				#&printer(".");
			}else{
				$mpos = $i;
				for (my $j=0; $j<20; $j++){  #this is more awkward than a do-while but safer since it will only count to a max of 20!
					$mpos = $map[$mpos+1];
					$multi[$mpos]="$multiloops";
					$multi_mask[$mpos] = "X";
					#&printer("$mpos ");
					if ($mpos eq ($strand_length-1)){$multi[$i]="0"; $j=21;}
					if ($mpos == $i){$multiloops++; $j=21;}			
				}
				#&printer("- ");
			}
		} else {
			$multi[$i]="-";
			$multi_mask[$mpos]="-";
		}
	}
	
	my @closed_structure = ();	
	for ($i=0; $i<$strand_length; $i++){  #the strand starts pre-blocked as 'O' unfolded
		$closed_structure[$i] = "O";
	}
	
	my $l=0;	
	for ($i=0; $i<$strand_length; $i++){
		if ( $dot_paren_array[$i]eq"]" && ( $timesteps > ($i + $KL_delay)  ) ) {
			for ($p=$map[$i];$p<$i;$p++){
				if ($multi[$p]ne"-"){
					for (my $l=0; $l<($strand_length-1); $l++){
						if ($multi[$p] eq $multi[$l]){$closed_structure[$l]="M";}
					}					
				}
				if ($multi[$map[$p]]ne"-"){
					for (my $l=0; $l<($strand_length-1); $l++){
						if ($multi[$map[$p]] eq $multi[$l]){$closed_structure[$l]="m";}
					}					
				}				
				$closed_structure[$p]="L";
			}
		}
	}
	
#	if ($timesteps>$strand_length){
#		my $end_steps = int( (($timesteps-$strand_length)*$strand_length)/($KL_delay*2) );
#		for ($i=0; $i<$end_steps; $i++){
#			$closed_structure[$i] = "E";
#		}
#		for ($i=$strand_length; $i>($strand_length-$end_steps); $i--){
#			$closed_structure[$i] = "E";
#		}
#	}
	
	for ($i=0; $i<$strand_length; $i++){  #now look at the pairs of anything marked closed, and close that and the neighbors of it too.
		if ( $closed_structure[$i] eq "O" && $closed_structure[$map[$i]] ne "O"){
			$closed_structure[$i] = "D";
		}
	}

####### bug testing outputs
#	
#	&printer("REMARK dotpar- $dot_paren \n");
#	
#	my $output = join "",(@closed_structure);
#	&printer("REMARK closed- $output \n");
#
#	$output = join "",(@multi_mask);
#	&printer("REMARK multi - $output \n\n");
#
########
			
	if($framecount==$keyframes){  #if it's the last time around, all sites are marked 'X' as closed
		#$name=$name."\_Final";
		for(my $z=0; $z<$strand_length; $z++){$closed_structure[$z]="X";}
	}

	#Parse the input structure into motif chunks.  This method is rather clumsy.
	# Initially find all the loops and KLs, these are terminal to hairpins, and thus easiest to identify.  KLs motif placement inlcludes the closing bp.
	# BEWARE: Loops of lengths that are not 4,7 or 9 will not be recognized and modeled as loops!!! We can add other loop lengths to the motif library to correct this.
	my $parsed_struc = $sub_structure;
	
	$find = "((((((.(.((((....)))))))))))";  #Tar/Tat RNA length28
	$replace = "R";
	$find = quotemeta $find; # escape regex metachars if present
	$parsed_struc =~ s/$find/$replace/g;

	$find = "(((((.((((......)))))))))";  #PP7  length25
	$replace = "N";
	$find = quotemeta $find;
	$parsed_struc =~ s/$find/$replace/g;

	$find = "((.((....))))";  #MS2 length13
	$replace = "M";
	$find = quotemeta $find;
	$parsed_struc =~ s/$find/$replace/g;

	$find = "(((..((............(.(";  #iSpinach A Length22
	$replace = "S";
	$find = quotemeta $find;
	$parsed_struc =~ s/$find/$replace/g;

	$find = ")))..))............).)";  #iSpinach A
	$replace = "S";
	$find = quotemeta $find;
	$parsed_struc =~ s/$find/$replace/g;

	$find = "(((..((..............(";  #iSpinach A
	$replace = "S";
	$find = quotemeta $find;
	$parsed_struc =~ s/$find/$replace/g;

	$find = ")))..))..............)";  #iSpinach A
	$replace = "S";
	$find = quotemeta $find;
	$parsed_struc =~ s/$find/$replace/g;

	$find = ").).........)))))";  #iSpinach B Length17
	$replace = "U";
	$find = quotemeta $find;
	$parsed_struc =~ s/$find/$replace/g;

	$find = "(.(.........(((((";  #iSpinach B
	$replace = "U";
	$find = quotemeta $find;
	$parsed_struc =~ s/$find/$replace/g;

	$find = "(((.......................)))";  #Mango  length29
	$replace = "Q";
	$find = quotemeta $find;
	$parsed_struc =~ s/$find/$replace/g;

	#Tetraloops 
	$find = "(....)";
	$replace = "T";
	$find = quotemeta $find; # escape regex metachars if present
	$parsed_struc =~ s/$find/$replace/g;

	#180KL loops
	$find = qr/\(\.\.(\[|.|\]){6}\.\)/;  # (..[[[[[[.) (..]]]]]].) or (.........)
	$replace = "K";
	$parsed_struc =~ s/$find/$replace/g;

	#120KL loops
	$find = "(.......)";
	$replace = "L";
	$find = quotemeta $find; # escape regex metachars if present
	$parsed_struc =~ s/$find/$replace/g;

	$find = "([[[[[[[)";
	$replace = "L";
	$find = quotemeta $find; # escape regex metachars if present
	$parsed_struc =~ s/$find/$replace/g;

	$find = "(]]]]]]])";
	$replace = "L";
	$find = quotemeta $find; # escape regex metachars if present
	$parsed_struc =~ s/$find/$replace/g;
	
	$find = "(...((";   #search for K-turn, note that this actually reads any 3nt or 6nt bulge as a Kturn pattern, but no other motifs use 3nt or 6nt at this stage.
	$replace = "B";
	$find = quotemeta $find;
	$parsed_struc =~ s/$find/$replace/g;

	$find = ")...))";
	$replace = "B";
	$find = quotemeta $find;
	$parsed_struc =~ s/$find/$replace/g;

	$find = "((......(";  # B is the 2nd side of the Kturn
	$replace = "C";
	$find = quotemeta $find;
	$parsed_struc =~ s/$find/$replace/g;

	$find = "))......)";
	$replace = "C";
	$find = quotemeta $find;
	$parsed_struc =~ s/$find/$replace/g;
	
	#Crossover Symbols are specially recognized
	$find = qr/.[\^]../;
	$replace = "E";
	$parsed_struc =~ s/$find/$replace/g;
			
	$find = ",";
	$replace = "A";
	$find = quotemeta $find;
	$parsed_struc =~ s/$find/$replace/g;

	#Everything remaining is modeled as a nt stacked on helix
	$find = qr/\W/;
	$replace = "H";
	$parsed_struc =~ s/$find/$replace/g;

	$nt_pos = 1; 

	###  Code conversion table
	###  T=6nts
	###  K=11nts
	###  L=9nts
	###  E=3nts
	###  H=1nt

	#&printer("REMARK - $parsed_struc \n");  ##Structures larger than 1000nts can crash SwissPDB by overflowing the remark buffer, so I am commenting this line out.   It is for debugging the parser on smaller designs.

	my @assembly_instructions = split(//, $parsed_struc);
	
	#&printer("REMARK - $parsed_struc\n");

	## Now this is clumsy, but we scan over the assembly instruction list and the map and add in 'O' everywhere we have a multi-junction.
	my $index=0; #index is where the build is currently, tracked separately from synthesized
	my $k =0; my $count=0;
	
	$cross_rand_a = rand(1)*$eccentricity-$eccentricity/2;  #randomize variables related to the relaxed crossover (use same random number for the whole construct each pdb iteration)
	$cross_rand_b = rand(1)*$eccentricity-$eccentricity/2;
	$cross_rand_c = rand(1)*$eccentricity-$eccentricity/2;


	for ( $k=0; $k< scalar(@assembly_instructions) ; $k++){
		if ($index< $synthesized ){
			if ($index > $timesteps-$nascent_chain){
				if ($assembly_instructions[$k] eq 'H') { &add_nascent; $index++;}
				if ($assembly_instructions[$k] eq 'A') { &add_nascent; $index++;}  #A is single-strand or unpaired
				if ($assembly_instructions[$k] eq 'T') { &add_tetraloop; $index+=6;}
				if ($assembly_instructions[$k] eq 'K') { &add_KL; $index+=11;}
				if ($assembly_instructions[$k] eq 'L') { &add_KLbend; $index+=9;}
				if ($assembly_instructions[$k] eq 'E') { for($count=0; $count<3; $count++){&add_nascent;} $index+=3;}
				if ($assembly_instructions[$k] eq 'O') { for($count=0; $count<3; $count++){&add_nascent;} $index+=3;}	
				
				if ($assembly_instructions[$k] eq 'R') { &add_tar; $index+=28;} #TAR
				if ($assembly_instructions[$k] eq 'N') { &add_pp7; $index+=25;} #PP7
				if ($assembly_instructions[$k] eq 'M') { &add_ms2; $index+=13;} #MS2
		
				if ($assembly_instructions[$k] eq 'S') { for($count=0; $count<22; $count++){&add_nascent;} $index+=22;} #SpinachA
				if ($assembly_instructions[$k] eq 'U') { for($count=0; $count<17; $count++){&add_nascent;} $index+=17;} #SpinachB

				if ($assembly_instructions[$k] eq 'B') { for($count=0; $count<6; $count++){&add_nascent;} $index+=6;} #KturnA
				if ($assembly_instructions[$k] eq 'C') { for($count=0; $count<9; $count++){&add_nascent;} $index+=9;} #KturnB

				if ($assembly_instructions[$k] eq 'Q') {&add_mango; $index+=29;} #Mango	

			} else {
				if ($closed_structure[$index]eq"O"){
					if ($assembly_instructions[$k] eq 'H') { &add_nts; $index++;}
					if ($assembly_instructions[$k] eq 'A') { &add_nascent; $index++;}
					if ($assembly_instructions[$k] eq 'T') { &add_tetraloop; $index+=6;}
					if ($assembly_instructions[$k] eq 'K') { &add_KL; $index+=11;}
					if ($assembly_instructions[$k] eq 'L') { &add_KLbend; $index+=9;}
					if ($assembly_instructions[$k] eq 'E' ) {
					 	if($index>($timesteps-$stacking_delay) ) {&opened; $index+=3;
					 	} else {
					 		&crossover_relaxed; $index+=3; 
					 	}
					 }
					if ($assembly_instructions[$k] eq 'O') { 
						if($map[$index+3]>$synthesized){
							&add_nascent; &add_nascent; &add_nascent; $index+=3;
						} else {
							if ($index > ($timesteps-$stacking_delay) ){
							 	 &add_nts; &add_nts; &add_nts; $index+=3;
							 } else {
							 	&opened;
							 	$index+=3;
							 }
						}
					}	
					
					if ($assembly_instructions[$k] eq 'R') { &add_tar; $index+=28;}
					if ($assembly_instructions[$k] eq 'M') { &add_ms2; $index+=13;}
					if ($assembly_instructions[$k] eq 'N') { &add_pp7; $index+=25; }
		
					if ($assembly_instructions[$k] eq 'S') { &add_spinacha; $index+=22;}
					if ($assembly_instructions[$k] eq 'U') { &add_spinachb; $index+=17;}

					if ($assembly_instructions[$k] eq 'B') { &add_kturna; $index+=6;}
					if ($assembly_instructions[$k] eq 'C') { &add_kturnb; $index+=9;}

					if ($assembly_instructions[$k] eq 'Q') { &add_mango; $index+=29;}	

				} else {	
					if ($assembly_instructions[$k] eq 'H') { &add_nts; $index++;}
					if ($assembly_instructions[$k] eq 'A') { &add_nts; $index++;}
					if ($assembly_instructions[$k] eq 'T') { &add_tetraloop; $index+=6;}
					if ($assembly_instructions[$k] eq 'K') { &add_KL; $index+=11;}
					if ($assembly_instructions[$k] eq 'L') { &add_KLbend; $index+=9;}
					if ($assembly_instructions[$k] eq 'E') { &crossover; $index+=3;}
					if ($assembly_instructions[$k] eq 'O') { &add_nts; &add_nts; &add_nts; $index+=3;}	
					
					if ($assembly_instructions[$k] eq 'R') { &add_tar; $index+=28;}
					if ($assembly_instructions[$k] eq 'M') { &add_ms2; $index+=13;}
					if ($assembly_instructions[$k] eq 'N') { &add_pp7; $index+=25; }
		
					if ($assembly_instructions[$k] eq 'S') { &add_spinacha; $index+=22;}
					if ($assembly_instructions[$k] eq 'U') { &add_spinachb; $index+=17;}

					if ($assembly_instructions[$k] eq 'B') { &add_kturna; $index+=6;}
					if ($assembly_instructions[$k] eq 'C') { &add_kturnb; $index+=9;}

					if ($assembly_instructions[$k] eq 'Q') { &add_mango; $index+=29;}	
										
				}
			}
				
		} else {
			if ($assembly_instructions[$k] eq 'H') { &add_hidden; $index++;}
			if ($assembly_instructions[$k] eq 'A') { &add_hidden; $index++;}
			if ($assembly_instructions[$k] eq 'T') { for($count=0; $count<6; $count++){&add_hidden;} $index+=6;}
			if ($assembly_instructions[$k] eq 'K') { for($count=0; $count<11; $count++){&add_hidden;} $index+=11;}
			if ($assembly_instructions[$k] eq 'L') { for($count=0; $count<9; $count++){&add_hidden;} $index+=9;}
			if ($assembly_instructions[$k] eq 'E') { &add_hidden; &add_hidden; &add_hidden;  $index+=3;}
			if ($assembly_instructions[$k] eq 'O') { &add_hidden; &add_hidden; &add_hidden;  $index+=3;}	
			
			if ($assembly_instructions[$k] eq 'R') { for($count=0; $count<28; $count++){&add_hidden;} $index+=28;} #TAR
			if ($assembly_instructions[$k] eq 'M') { for($count=0; $count<13; $count++){&add_hidden;} $index+=13;} #MS2
			if ($assembly_instructions[$k] eq 'N') { for($count=0; $count<25; $count++){&add_hidden;} $index+=25;} #PP7
	
			if ($assembly_instructions[$k] eq 'S') { for($count=0; $count<22; $count++){&add_hidden;} $index+=22;} #SpinachA
			if ($assembly_instructions[$k] eq 'U') { for($count=0; $count<17; $count++){&add_hidden;} $index+=17;} #SpinachB

			if ($assembly_instructions[$k] eq 'B') { for($count=0; $count<6; $count++){&add_hidden;} $index+=6;} #KturnA
			if ($assembly_instructions[$k] eq 'C') { for($count=0; $count<9; $count++){&add_hidden;} $index+=9;} #KturnB

			if ($assembly_instructions[$k] eq 'Q') { for($count=0; $count<29; $count++){&add_hidden;} $index+=29;} #Mango	
			
		}
	}

	my $t = 1;
	my $seam = 0;
	my $filecount=1;
	my $generated_seq = "";  #for bugchecking
	
	@input_seq = split(//,$input_sequence);
	
	foreach $PDB ( @PDB ) {  #print out the PDB file
		if ($t<99000){$seam = $PDB->{i};}
		
		if ($PDB->{i}<=$seam){
			if($pmode==0 || $PDB->{a} eq 'P  '){
				printf $output_spool "%-6s", "ATOM";
				printf $output_spool "%5s", $t;
				print $output_spool "  ";
				print $output_spool $PDB->{a};
				print $output_spool "   ";
				#print $output_spool $PDB->{n};
				#$generated_seq=$generated_seq.$PDB->{n};
				
				print $output_spool $input_seq[$t-1];  #print the fixed input sequence, rather than the sequence generated for each frame that could vary based on locked motifs
				
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
			open( $output_spool, '>', "$name\_$synthesized\_$filecount.pdb"); #open a new file buffer
			$t=1;
			if($pmode==0 || $PDB->{a} eq 'P  '){		
				printf $output_spool "%-6s", "ATOM";
				printf $output_spool "%5s", $t;
				print $output_spool "  ";
				print $output_spool $PDB->{a};
				print $output_spool "   ";
				#print $output_spool $PDB->{n};
				
				print $output_spool $input_seq[$t-1];  #print the fixed input sequence, rather than the sequence generated for each frame that could vary based on locked motifs
				
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
	
	#$generated_seq=$generated_seq."\n";  #for bugchecking
	#&printer2("S",$generated_seq);       #output sequence at each step for bugchecking
	
	close $output_spool;
	@PDB = ();  ##very important!  Flush out the PDB buffer before the next file is started.

}

close $output_spool2;  #closes the spool with the substructure outputs

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
	if ($map[$nt_pos-1]<($nt_pos-1) && $motif[$map[$nt_pos-1]]==0 ) {  #if map<nt_pos then we are a closing nt, and we calculate from the partner's ref-frame (later we may average the two options)

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

######

sub add_nascent {
	my $no = 1;

	$ref_frame_position = $nt_pos;
	@ref_frame_angle = &update_ref_frame($ref_frame_position);
	

	my $rand_a = rand(1)*$eccentricity-$eccentricity/2;
	my $rand_b = rand(1)*$eccentricity-$eccentricity/2;
	my $rand_c = rand(1)*$eccentricity-$eccentricity/2;
	
	if ($nt_pos < $frozen) {
		$rand_a = $rand_a/5;
		$rand_b = $rand_b/5;
		$rand_c = $rand_c/5;
	}
						
	foreach $ATOM ( @ATOM ) {			#we scan the datafile for the nt type and paste it into the active site
		if ( ($ATOM->{n} eq $seq[$nt_pos-1]) && ($ATOM->{h} eq 'A')) {   #$c is the current resi, atom(i) is the resi number element in the PDB file.							

		#	$nt_step[0]=4.884777778;	$nt_step[1]=-0.976222222;	$nt_step[2]=-2.029777778;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from a section of single-strand model
		#	$angles[0]=-0.785139878+$rand_a;	$angles[1]=0.313181556+$rand_b;	$angles[2]=0.732568308+$rand_c;			

			$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
			$angles[0]=-1.05152505+$rand_a;	$angles[1]=0.46918242+$rand_b;	$angles[2]=1.37954160+$rand_c;			


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
	
	
	$motif[$nt_pos-2] = 0;
}

#######

sub add_hidden {
	my $no = 1;

		$ref_frame_position = $synthesized;
		@ref_frame_angle = &update_ref_frame($ref_frame_position);
							
		foreach $ATOM ( @ATOM ) {			#we scan the datafile for the nt type and paste it into the active site
			if (defined $seq[$nt_pos-1] ){
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

					$moved[0] += ($hidden_offset*($nt_pos+1-$synthesized));
					$moved[1] += ($hidden_offset*($nt_pos+1-$synthesized));
					$moved[2] += ($hidden_offset*($nt_pos+1-$synthesized));


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
		}
		$nt_pos += 1;	#increment nt_counter
	
	
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
		&printer("REMARK     Loop is a GNRA shown as GAAA\n");
		
		
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
						#"n" => $ATOM->{n},  #residue type
						"n" => $seq[$nt_pos-2+$i], ###ADDED
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
		&printer("REMARK     Loop is a UNCG shown as UUCG\n");
		
		
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
						#"n" => $ATOM->{n},  #residue type
						"n" => $seq[$nt_pos-2+$i], ###ADDED
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

sub crossover_relaxed {	
	$ref_frame_position = $nt_pos-1;
	@ref_frame_angle = &update_ref_frame($ref_frame_position);
	
	$cross_rand_a = rand(1)*$eccentricity/2-$eccentricity/4;
	$cross_rand_b = rand(1)*$eccentricity/2-$eccentricity/4;
	$cross_rand_c = rand(1)*$eccentricity/2-$eccentricity/4;
	
	
	foreach $ATOM ( @ATOM ) {			#Copy the backbone of resi 1 of this motif,  o2' has to come at the end, so we pop it off at the end.
		if ( ($ATOM->{h} eq 'F')  &&  ($ATOM->{i} == 2) )	{	#we start with nt 2,and moved back 1 frame						
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
		if ( ($ATOM->{h} eq 'F')  &&  ($ATOM->{i} == 3) )	{							
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
		if ( ($ATOM->{h} eq 'F')  &&  ($ATOM->{i} == 4) )	{							
			$point[0]=$ATOM->{x};	$point[1]=$ATOM->{y};	$point[2]=$ATOM->{z};	$point[3]=1;	

			$nt_step[0]=5.05026531;	$nt_step[1]=0.63351020;	$nt_step[2]=-2.27143878;	$nt_step[3]=1;	#these were measured by nt_diff.pl, averaged from 50bp of A-form generated in Assemble/Chimera
			$angles[0]=-1.05152505+$cross_rand_a;	$angles[1]=0.46918242+$cross_rand_b;	$angles[2]=1.37954160+$cross_rand_c;	#		

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

#########

sub opened {	
	$ref_frame_position = $nt_pos-1;
	@ref_frame_angle = &update_ref_frame($ref_frame_position);
	
	foreach $ATOM ( @ATOM ) {			#Copy the backbone of resi 1 of this motif,  o2' has to come at the end, so we pop it off at the end.
		if ( ($ATOM->{h} eq 'O')  &&  ($ATOM->{i} == 2) )	{	#we start with nt 2,and moved back 1 frame						
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
		if ( ($ATOM->{h} eq 'O')  &&  ($ATOM->{i} == 3) )	{							
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
		if ( ($ATOM->{i} == 4) && ($ATOM->{h} eq 'O') ) {   #$nt_pos is the current position that we grab the ref frame of here

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
		if ( ($ATOM->{h} eq 'O')  &&  ($ATOM->{i} == 4) )	{							
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



#########

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
					#"n" => $ATOM->{n},  #residue type
					"n" => $seq[$nt_pos-2+$i], ###ADDED
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
					#"n" => $ATOM->{n},  #residue type
					"n" => $seq[$nt_pos-2+$i], ###ADDED
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
					#"n" => $ATOM->{n},  #residue type
					"n" => $seq[$nt_pos-1+$i], ###ADDED
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
		elsif ($_ eq "^"){
			$map[$counter] = $counter;	 ##crossover also map to themselves.
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
			my $pick = int rand(4);			#try to put random sequence in where the sequence is not specified.

			if ($pick == 0){$seq[$counter]= 'A'; $seq[$map[$counter]]='U'}
			if ($pick == 1){$seq[$counter]= 'U'; $seq[$map[$counter]]='A'}
			if ($pick == 2){$seq[$counter]= 'G'; $seq[$map[$counter]]='C'}
			if ($pick == 3){$seq[$counter]= 'C'; $seq[$map[$counter]]='G'}
			#die "Error: To create RNApath animations you must fully specify the sequence in the design with A,U,C,G only. \n";  
		}
		$counter ++;
	}
	
}

#####################################

sub trace{
	my ( $line,@lines, @cols, $seq );
	my @pri = ( );
	
	#if ( not open FILE, "< $file1" ) { die "file not found!"; }
	@lines = &read_file( $file1 );

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

	my @t = ( );	@cols = ( );
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
			if ( $cols eq "\/" &&  substr("$line", $i+1, 1)!~ /[+xNXGACURYKMSWVHBDT35-]/ ) { $m[$j][$k] = "J"; $k++; $wide=$k; }     
			if ( $cols eq "\200" ) { $m[$j][$k] = "-"; $k++; $wide=$k;}
			if ( $cols eq "-" ) { $m[$j][$k] = "-"; $k++; }     
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
						$strand_dir[$r][$c]=$d; #store the strand directionallity in the strand_dir array.
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
		$strand_dir[$r][$c]=$d;  #store the current strand direction into the matrix


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
	
		$strand_dir[$r][$c]=$d;
		if ( $d eq "right" ) {
			if ( $m[$r][$c+1] eq "7" ) { $d = "down"; $c++; $strand_dir[$r][$c+1]=$d; }#store the strand directionallity in the strand_dir array.
			if ( $m[$r][$c+1] eq "J" ) { $d = "up"; $c++; $strand_dir[$r][$c+1]=$d; }
		}    
		if ( $d eq "left" ) {
			if ( $m[$r][$c-1] eq "L" ) { $d = "up"; $c--; $strand_dir[$r][$c-1]=$d; }
			if ( $m[$r][$c-1] eq "r" ) { $d = "down";  $c--; $strand_dir[$r][$c-1]=$d; }
		}

		if ( $d eq "down" ) {
			if ( $m[$r+1][$c] eq "L" ) { $d = "right"; $r++; $strand_dir[$r+1][$c]=$d;}
			if ( $m[$r+1][$c] eq "J" ) { $d = "left"; $r++; $strand_dir[$r+1][$c]=$d; }
		}    
		if ( $d eq "up" ) {
			if ( $m[$r-1][$c] eq "7" ) { $d = "left"; $r--; $strand_dir[$r-1][$c]=$d; }
			if ( $m[$r-1][$c] eq "r" ) { $d = "right"; $r--; $strand_dir[$r-1][$c]=$d;}
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
					$strand_dir[$r+1][$c+1]=$d 
				} elsif ( $m[$r-1][$c+1] =~ /[!p\*]/ ) {
					push @a, $n[$r][$c+1];
					push @b, $n[$r-2][$c+1]; 
					push @p, $m[$r-1][$c+1];
					$strand_dir[$r-1][$c+1]=$d
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
					$strand_dir[$r+1][$c-1]=$d
				} elsif ( $m[$r-1][$c-1] =~ /[!p\*]/ ) {
					push @a, $n[$r][$c-1];
					push @b, $n[$r-2][$c-1]; 
					push @p, $m[$r-1][$c-1];
					$strand_dir[$r-1][$c-1]=$d					 
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
					$strand_dir[$r+1][$c+1]=$d 
				} elsif ( $m[$r+1][$c-1] =~ /[b\*]/ ) {
					push @a, $n[$r+1][$c];
					push @b, $n[$r+1][$c-2]; 
					push @p, $m[$r+1][$c-1]; 
					$strand_dir[$r+1][$c-1]=$d
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
					$strand_dir[$r-1][$c+1]=$d
				} elsif ( $m[$r-1][$c-1] =~ /[b\*]/ ) {
					push @a, $n[$r-1][$c];
					push @b, $n[$r-1][$c-2]; 
					push @p, $m[$r-1][$c-1]; 
					$strand_dir[$r-1][$c-1]=$d 
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
	$barriers_input_structure = "";
	$input_sequence = "";

	my $a = 0;
	my $b = 0;
	$i = 0;
	foreach $a ( @a ) {
		if ( defined $b[$i] && defined $a ) {
			if ( $a > $b[$i] and $b[$i] != 0 ) { 
				if ( $p[$i] eq "p" ) { $input_structure = $input_structure.")"; $barriers_input_structure = $barriers_input_structure.")";} 
				if ( $p[$i] eq "b" ) { $input_structure = $input_structure.")"; $barriers_input_structure = $barriers_input_structure.")";} 
				if ( $p[$i] eq "!" ) { $input_structure = $input_structure.")"; $barriers_input_structure = $barriers_input_structure.")";} 
				if ( $p[$i] eq "*" ) { $input_structure = $input_structure."]"; $barriers_input_structure = $barriers_input_structure."]";} 
			}
			if ( $a < $b[$i] and $b[$i] != 0 ) { 
				if ( $p[$i] eq "p" ) { $input_structure = $input_structure."(";  $barriers_input_structure =  $barriers_input_structure."(";} 
				if ( $p[$i] eq "b" ) { $input_structure = $input_structure."(";  $barriers_input_structure =  $barriers_input_structure."(";} 
				if ( $p[$i] eq "!" ) { $input_structure = $input_structure."(";  $barriers_input_structure =  $barriers_input_structure."(";} 
				if ( $p[$i] eq "*" ) { $input_structure = $input_structure."[";  $barriers_input_structure =  $barriers_input_structure."[";} 
			}
			if ( $b[$i] == 0 ) {
				if ($seq[$i] eq "^" ) {
					$input_structure = $input_structure."^";
				} else {
				 $input_structure = $input_structure."."; $barriers_input_structure = $barriers_input_structure."."; 
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


sub printer {  #Main Output Printer
	print $output_spool "$_[0]";
}
sub printer2 {  #Secondary Printer is for printing the substructures list
	printf $output_spool2 ("%8s",$_[0]);
	print $output_spool2 " - $_[1]";
}
