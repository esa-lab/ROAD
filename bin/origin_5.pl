use strict;
use warnings FATAL => qw ( all );
use Data::Dumper;
use Math::Trig ':radial';
use Math::Trig;

# >>>>>>>>>>>>>>>>>> RUN PROGRAM <<<<<<<<<<<<<<<<<<<<

my ( $file, $line, @cols, $ATOM, @ATOM );

( $file ) = @ARGV;

if ( not open FILE, "< $file" ) {
    die "file not found!";
}

$ATOM = 0;
@ATOM = ( );
while ( $line = <FILE> ) {
    @cols = ( );
    @cols = split(/\s+/, $line);
    if ( defined $cols[0] ) {
        if ( $cols[0] eq "ATOM") { 
            push @ATOM, {
                "a" => substr("$line", 13, 3),  #atom type
                "n" => substr("$line", 17, 3),  #residue type
                "i" => int(substr("$line", 22, 4)),  #residue number
                "h" => substr("$line", 21, 1), #chain name
                "x" => substr("$line", 30, 8),
                "y" => substr("$line", 38, 8),
                "z" => substr("$line", 46, 8),
                "j" => "N",    #set each atom to N, which means it might be an O3' position
            };
        }
        if ( $cols[0] eq "HETATM") { 
            push @ATOM, {
                "a" => substr("$line", 13, 3),  #atom type
                "n" => substr("$line", 17, 3),  #residue type  ##for HETATM the residue type always gets printed 2 positions earlier in swisspdb!
                "i" => int(substr("$line", 22, 4)),  #residue number
                "h" => substr("$line", 21, 1), #chain name
                "x" => substr("$line", 30, 8),
                "y" => substr("$line", 38, 8),
                "z" => substr("$line", 46, 8),
                "j" => "N",    #set each atom to N, which means it might be an O3' position
            };
        }        
    }
}

#print Dumper ( @ATOM );
#exit;
# PRINT PDB FILE

my $i = 0;
my $x = 1;
my $y = 1;
my $z = 1;
my $c = $ATOM[0]->{i}; 
my $W = 1;
my $AT = "";
my $t = 1;
my $residue=1;
my $thresh=0;
my $distance = 100;
my $cur_chain = "";
my $verbose = 1; #toggle verbose mode  1=verbose, 0=quiet

my $c_one_x=0; my $c_one_y=0; my $c_one_z=0;  #coords of all the C1-primes
my $c_four_x=0; my $c_four_y=0; my $c_four_z=0;  #coords of all the C4-primes
my $n_x=0; my $n_y=0; my $n_z=0;  #coords of all the N1 or N9

#------
#Calculate these from the PDB:
my $sum_zz = 0; my $n_sum = 0;  #number of atoms
my $com_x = 0;
my $com_y = 0;
my $com_z = 0;  #these are sums of the X,Y,Z coords that will set the translation vector for centering the PDB
#------

# first find P and C3* of first nucleotide
# store P and C3* xyz
# then look at next nucleotide

print "REMARK ---File Generated by PLaneFitPDB.pl - written by Cody Geary \n";

# Count residues

foreach $ATOM ( @ATOM ) {      #counts assuming the residues are all numbered differently
    if ( $ATOM->{i} != $c ) {
        $residue++;
        $c = $ATOM->{i};
    }
}   
if ($verbose ==1) {print "REMARK "; print $residue; print " Residues in the PDB file \n";}

$c = $ATOM[0]->{i};   #start at the first residue in the file
#$c = 1;  #start at the residue numbered 1

$cur_chain = $ATOM[0]->{h};
if ($verbose ==1) {print "REMARK starting on chain: "; print $cur_chain; print "\n";}

foreach $ATOM ( @ATOM ) {
	if ( $ATOM->{i} == $c ) {   #$c is the current resi, atom(i) is the resi number element in the PDB file.
								#here, for each residue we find some common atoms to use
		if ( $ATOM->{a} eq "O3*" || $ATOM->{a} eq "O3'") {   #look for two different spellings of O3'
			$x = $ATOM->{x};
			$y = $ATOM->{y};
			$z = $ATOM->{z};
		}
		if ( $ATOM->{a} eq "C1*" || $ATOM->{a} eq "C1'") {   #look for two different spellings of C1'
			$c_one_x = $ATOM->{x};
			$c_one_y = $ATOM->{y};
			$c_one_z = $ATOM->{z};
		}
		if ( $ATOM->{a} eq "C4*" || $ATOM->{a} eq "C4'") {   #look for two different spellings of C4'
			$c_four_x = $ATOM->{x};
			$c_four_y = $ATOM->{y};
			$c_four_z = $ATOM->{z};
		}
		if ( $ATOM->{a} eq "C2*" || $ATOM->{a} eq "C2'") {   #look for atoms
			$n_x = $ATOM->{x};
			$n_y = $ATOM->{y};
			$n_z = $ATOM->{z};
		}			
	}
}
    
my @angles = &calc_ref_frame ($c_one_x, $c_one_y, $c_one_z, $c_four_x, $c_four_y, $c_four_z, $n_x, $n_y, $n_z);  #calculate the ref frame of the current nt

my @vector1 = ();   #translates C1' of current nt to the origin
$vector1[0]=-$c_one_x;
$vector1[1]=-$c_one_y;
$vector1[2]=-$c_one_z;
$vector1[3]=1;

		
foreach $ATOM ( @ATOM ) {
	if ( $ATOM->{j} eq "N" ) { 
		printf "%-6s", "ATOM";
		printf "%5s", $t;
		print "  ";
		print $ATOM->{a};
		print " ";
		print $ATOM->{n};
		print " ";
		print $ATOM->{h};  #Retain chain names
		##print "A";           #or, force all chains to be A
		printf "%4s", $ATOM->{i};
		print "    ";
		
					
		my $coord_x = $ATOM->{x};
		my $coord_y = $ATOM->{y};
		my $coord_z = $ATOM->{z};
		
		## Here we translate and rotate our final data set:
		
		my @point = ();
		$point[0]=$coord_x;		$point[1]=$coord_y;		$point[2]=$coord_z;		$point[3]=1;
		
		my @moved_coords = &translate_matrix ("@point","@vector1");  #POINT followed by translation vector
		@moved_coords = &rotate_z ("@moved_coords",$angles[0]);  #first rotation
		@moved_coords = &rotate_y ("@moved_coords",$angles[1]); #second rotation	
		@moved_coords = &rotate_z ("@moved_coords",$angles[2]); #final rotation
	
	
		my $roundedx =0; my $roundedy=0; my $roundedz=0;
	
		$roundedx = sprintf("%4.3f", $moved_coords[0]);
		$roundedy = sprintf("%4.3f", $moved_coords[1]);
		$roundedz = sprintf("%4.3f", $moved_coords[2]);

		printf("%8s", $roundedx);
		printf("%8s", $roundedy);
		printf("%8s", $roundedz);
		
				print "  1.00100.00   ";
		print "\n";
		$t++;

		$ATOM->{j} = "X";  #Mark the position as found
	}
}   

print "END 1\n";


#-------------------------------subroutines


sub calc_ref_frame {

	my @point2 = ();   #C1'
	$point2[0]=$_[0];
	$point2[1]=$_[1];
	$point2[2]=$_[2];
	$point2[3]=1;

	my @point3 = ();   #C4'
	$point3[0]=$_[3];
	$point3[1]=$_[4];
	$point3[2]=$_[5];
	$point3[3]=1;
	
	my @point4 = ();   #C3
	$point4[0]=$_[6];
	$point4[1]=$_[7];
	$point4[2]=$_[8];
	$point4[3]=1;

	#Calculate the translation vector we need
	
	my @vector1 = ();   #translates C1' to the origin
	$vector1[0]=-$_[0];
	$vector1[1]=-$_[1];
	$vector1[2]=-$_[2];
	$vector1[3]=1;
	
	my @c4vect = &translate_matrix ("@point3","@vector1");	#look at C4' relative to C1'	
	my($rho, $theta, $phi)  = cartesian_to_spherical($c4vect[0], $c4vect[1], $c4vect[2]); #calculate the rotation angles we need to put C4 on the x-axis

	my @nvect = &translate_matrix ("@point4","@vector1"); #look at N(now c3') relative to C1, so we can follow and calculate the y-axis rotation	
	@nvect = &rotate_z ("@nvect",$theta);  #first rotation
	@nvect = &rotate_y ("@nvect",$phi); #second rotation
	if($nvect[0]==0){$nvect[0] = 0.000000000000001;  print("REMARK - ERROR DIVIDE BY ZERO \n");} #avoid divide by zero - 
	my $last_angle = atan($nvect[1]/$nvect[0]);
	@nvect = &rotate_z ("@nvect",$last_angle); #check the position of the c3' after the last move
	if ($nvect[0]<0){$last_angle += 3.141592654;}  #there are 2 possible orientations in the plane, so we flip all the negative ones

	($theta, $phi, $last_angle);  #these angles correspond to rotate_z, rotate_y, and rotate_z in that order
		
}


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
	my $vout=$m1[0]*$m2[3] + $m1[1]*$m2[7] + $m1[2]*$m2[11] + $m1[3]*$m2[15];  #we don't use the last digit, so no need to calculate it
	($xout, $yout, $zout, $vout);
}