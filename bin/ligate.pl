#!/usr/bin/env perl   
#  -*- perl -*-

# Ebbe Sloth Andersen, esa@inano.au.dk, Mar 2013
# updated by Cody Geary, April 2015
#    - Starts at the residue numbered 1
#    - requires that there's only 1 residue 1, and that all the nts are numbered differently
#    - searches for P-O3' distances within a sphere of user-set radius (rather than a cube!!)
#    - now reads both O3' and O3*, which are used interchangeably in some PDB files#    - ligates preferentially along the same chain, unless it cannot find a link within a threshold of 4 angstroms
use strict;
use warnings FATAL => qw ( all );
use Data::Dumper;

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
                "n" => substr("$line", 19, 1),  #residue type
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
                "n" => substr("$line", 17, 1),  #residue type  ##for HETATM the residue type always gets printed 2 positions earlier in swisspdb!
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
my $cur_chain = "";my $verbose = 0; #toggle verbose mode  1=verbose, 0=quiet

# first find P and C3* of first nucleotide
# store P and C3* xyz
# then look at next nucleotide

# Count residues

foreach $ATOM ( @ATOM ) {      #counts assuming the residues are all numbered differently
    if ( $ATOM->{i} != $c ) {
        $residue++;
        $c = $ATOM->{i};
    }
}   
if ($verbose ==1) {print "REMARK "; print $residue; print " Residues in the PDB file \n";}


#$c = $ATOM[0]->{i};   #start at the first residue in the file
$c = 1;  #start at the residue numbered 1$cur_chain = $ATOM[0]->{h};
if ($verbose ==1) {print "REMARK starting on chain: "; print $cur_chain; print "\n";}
for ($i=0;$i<$residue;$i++) {
    
    while ($W<$residue+1){

    # print the chosen nucleotide

    foreach $ATOM ( @ATOM ) {
        if ( $ATOM->{i} == $c && $ATOM->{j} eq "N" ) { 
            printf "%-6s", "ATOM";
            printf "%5s", $t;
            print "  ";
            print $ATOM->{a};
            print "   ";
            print $ATOM->{n};
            print " ";
            ##print $ATOM->{h};  #Retain chain names
            print "A";          #or, force all chains to be A
            printf "%4s", $W;
            print "    ";
            print $ATOM->{x};
            print $ATOM->{y};
            print $ATOM->{z};
            print "  1.00100.00   ";
            print "\n";
            $t++;
            if ( $ATOM->{a} eq "O3*" || $ATOM->{a} eq "O3'") {   #look for two different spellings of O3'
                $x = $ATOM->{x};
                $y = $ATOM->{y};
                $z = $ATOM->{z};
            }
            $ATOM->{j} = "X";  #Mark the position as found
            if ($W == $residue) { $residue--;}   #clumsy way to exit the while cycle once the last residue is reached and output printed.
        }
    }
    
    # find the next nucleotide 

    foreach $ATOM ( @ATOM ) {
        if ( $ATOM->{a} eq "P  " && $ATOM->{j} eq "N" ) { # many instances  
            $distance = sqrt(($x - $ATOM->{x})**2 +  ($y - $ATOM->{y})**2 + ($z - $ATOM->{z})**2 );
            if ( $distance < 1.8 + $thresh ) {   #sets the initial phosphate bond distance limit
               if ( $ATOM->{h} eq $cur_chain ) { #check of the potential link is the same chain
                $c = $ATOM->{i};              
                $W++;
                $thresh=0.0;   #reset threshold back to 0
                if ($verbose ==1) {print "REMARK  Distance = "; print $distance; print "\n";}
               } elsif ($thresh > 4) {  #Set the threshold for switching to a new chain
                $c = $ATOM->{i};              
                $W++;
                $thresh=0.0;   #reset threshold back to 0
                $cur_chain = $ATOM->{h};
                if ($verbose ==1) {
                     print "REMARK  Switching Chains to "; print $cur_chain; print "\n";
                     print "REMARK  Distance = "; print $distance; print "\n";
                }
               }

            }
        }
    }
    if ($W < $residue+1) {$thresh += 0.2;} #if it didn't finish searching, then increase the distance sphere by .2A until it finds a residue
    if ($thresh > 30) {$W = ($residue+2); print "REMARK    Warning: Chain end not found. "; } #if it cannot find a connection witin 20angstroms, then end.
    }
}

print "END\n";
