#!/usr/bin/env perl   

#  -*- perl -*-

# RNA-trace
# written by Ebbe S. Andersen <esa@inano.au.dk>, 2013
#
# Updated by Cody Geary, April 2015
#  Added support for NGACURYKMSWVHBDT
#  Added support for -,|,+ in plain ascii, which then converts to special characters
#  Removed space-padding requirement at the end of most lines (only the first line must be padded enough)
#
#  Converts T to U
#
#  June20, 2015 edit: bug fixes in sequence read-in parser
#
#  Fixed support for "!" extra base pairs.
#  Added strand-path diagram as extra output at the end
#
#  June30 2017 - Update to generate target.txt outputs  use: perl trace_pattern.pl pattern.txt > target.txt

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
use Data::Dumper;
use Encode;
use utf8;


# >>>>>>>>>>>>>>>>>> RUN PROGRAM <<<<<<<<<<<<<<<<<<<<

my ( $file1, $file2, $line, @lines, @cols, $seq );

( $file1, $file2 ) = @ARGV;


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



my @pri = ( );
if ( defined $file2 ) {
	@lines = &read_file( $file2 );
    foreach $line( @lines ) {
        print "$line";
        print "\n";
        @pri = split(//, $line);
    }
} else {
    $file2 = "no";
}

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

my @m = ( );
my @n = ( );
my @t = ( );
@cols = ( );
my $i = 0;
my $j = 0;
my $cols = 0;
my $l = 0;
my $k = 0;
my $maxlength = 0;
my $endpadding = 0;
my $numtees = 0;
my $name = 'Untitled';
my $KL_pattern = '';
my $wide =0;
my $widest = 0;
my @cross = ();
	my $find = "";
	my $replace = "";


	foreach $line ( @lines ) {
		$l = length $line;
		
		$k = 0;
		$wide=0;
		
		if ($j<1){$maxlength = $l;}


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
			if ( $cols eq "^" ) { $m[$j][$k] = "i"; $cross[$j][$k] = "^"; $k++;$wide=$k;}  
			   
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


# print back translation

#print "Input file:\n";

for ($i=0;$i<1000;$i++) {
    for ($j=0;$j<1000;$j++) {
        if ( defined $m[$i][$j] ) { 
            if ( $m[$i][$j] eq "T" ) { $numtees++; }
        }
    }
}


# Find 5 prime end

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
                    if ( $m[$r+1][$c] =~ /[NXGACURYKMSWVHBDTi]/ ) { $d = "down"; }
                    if ( $m[$r-1][$c] =~ /[NXGACURYKMSWVHBDTi]/ ) { $d = "up"; }
                    #print "\nThe 5p end is found at row $r, column $c and is running $d.\n"; 
                }
            }
        }
    }
}

# Find number of nucleotides in blueprint

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


#Trace the structure

my $r2 = scalar ($r);
my $c2 = scalar ($c);
my $d2 = $d;
my $num = 0;
my @seq = ( );
my $test = "test";
for ($k=0;$k<$nt+10000;$k++) { 
    # TRACE HORIZONTAL STRAND
    if ( $file2 eq "no" ) {
        if ( $d eq "right" && $m[$r][$c+1] =~ /[xNXGACURYKMSWVHBDT-]/ ) { 
            if ( $m[$r][$c+1] =~ /[NXGACURYKMSWVHBDT]/ ) { $num++; $n[$r][$c+1] = $num; push @seq, $m[$r][$c+1]; } 
            $c++;
        }
        if ( $d eq "left"  && $m[$r][$c-1] =~ /[xNXGACURYKMSWVHBDT-]/ ) { 
            if ( $m[$r][$c-1] =~ /[NXGACURYKMSWVHBDT]/ ) { $num++; $n[$r][$c-1] = $num; push @seq, $m[$r][$c-1]; }
            $c--;
        }

        if ( $d eq "up" && $m[$r-1][$c] =~ /[xNXGACURYKMSWVHBDTi]/ ) { 
            if ( $m[$r-1][$c] =~ /[NXGACURYKMSWVHBDT]/ ) { $num++; $n[$r-1][$c] = $num; push @seq, $m[$r-1][$c]; } 
            $r--; 
        }
        if ( $d eq "down" && $m[$r+1][$c] =~ /[xNXGACURYKMSWVHBDTi]/ ) { 
            if ( $m[$r+1][$c] =~ /[NXGACURYKMSWVHBDT]/ ) { $num++; $n[$r+1][$c] = $num; push @seq, $m[$r+1][$c]; } 
            $r++; 
        }



    } else { # when a sequence is available as file2 then we add it on the m-grid here
        if ( $d eq "right" && $m[$r][$c+1] =~ /[xNXGACURYKMSWVHBDT-]/ ) { 
            if ( $m[$r][$c+1] =~ /[NXGACURYKMSWVHBDT]/ ) { $m[$r][$c+1] = $pri[$num]; $num++; $n[$r][$c+1] = $num; }
            $c++;
        }
        if ( $d eq "left"  && $m[$r][$c-1] =~ /[xNXGACURYKMSWVHBDT-]/ ) { 
            if ( $m[$r][$c-1] =~ /[NXGACURYKMSWVHBDT]/ ) { $m[$r][$c-1] = $pri[$num]; $num++; $n[$r][$c-1] = $num; }
            $c--; 
        }


        if ( $d eq "up" && $m[$r-1][$c] =~ /[xNXGACURYKMSWVHBDTi]/ ) { 
            if ( $m[$r-1][$c] =~ /[NXGACURYKMSWVHBDT]/ ) { $m[$r-1][$c] = $pri[$num]; $num++; $n[$r-1][$c] = $num; }
            $r--;
        }
        if ( $d eq "down" && $m[$r+1][$c] =~ /[xNXGACURYKMSWVHBDTi]/ ) { 
            if ( $m[$r+1][$c] =~ /[NXGACURYKMSWVHBDT]/ ) { $m[$r+1][$c] = $pri[$num]; $num++; $n[$r+1][$c] = $num; }
            $r++; 
        }


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


    if ( $m[$r][$c+1] =~ /\d/ ) { if ( $m[$r][$c+1] == 3 ) { $test = "success"; last; } }
    if ( $m[$r][$c-1] =~ /\d/ ) { if ( $m[$r][$c-1] == 3 ) { $test = "success"; last; } }

    if ( $m[$r+1][$c] =~ /\d/ ) { if ( $m[$r+1][$c] == 3 ) { $test = "success"; last; } }
    if ( $m[$r-1][$c] =~ /\d/ ) { if ( $m[$r-1][$c] == 3 ) { $test = "success"; last; } }


}
if ( $test eq "success" ) {
    #print "The structure has been successfully traced (3p end found).\n";
} else {
    print "The trace through the structure failed (3p end not found). Ended at row $r, column $c.\n";
}
#print "The are " . scalar (@seq) . " nts from 5p to 3p.\n";

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
            if ( $m[$r+1][$c+1] =~ /[!p\*]/ ) {  #check above for pairs
                push @a, $n[$r][$c+1];
                push @b, $n[$r+2][$c+1]; 
                push @p, $m[$r+1][$c+1]; 
            } elsif ( $m[$r-1][$c+1] =~ /[!p\*]/ ) {  #check below for pairs
                push @a, $n[$r][$c+1];
                push @b, $n[$r-2][$c+1]; 
                push @p, $m[$r-1][$c+1]; 
            } else {
            	if ( $m[$r][$c+1] =~ /[X]/) {
            		push @a, $n[$r][$c+1];
            		push @b, 0;
            		push @p, $KL_pat[$current_KL];
            		if ($m[$r][$c+2] !~ /[X]/) {$current_KL++;}
            	} else {
            	    push @a, $n[$r][$c+1];
               		push @b, 0; 
               		push @p, "-";
 				}               
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
            	if ( $m[$r][$c-1] =~ /[X]/) {
            		push @a, $n[$r][$c-1];
            		push @b, 0;
            		push @p, $KL_pat[$current_KL];
            		if ($m[$r][$c-2] !~ /[X]/) {$current_KL++;}
            	} else {
            	    push @a, $n[$r][$c-1];
               		push @b, 0; 
               		push @p, "-";
 				}               
            }
        }
        $c--; 
    }

    if ( $d eq "down" && $m[$r+1][$c] =~ /[xNXGACURYKMSWVHBDTi]/ ) { 
        if ( $m[$r+1][$c] =~ /[NXGACURYKMSWVHBDT]/ ) { 
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

    if ( $d eq "up" && $m[$r-1][$c] =~ /[xNXGACURYKMSWVHBDTi]/ ) { 
        if ( $m[$r-1][$c] =~ /[NXGACURYKMSWVHBDT]/ ) { 
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


# print filename
print "$name\n";


# print structure

my $a = 0;
my $b = 0;
$i = 0;
foreach $a ( @a ) {
    if ( defined $b[$i] && defined $a ) {
        if ( $a > $b[$i] and $b[$i] != 0 ) { 
            if ( $p[$i] eq "p" ) { print ")"; } 
            if ( $p[$i] eq "b" ) { print ")"; } 
            if ( $p[$i] eq "!" ) { print "}"; } 
            if ( $p[$i] eq "*" ) { print "]"; }
        }
        if ( $a < $b[$i] and $b[$i] != 0 ) { 
            if ( $p[$i] eq "p" ) { print "("; } 
            if ( $p[$i] eq "b" ) { print "("; } 
            if ( $p[$i] eq "!" ) { print "{"; } 
            if ( $p[$i] eq "*" ) { print "["; } 
        }
        if ( $b[$i] == 0 ) {
        	if ( $p[$i] =~ /[ABCDEFGHIJ1234567890]/ ) { print $p[$i]; } else { print "."; }
        }
    } else {
        if ( $p[$i] =~ /[ABCDEFGHIJ1234567890]/ ) { print $p[$i]; } else { print "."; }
    } 
    $i++;
}
print "\n";


# print sequence

my $pri = 0;
if ( $file2 eq "no" ) {
    foreach $seq ( @seq ) {
        if ( $seq =~ /[NGACURYKMSWVHBD]/ ) { print "$seq"; }
        if ( $seq eq "T" ) { print "U"; }
        if ( $seq eq "X" ) { print "N"; }
    }
    print "\n";
} else { 
    foreach $pri ( @pri ) {
        if ( $pri =~ /[NGACURYKMSWVHBD]/ ) { print "$pri"; }
        if ( $pri eq "T" ) { print "U"; }
        if ( $seq eq "X" ) { print "N"; }
    }
    print "\n";
}


