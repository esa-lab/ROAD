#!/usr/bin/env perl   

#  -*- perl -*-

# RNA-flipper
# written by Ebbe S. Andersen <esa@inano.au.dk>, 2013
#
# Updated by Cody Geary, April 2015
#  Added support for NGACURYKMSWVHBD
#  Added support for -,|,+ in plain ascii, which then converts to special characters
#  Removed space-padding requirements

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

my ( $file, @lines, $line, @cols, $seq );

( $file ) = @ARGV;

	@lines = &read_file( $file ); #use the read_file sub to read in the input file


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
my $wide =0;
my $widest = 0;
my @cross = ();
my $KL_pattern = '';
	my $find = "";
	my $replace = "";
my $name = 'Untitled';


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
			if ( $cols eq "-" ) { $m[$j][$k] = "-"; $k++; $wide=$k; }     
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



# print back translation

print "Input file:\n";
for ($i=0;$i<1000;$i++) {
    for ($j=0;$j<1000;$j++) {
        if ( defined $m[$i][$j] ) { 
            if ( $m[$i][$j] =~ /[NXGACURYKMSWVHBDT35\*]/ ) { print "$m[$i][$j]"; }
            if ( $m[$i][$j] eq "\n" ){ print "$m[$i][$j]"; }   
            if ( $m[$i][$j] eq " " ) { print "$m[$i][$j]"; }   
            if ( $m[$i][$j] eq "7" ) { print "\342\225\256"; } 
            if ( $m[$i][$j] eq "-" ) { print "\342\224\200"; } 
            if ( $m[$i][$j] eq "r" ) { print "\342\225\255"; } 
            if ( $m[$i][$j] eq "L" ) { print "\342\225\260"; } 
            if ( $m[$i][$j] eq "J" ) { print "\342\225\257"; } 
            if ( $m[$i][$j] eq "i" ) { print "\342\224\202"; } 
            if ( $m[$i][$j] eq "p" ) { print "\342\224\212"; } 
            if ( $m[$i][$j] eq "b" ) { print "="; }
            if ( $m[$i][$j] eq "!" ) { print "!"; }   
            if ( $m[$i][$j] eq "^" ) { print "^"; }            
            if ( $m[$i][$j] eq "x" ) { print "\342\224\274"; } 
        }
    }
}
print "\n\n";

print "Flipped file:\n";
for ($i=1000;$i>0;$i+=-1) {
    for ($j=0;$j<1000;$j++) {
        if ( defined $m[$i][$j] ) { 
            if ( $m[$i][$j] =~ /[NXGACURYKMSWVHBDT35\*]/ ) { print "$m[$i][$j]"; }
            if ( $m[$i][$j] eq "\n" ){ print "$m[$i][$j]"; }   
            if ( $m[$i][$j] eq " " ) { print "$m[$i][$j]"; }   
            if ( $m[$i][$j] eq "J" ) { print "\342\225\256"; } 
            if ( $m[$i][$j] eq "-" ) { print "\342\224\200"; } 
            if ( $m[$i][$j] eq "L" ) { print "\342\225\255"; } 
            if ( $m[$i][$j] eq "r" ) { print "\342\225\260"; } 
            if ( $m[$i][$j] eq "7" ) { print "\342\225\257"; } 
            if ( $m[$i][$j] eq "i" ) { print "\342\224\202"; } 
            if ( $m[$i][$j] eq "p" ) { print "\342\224\212"; } 
            if ( $m[$i][$j] eq "b" ) { print "="; }
            if ( $m[$i][$j] eq "!" ) { print "!"; }  
            if ( $m[$i][$j] eq "^" ) { print "^"; }              
            if ( $m[$i][$j] eq "x" ) { print "\342\224\274"; } 
        }
    }
}

print "Reversed file:\n";
for ($i=0;$i<1000;$i++) {
    for ($j=1000;$j>-1;$j+=-1) {
        if ( defined $m[$i][$j] ) { 
            if ( $m[$i][$j] =~ /[NXGACURYKMSWVHBDT35\*]/ ) { print "$m[$i][$j]"; }
            if ( $m[$i][$j] eq "\n" ){ print "$m[$i][$j]"; }   
            if ( $m[$i][$j] eq " " ) { print "$m[$i][$j]"; }   
            if ( $m[$i][$j] eq "r" ) { print "\342\225\256"; } 
            if ( $m[$i][$j] eq "-" ) { print "\342\224\200"; } 
            if ( $m[$i][$j] eq "7" ) { print "\342\225\255"; } 
            if ( $m[$i][$j] eq "J" ) { print "\342\225\260"; } 
            if ( $m[$i][$j] eq "L" ) { print "\342\225\257"; } 
            if ( $m[$i][$j] eq "i" ) { print "\342\224\202"; } 
            if ( $m[$i][$j] eq "p" ) { print "\342\224\212"; } 
            if ( $m[$i][$j] eq "b" ) { print "="; }
            if ( $m[$i][$j] eq "!" ) { print "!"; } 
            if ( $m[$i][$j] eq "^" ) { print "^"; }               
            if ( $m[$i][$j] eq "x" ) { print "\342\224\274"; } 
        }
    }
}
print "\n\n";

print "Reversed and Flipped file:\n";
for ($i=1000;$i>-1;$i+=-1) {
    for ($j=1000;$j>-1;$j+=-1) {
        if ( defined $m[$i][$j] ) { 
            if ( $m[$i][$j] =~ /[NXGACURYKMSWVHBDT35\*]/ ) { print "$m[$i][$j]"; }
            if ( $m[$i][$j] eq "\n" ){ print "$m[$i][$j]"; }   
            if ( $m[$i][$j] eq " " ) { print "$m[$i][$j]"; }   
            if ( $m[$i][$j] eq "L" ) { print "\342\225\256"; } 
            if ( $m[$i][$j] eq "-" ) { print "\342\224\200"; } 
            if ( $m[$i][$j] eq "J" ) { print "\342\225\255"; } 
            if ( $m[$i][$j] eq "7" ) { print "\342\225\260"; } 
            if ( $m[$i][$j] eq "r" ) { print "\342\225\257"; } 
            if ( $m[$i][$j] eq "i" ) { print "\342\224\202"; } 
            if ( $m[$i][$j] eq "p" ) { print "\342\224\212"; }
            if ( $m[$i][$j] eq "b" ) { print "="; } 
            if ( $m[$i][$j] eq "!" ) { print "!"; }   
            if ( $m[$i][$j] eq "^" ) { print "^"; }             
            if ( $m[$i][$j] eq "x" ) { print "\342\224\274"; } 
        }
    }
}
