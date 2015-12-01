#! /usr/bin/perl
# This script does the following tasks:
#    1. Reads a FRIC2D outputfile and 
#    2. calculates both the tectonic and gravity contributions
#	to external work	
#    3. Considers the non-linear influence of loading
#
#	Under most conditions where the side boundaries remain upright,
#	the input files should use the bilateral constraint.
#
#	Possible improvement to this script could include special consideration at
#	model corners
#
#   This scripts works for rectangular models where the top surface is at y=0
#   any variations on this will require modifications

($infile) = @ARGV;
   
# set some constants
$grav = 9.81;

#some initial values
$i=1;
$load = 0;
$gravwork = 0;
$tectWork = 0;

#open(IN, $infile);

############### JAM modification
# make sure that the input file exists 
open IN, $infile, or die "$! when calculating external work: $infile \n";

READING: while (<IN>) {
	chop;
	($first, $second, $third, $fourth, $fifth, $sixth, $seventh, $eight, $ninth, $tenth, $eleven, $twelve, $thirt, $fourteen) = split(" ", $_);	

	# Read lines until you get to the header
		while($first ne "BOUNDARY") {
			$a = <IN>;
			chop($a);
			($first, $second, $third, $fourth, $fifth, $sixth, $seventh, $eight, $ninth, $tenth, $eleven, $twelve, $thirt) = split(" ", $a);
	
			if($fourth eq "Density") {
				$rho = $first/1000000;
			}
			if($third eq "Ratio") {
			        $kratio = $first;
			}
	
			if($third eq "LOADING") {
				$load = $load+1; 
				$tect[$load]=0;
#				print "Loading step # $load\n";
				$endload = $seventh;
			}

			#escape loop once you get to the end of file
			if (eof ){
				last READING;
			}
		} #while

		#read in four lines and throw away
		$a = <IN>;
		$a = <IN>;
		$a = <IN>;
		$a = <IN>;

		#read in the data lines
		$a = <IN>;
		chop($a);
		($first, $second, $third, $fourth, $fifth, $sixth, $seventh, $eight, $ninth, $tenth, $eleven, $twelve, $thirt) = split(" ", $a);


		#read in the length and depth that will be needed for later calculations 	
		if($load == 1) {
			# $i counts the boundary elements that have nonzero forces or displacements
			#	as one of their boundary conditions
		while($first =~ /\d/){
			if( $ninth !=0 || $tenth !=0 || $eleven != 0 || $twelve != 0) { 
				$length[$i] = $sixth;
				$angle[$i] = $seventh;
				$elem[$i] = $first;
				$i=$i+1;
			} #if

		$a = <IN>;
		chop($a);
		($first, $second, $third, $fourth, $fifth, $sixth, $seventh, $eight, $ninth, $tenth, $eleven, $twelve, $thirt) = split(" ", $a);
		} #while

		#read in a couple lines before the data begin
		$a = <IN>;
		$a = <IN>;
		$a = <IN>;
		chop($a);
                ($first, $slip, $third, $fourth, $opening, $sixth, $seventh, $eighth, $ninth) = split(" ", $a);
		} #if
	
		#Now read in the results and calculate the tectonic load	
		$j = 1;
 		while($first =~ /\d/){

			# want to start on the second element because the first
			# has some undesirable corner effects.  Will count the 
			# second and second to last elements twice  <-- implement later

			# Setting stress and displacments into arrays so that we can integrate over
			# the loading path
			$normdisp[$load][$j] = $sixth;
			$normstress[$load][$j] = $ninth;
			$sheardisp[$load][$j] = $third;
			$shearstress[$load][$j] = $eighth;
			$oldload  = $load -1;	
			
			# For bilateral constraint the gravity stress on the dsides is different than the top and
			# bottom
			
			 if($first == $elem[$j]) {
				
			 	# Calculate tectonic work on this element
				# Normal force times displacement is calculated and then shear forces times displacement
				# This work is path dependent and must be integrated over laoding to be the area under
				# the laoding curve. I do this by summing the area of the rectangle and the square under
				# each new segment of the loading path. 
				 if($load==1) {
				       	$tectWork +=  0.5 * $length[$j] * $sixth * ($ninth);
					$tectWork +=   0.5 * $third*$length[$j] * ($eighth);
				}else{
					$diffnormdisp = $sixth-$normdisp[$oldload][$j];
					$diffsheardisp = $third-$sheardisp[$oldload][$j];
					$diffnormstress = $normstress[$load][$j]-$normstress[$oldload][$j];
					$diffshearstress = $shearstress[$load][$j]-$shearstress[$oldload][$j];
					
					$normsquare =  $length[$j] * $diffnormdisp * ($normstress[$oldload][$j]);
					$shearsquare =  $length[$j] * $diffsheardisp * $shearstress[$oldload][$j];	
					
					$normtri =  0.5 * $length[$j] * $diffnormdisp * ($diffnormstress);
					$sheartri=  0.5 * $length[$j] * $diffsheardisp * ($diffshearstress);
					
					$tectWork += $normsquare + $shearsquare + $normtri + $sheartri;
				}
			
				                                                                                                                                         
				#print "$first $j $tectWork \n";
			       	$j = $j + 1;
                        }
        	        $a = <IN>;
                	chop($a);
			($first, $slip, $third, $fourth, $opening, $sixth, $seventh, $eighth, $ninth) = split(" ", $a);

                } #while

	# Calculate the tectonic external work 
	#	work at each step.
		$work = $tectWork *1000000;	
	#print"at load $load tectonic work $work Joules\n";
}

	#print "$infile, $work\n";
#print "WORK CALCED IN OTHER SCRIPT: $work\n";

# JAM modification 
# for the wrapper version of this, write the work to an output file
$new_filename = $infile.'-extWork';
open FILE_HANDLER, ">", $new_filename or die $!;
print FILE_HANDLER "$work";
close FILE_HANDLER;
