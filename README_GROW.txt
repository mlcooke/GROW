######################################################################################################################
# Usage and general algorithm text file for GROW
# Dec. 2014
# By: Jess McBeck
######################################################################################################################

# GENERAL USAGE ######################################################################################################
 GROW.pl <filename.in>  <increment angle> <start angle>  <end angle> 
        <filename.in> fric2d input file, must exist in the current directory
        <increment angle> the number of degrees between each randomly oriented crack
        <start angle> the angle at which to start searching efficient orietations
        <end angle> the angle at which to end searching

# GENERAL ALGORITHM ##################################################################################################
 When multiple faults growing, test each angle sequentially first, 
    then narrow range and calculate Wext/A for element oriented efficient_angle +/- angle_increment/2 
 For each propagation of crack growth
    Find the names of all the faults growing (@fault)
    Find the range of angles to search (@angles)
    For each fault
       Create hash table describing the fault geometry to be tested for efficiency (fault1_end1 -> angle1)
       Create input file with the geometry
       Test if the growing fault intersects itself, another fault, or a boundary
          Correct intersections if fault intersects another fault or boundary, 
              calculate work for intersecting geometry and set flag in input file 
                that the end of the fault that is intersecting is no longer growing
          Do not correct intersections or calculate work if fault intersects itself
       Calculate external work normalized by fault area for each of the geometries
       Find the fault geometry that optimizes external work normalized by fault area (the most efficient file)
           if non-zero boundary conditions are tractions, most efficient Wext/A is maximum
	   if non-zero boundary conditions are displacements, most efficient Wext/A is minimum
       Rename the most efficient input file <inputname>.eff and the last input.eff as input.prev
              so then in the next iteration of the for loop, GROW will propagate the faults in the updated file
    After this loop executed for all faults, do a "tuning search"
      Keep all faults oriented at most efficient orientation found in sequential search
         Except change the orientation of one fault to (that fault's efficient orientation +/- angle_increment/2)
         Calculate normalized work for each of the "tuning" geometries 
      If the geometry where the orientation of a fault is (that fault's efficient orientation +/- angle_increment/2)
         produces a greater change in work than the efficient orientation found in the first,
         lower resolution search, then the new more efficient orientation of that fault
	 will be the orientation of that fault tested in the remaining tuning steps, and 
	 incorporated into the geometry of the most efficient file found in this propagation of crack growth

# COULOMB FAILURE CRITERION ##########################################################################################
GROW uses the Coulomb criterion to determine if an fracture should propagate (i.e. if an element should be added)
    1) For all orientations of potential elements that do not intersect its own fault, calculate external work normalized by fault area.
    2) For all orientations for which work was calculated, read Fric2D output file to determine if elements
       that were just added slipped or opened according to the tensile criterion and/or Coulomb shear criterion.
       2a) If the element(s) slipped, then allow orientation to be considered for propagation.
       2b) If the element did not slip, then do not allow orientation to be considered.
    3) Of all the remaining orientations (with slipping/opening elements) search for minimum/max work. 
       3a) If one or more orientations failed by the tensile or shear Coulomb criterion 
              identify the most efficient scenario with slipping/opening elements, and use this geometry for the next propagation. 
       3b) If none of the elements slipped or opened then set flag in input file that this fault 
              is no longer propagating, continue propagating other faults if multiple fault tips propagating
 

# EXAMPLE OF INPUT FILE FORMATTING
# *Crack INPUT LINE USAGE ############################################################################################
 *Crack properties represent intact rock values.
 *Fault properties represent fault values.
In FRIC2D input file:
Line that specifies the new fault characteristics of the potential element:
    Must contain the tag *Crack at the beginning of the input file line
    This line is placed after the fault header that contains the fault name, 
    but before the list of individual fault segments 
    For example:
*----------------------
*Fault Conditions
*----------------------
*fault	name	grow_tails?	from_end1?	from_end2?
*num	xbeg	ybeg	xend	yend	stiffS	   stiffN	ten-str	init-coh slid-coh stat-fric dy-fric crit-slip-dist
*----	----	----	----	----	------	   ------	---     ---      -----    -------   ------  ---------
fault	left	no	no	yes					
*tag	                                stiffS	   stiffN	ten-str	shear-str   coh	int-fric dy-fric crit-slip-distance	  
*Crack					1.00E+10   1.00E+10	100	0.00008	    0	 0.6	 0.6	 0 		
41	0	10	20.670	  10	1.00E+10   1.00E+10	100.0	0.0	    0.0	 1.35	 1.35	 0

After one element is added to this fault the input file reads:
*----------------------
*Fault Conditions
*----------------------
*fault	name	grow_tails?	from_end1?	from_end2?
*num	xbeg	ybeg	xend	yend			stiffS	   stiffN	ten-str	init-coh slid-coh stat-fric dy-fric crit-slip-dist
*----	----	----	----	----			------	   ------	---     ---      -----    -------   ------  ---------
fault	left	no	no	yes	
*tag	                                                stiffS	   stiffN	ten-str	shear-str   coh	int-fric dy-fric crit-slip-distance	  
*Crack					        	1.00E+10   1.00E+10	100	0.00008	    0	 0.6	 0.6	 0      	
41	0	10	20.670	        10	        1.00E+10   1.00E+10	100.0	0.0	    0.0	 1.35	 1.35	 0
1	20.670	10	20.31351470	10.35648530	1.00E+10   1.00E+10	100	0.00008	    0	 0.6	 0.6	 0

# WHEN PROPAGATING FRACTURES WHERE TOPOGRAPHY AND GRAVITY IS IMPORTANT ######################################################################
GROW.pl <input.in> <angle> <start angle> <end angle> <topo.topo>
     <input.in> = fric2d input filename  
     <angle> = angle increment (resolution of search for orientations)
     <start angle> the angle at which to start searching efficient orietations
     <end angle> the angle at which to end searching
     <topo.topo> = fric2d topography filename
In FRIC2D input file:
Line that specifies the point from which a new flaw will grow:
    must contain *Flaw at the beginning of the line, 
         then name, x coordinate of point, y coordinate of point, length of 1 element, 
         and remaining characteristics of intact rock through which crack propagates
This line should be placed after the Fault Conditions header

For example:  

*----------------------													
*Fault	Conditions												
*----------------------	
*Tag 	fault 		xcoor 	ycoor 		length 	grow_both? stiffS stiffN 	cohes	friction - s friction -d L 	
*Flaw backthrust	0.04645 0.02037		0.002	yes	   1e10	  1e10		0	0.96	     0.72	0.00025

# DESCRIPTION OF FILES GENERATED #########################################################################################################
Input filename   : foo.in
foo.eff          : fric2d input file of most efficient geometry found so far
foo.raw          : summary of work calculated for various orientations, 
                    nearly identical to standard output
foo_cont<num>.in : If any faults propagate for more than 5 iterations, this program 
                  restarts the process with the same input parameters, except the new input
                  file, which contains the contents of the most efficient geometry found in the 
                  last iteration of crack growth. This automatic restart frees memory.  
standard output: fric2d output statements, work calculated for each scenario, 
                  scenario number of most efficient geometry for each propagation
                  if any scenario did not produce slip on the pupative element just created
Testing filenames: i.e. files generated to test all possible orientations
  foo_prop_fault_end_angle.in : prop = propagation number, fault = fault name, 
                              end = end growing from (1 or 2), 
                              angle = clockwise orientation of potential element
  foo_tune_prop_fault_end_angle.in : same values as file described above, but "tune"
                   indicates that the faults are being tuned
  foo.prev_seq     : the input file of the geometry before a sequential growth of the faults
                   this input file is preserved so the linear search of the most efficient orientations
                   can use this input file to generate correct geometries for new geometries to test
