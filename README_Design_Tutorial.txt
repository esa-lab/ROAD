##################################################
Tutorial on RNA Design with RNABuild, RNAVis, and Revolvr
##################################################

The following Perl scripts are included in the software package and licensed under the MIT License. They all work with RNA blueprints (text-formatted RNA structure diagrams).

RNAbuild      	
- generates 3D PDB coordinate files from a blueprint .

Revolvr         
- repeatedly mutates a randomized RNA starting sequence until it folds into a target 2D structure under the ViennaRNA folding model.

Batch_revolvr         
- runs Revolvr multiple times to create different candidate designs for the same blueprint.

RNAvis	
- analyzes and animates one potential folding path from a blueprint.

flip_trace
- outputs a blueprint in all four possible orientations.

trace_analysis 
- parses a blueprint and outputs analyses and helpful variants of the blueprint. 

trace_pattern  
- parses a blueprint and outputs the dot-paren notation for the structure in the blueprint, as well as the sequence constraints. It creates an input for Revolvr.

trace  
- threads a sequence onto a blueprint. Only used by Revolvr.

#########################################
Installation requirements:
To run RNAbuild you will need to first install a Perl environment.
Next, to run Revolvr, you will additionally need to install the ViennaRNA package.
Finally, you will need a plaintext editor that can copy/paste blocks of fixed-width text (using so-called “block” or “column” mode) such as TextEdit on MacOS, or Notepad++ and Vim on Windows.

This software was written and tested with Perl v5.18.2 for MacOS, and ViennaRNA versions 2.1.9 and 2.4.10. Installation for Perl and ViennaRNA should be a few minutes. Software was tested on a Mid 2012 Macbook Pro  running El Capitan 10.11.6 and a Mid 2015 Macbook Pro running Sierra 10.12.6.

#########################################

To Run RNABuild you will need to create an input file in a text editor and save it in the same directory as the RNABuild.pl script file.

The first line of the input file contains a ‘>’ followed by the name of the design:
===============================================================================
>My_Pattern                                             

===============================================================================

Next, the RNA is drawn as a secondary structure diagram (blueprint).  The blueprint is composed of text symbols in a grid of blank spaces.  A text editor capable of displaying fixed-width characters and copy/pasting columns of text is needed.  In MacOS the TextEdit app with the Menlo font works for this purpose, and holding the 'option' key allows rectangles of text characters to be copied and moved like modular blocks.

The 5' end of the strand is indicated by a ‘5’, and a ‘3’ marks the 3' end of the strand.  A valid blueprint will have exactly one ‘5’ and one ‘3’, and should not contain multiple RNA strands.

A dash ‘-’ symbol is used to indicate the path of the backbone and connectivity of the strand in the horizontal direction.
Likewise, the ‘|’ symbol is used to indicate vertical paths.
Bends can be drawn with ‘/’ and ‘\’, and will be parsed by the script into the correct bend based on the surrounding context.

The RNA sequence is indicated by a string of ‘N’s indicating nucleotide positions.

Base pair interactions can be specified by the ‘:’ symbol.

Here is a blueprint diagram for a simple hairpin:

===============================================================================
>My_Pattern                              
                                         
         5--NNNNNNNNNNN-NNNNNNNNNNN-NN-\ 
            ::::::::::: :::::::::::    | 
         3--NNNNNNNNNNN-NNNNNNNNNNN-NN-/ 
                                         
                                             
===============================================================================

To build this into a 3D PDB model, run RNAbuild.pl:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
perl RNAbuild.pl pattern.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This will generate a new PDB file called My_Pattern.pdb (as specified by the first line in the pattern.txt input file)

Positions in the pattern file that are not explicitly specified by a nucleotide (A,U,C,G) will be represented in the model by a random base-pair according to the specified constraints.

You can view the PDB file in any standard PDB viewer.

Designs that have more than 10,000 atoms are split into multiple smaller PDB files.


Edit the pattern.txt file to chance the unspecified positions to specified ones:

===============================================================================
>My_Pattern                                                                         
                                                                                    
         5--GGAACUCUGCG-CCCCCGAGUAG-UU-\                                                                           
            ::::::::::: :::::::::::    |                                                                        
         3--CCUUGAGAUGC-GGGGGCUUAUC-GC-/                                                                           
                                                                                    
                                                                                    
===============================================================================

Now run RNAbuild.pl again:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
perl RNAbuild.pl pattern.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Notice that the My_Pattern.pdb file now has the specified sequence.

The next program that will be useful to use is trace_analysis.pl, which reads and parses RNA blueprint inputs. This script is useful for quickly checking that a blueprint is valid. It also outputs a new blueprint with any plain-text symbols translated into nicer-looking UNICODE characters. The program is additionally useful for analyzing finalized RNA designs, as will be described later.

Run trace_analysis.pl in the command line on the pattern input:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
perl trace_analysis.pl pattern.txt > output.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This will read the pattern you drew in pattern.txt, and produce the following output saved as output.txt:


<contents of output.txt>
===============================================================================
My_Pattern
((((((((((((((((((((((....))))))))))))))))))))))
GGAACUCUGCGCCCCCGAGUAGUUCGCUAUUCGGGGGCGUAGAGUUCC




My Structure map:  
((((((((((((((((((((((....)))))))))))))))))))))) 

                                                                                     
                                                                                     
         5──GGAACUCUGCG─CCCCCGAGUAG─UU─╮                                                                            
            ┊┊┊┊┊┊┊┊┊┊┊ ┊┊┊┊┊┊┊┊┊┊┊    │                                                                         
         3──CCUUGAGAUGC─GGGGGCUUAUC─GC─╯                                                                            
                                                                                     
                                                                                     


Strand Path
                                                                                     
                                                                                     
         5─────────────────────────────╮                                                                            
            │││││││││││ │││││││││││    │                                                                         
         3─────────────────────────────╯                                                                            
                                                                                     
                                                                                     
Sequence good for pattern search..


Highlighting Repeat Sequences 
 8S|8W Repeats: 3	  5G|5C|5A|5U Repeats: 2   Restriction: 0   Complements: 18 
                                                                                     
                                                                                     
         5──◦◦◦◦◦◦◦◦SSS─CCCCCSP◦◦◦◦─◦◦─╮                                                                            
            │││││││││││ │││││││││││    │                                                                         
         3──◦◦◦◦◦◦◦◦◦SS─GGGGGSP◦◦◦◦─◦◦─╯                                                                            
                                                                                     
                                                                                     


Highlighting GU Pairs
                                                                                     
                                                                                     
         5──◦◦◦◦◦◦◦◦G◦◦─◦◦◦◦◦◦◦G◦◦◦─◦◦─╮                                                                            
            │││││││││││ │││││││││││    │                                                                         
         3──◦◦◦◦◦◦◦◦U◦◦─◦◦◦◦◦◦◦U◦◦◦─◦◦─╯                                                                            
                                                                          

Highlighting Structural Barriers

                                                                                     
                                                                                     
         5──◦◦◦◦◦◦◦◦◦◦◦─◦◦◦◦◦◦◦◦◦◦◦─◦◦─╮                                                                            
            │││││││││││ │││││││││││    │                                                                         
         3──◦◦◦◦◦◦◦◦◦◦◦─◦◦◦◦◦◦◦◦◦◦◦─◦◦─╯                                                                            
                                                                                     
===============================================================================

Observe that trace_analysis.pl has read the input blueprint and produced a series of different outputs that might be useful.

The first three lines of output contain the name, extended-alphabet dot-parens, and sequence constraints for the design - conveniently this is the input format that Revolvr will read.

Next, the "Structure map" contains just the normal dot-parens (the difference between normal and extended-alphabet dot-parens will be explored in more detail later, so stay tuned...)

Next, the program outputs a 'beautified' strand path that converts the standard text characters into UNICODE symbols.  These diagrams can be copy-pasted and used as inputs to all of the different scripts: RNAbuild, RNAvis and Revolvr.

For analyzing structures that contain sequences, trace_analysis.pl will also search for patterns and GU  pairs and highlight these on the strand path.  S/W mark runs of 8 or more Strong (C,G) or Weak (A,U) nts, respectively.  CCCCC, GGGGG, UUUUU, AAAAA mark places where 5 in a row are found.  ‘P’ marks where a sequence and a complementary Watson-Crick sequence of 10 bps are found.  Restriction sites that are found are marked ‘X’.  

The Structural Barriers section of the output will be described later, as it only applies to structures that contain at least 1 KL interaction - for a single hairpin no barriers are observed. (Structural barriers highlight portions of the 2D sequence that become nested beneath KL interactions by more than a 1/2 helix turn of stem and are thus potentially topologically trapped).
------------
Now let us focus on how to represent the structure of an RNA origami in blueprint format. For this example we will use the 2AE structure from our 2014 paper8, which creates a hexagonal tiling:

Diagraming this structure requires the use of two tertiary motifs, namely the KL interactions (of which there are two):

1) 180 KL; this requires a 9 nt loop sequence to be specified at the 2D level, and the pattern for the motif is 5'-AA-NNNNNN-A-3' such that only the middle 6 nts pair with Watson-Crick (WC) specificity.  RNAbuild and Revolvr are designed to specifically recognize the 9 nt loop as a 180KL.

2) 120KL; this requires a 7 nt loop sequence to be specified at the 2D level.  The pattern for 120KL is 5'-NNNNNNN-3', where all 7 nts pair with WC specificity.

On the 2D diagram, intramolecular 180KLs that provide internal links in the structure must have their WC pairings marked by ‘*’ rather than ‘:’, to specify a pseudoknot (PK).
Additionally, and very importantly, for each crossover position we must delineate the non-stacked part of the junction with a ‘^’.  The ‘^’ symbol is how the program will orient the antiparallel crossings correctly.

Loops coding for up to 10 different intermolecular KLs are encoded with an ‘X’.  The string of numbers and letters after the ‘@’ symbol on Line 2 defines the connectivity of intermolecular KL as they occur on the strand from 5' to 3', where the kth string of ‘X’es will later be automatically be replaced by the kth symbol in the string by trace_pattern.pl. Numbered loops pair with lettered loops to form KLs according to the rule: A-1, B-2, C-3, D-4, E-5, F-6, G-7, H-8, I-9, and J-0.

===============================================================================
>2AE                                                                        
@A1B2                                                                        
                                ╭──AA────╮                                  
╭XXXXXXX─NNNNKNNNNN──NN──NNNN─A╮╰NNNNNN─╮╰─NNNNNNN──NN──NNN────────╮        
│        ┊┊┊┊┊┊┊┊┊┊  ┊┊  ┊┊┊┊  │ ****** │  ┊┊┊┊┊┊┊  ┊┊  ┊┊┊        │        
╰────────NNNNKNNNNN╮╭NN──NNNN─╮╰─NNNNNN╮╰A─NNNNNNN╮╭NN──NNN─XXXXXXX╯        
                   ^╰───╮     ╰────AA──╯          ^╰───╮                    
                   ╰───╮^                         ╰───╮^                    
       ╭XXXXXXX─NNN──NN╯╰NNNNNKNNNNN─────NNNNNKNNN──NN╯╰NNNNKNNNNN────────╮ 
       │        ┊┊┊  ┊┊  ┊┊┊┊┊┊┊┊┊┊┊     ┊┊┊┊┊┊┊┊┊  ┊┊  ┊┊┊┊┊┊┊┊┊┊        │ 
       ╰────────NNN──NN──NNNNNKNNNNN3   5GGAANKNNN──NN──NNNNKNNNNN─XXXXXXX╯ 
                                                                            
                                                                             
===============================================================================
To build this into a 3D model, run RNAbuild:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
perl RNAbuild.pl pattern.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The script will generate 2AE.pdb

Now, to generate a sequence for the 2AE design, we will use Revolvr.
Revolvr requires two input files to run properly, the first file is called pattern.txt, which will hold the blueprint design, and the second file is called target.txt, that will encode the name, 2D structure and constraints. Importantly, these files must have exactly the names pattern.txt and target.txt.

You can generate target.txt from the pattern file by using trace_pattern.pl:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
perl trace_pattern.pl pattern.txt > target.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The script will analyze pattern.txt and write the output to target.txt

===============================================================================
2AE
(((((((((((((((((((((AAAAAAA))))))))))(((((1111111)))))(((((((..[[[[[[.))))))))))))))))))(((((((((((((((((..]]]]]].))))))((((((((((BBBBBBB))))))))))(((((2222222))))))))))))))))
GGAANKNNNNNNNNNKNNNNNNNNNNNNNNNNNKNNNNNNNNNNNNNNNNNNNNNNNNNNNNAANNNNNNANNNNNNNNNNNNKNNNNNNNNNNKNNNNNNNNNNNAANNNNNNANNNNNNNNNNNKNNNNNNNNNNNNNNNKNNNNNNNNNNNNNNNNNNNNNNNNNNNKNNNNN
===============================================================================
target.txt contains 3-4 lines that encode the following:
  Line 1) Name of the design
  Line 2) Dot-parens description of the structure
     In this dot-parens diagram, ‘( )’ indicate base-pairs and ‘.’ indicate single strand.
     ‘[ ]’ indicate intramolecular KL interactions.
   A-1, B-2, C-3... indicate pairings of intermolecular KLs marked ‘X’
  Line 3) Sequence constraints (N, AUCG, S, K, R, Y)
  Line 4) (Optional) Initial Sequence

After generating a target.txt it is possible to run Revolvr a single time using batch_revolvr.pl and specifying the number of desired output designs (1):
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
perl batch_revolvr.pl 1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Revolvr outputs a folder named data_2AE filled with design files.

The first file, 2AE-1_trace.txt will contain three things: a sequence design, a strand diagram filled with the sequence, and a 2D strand cartoon.

===============================================================================
GGAAAGCCAGCAGCGUCAGCGUUCUUCCCGCUGGCGCUAUACAGGAAGAAUGUAUGGCGAACAACGCGAAAGUUCGCCGCUGGUUUUCCAGGCCUUAAUCUACCGCAAUUCGCGAGCGGUAGUCGCUUUCAUGCCGAUUGAAGGCGACCUAAGAUCGGCACUUAGGAUUAGGGCCU

                                                                                                                                                                                                                                                                             
                                ╭──AA────╮                                             
╭UAGCCGU─ACUUUCGCUG──AU──GGCG─A╮╰CGCGAA─╮╰─CAAGCGG──UA──UGU────────╮                   
│        ┊┊┊┊┊┊┊┊┊┊  ┊┊  ┊┊┊┊  │ ****** │  ┊┊┊┊┊┊┊  ┊┊  ┊┊┊        │                   
╰────────UGAAGGCGAC╮╭UA──CCGC─╮╰─GCGCUU╮╰A─GUUCGCC╮╭AU──ACA─GGAAGAA╯                   
                   │╰───╮     ╰────AA──╯          │╰───╮                               
                   ╰───╮│                         ╰───╮│                               
       ╭ACGGCUA─GAA──UC╯╰CUAAUUCCGGA─────CCUUUUGGU──CG╯╰UCGCGGUCGC────────╮            
       │        ┊┊┊  ┊┊  ┊┊┊┊┊┊┊┊┊┊┊     ┊┊┊┊┊┊┊┊┊  ┊┊  ┊┊┊┊┊┊┊┊┊┊        │            
       ╰────────CUU──AG──GAUUAGGGCCU3   5GGAAAGCCA──GC──AGCGUCAGCG─UUCUUCC╯            
                                                                                       
                                                                                       
===============================================================================
    
The second file, 2AE-1_spool.txt contains a record of all mutations applied to the sequence to generate the output.

The third file, 2AE-1_design.txt contains a summary of the design, GC-content and other stats, as well as a survey of all the KL interactions and their energies.

To have the program generate more than one design, reset the target.txt and change the number of designs.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
perl trace_pattern.pl pattern.txt > target.txt
perl batch_revolvr.pl 10
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Revolvr will now put 10 designs into the folder it generates.
A summary file 2AE_summary.txt is also included, and it simply lists all the final design details for each design along with the ‘ensemble diversity’, which is the average base-pair distance between all structures in the Boltzmann ensemble.  Revolvr only outputs a design, if it finds a sequence which folds into the target as a minimum free energy (MFE) structure. However, the sequence it finds may also be consistent with a large number of other folds, which are not the free energy structure, but contribute to the thermodynamic ensemble. Thus even though the target structure might be the MFE structure for the sequence found, the sequence might spend very little time in that structure, instead visiting the other low-energy folds. The partition function folding performed by Vienna takes all of these structures into account, and so it is possible to rank designs based on how much time they spend in the MFE target structure rather than other folds. The lower the ensemble diversity, the better. 

Ensemble diversity is related to another commonly used measure of how well the ensemble approximates the target structure, the ‘ensemble defect’.  Ensemble defect captures the average number of incorrect bases of the thermodynamic ensemble, both bases which should be paired but aren’t and bases which aren’t paired but should be. The ensemble defect is often normalized by the length of the structure to give a ‘normalized ensemble defect’ (NED), the average percentage of nucleotides that are incorrectly paired at equilibrium relative to the target secondary structure. We have not added a facility to Revolvr to output NED, but structures can often be analyzed with NUPACK to calculate this measure. Where NUPACK’s energy model yields an identical MFE structure for the design, one can simply use the utilities function of NUPACK to calculate the NED. Where NUPACK’s energy model yields a slightly different MFE, the NED calculated by NUPACK utilities should be close to the actual NED for the targets structure. To calculate the NED perfectly (under NUPACK’s energy model), when the target structure differs significantly from NUPACK’s MFE, download the base pair probabilities from NUPACK’s utilities function and calculate the NED using the desired target structure. 

For the larger structures described in this paper (e.g. ZigZag-B-4x and Ribbon-5H-3X) Revolvr typically outputs sequences having a NED of ~5% but without having been specifically designed to the thermodynamic ensemble. NED for sequences found for smaller designs are more variable, ranging from 8 15%. Importantly, aptamers typically have empirically determined structures that do not correlate well with their MFE structures under any energy model and thus structures incorporating aptamers will often have overestimated ensemble defects. For such structures it is important, if you should desire an accurate NED, to perform a restricted calculation limited to the non-aptamer part of the target structure. Importantly, the measure Revolvr outputs—the ensemble diversity—does not suffer from this problem.

Run times for Revolvr depend on both the length of the RNA design, Tthe structural complexity and the speed of your processor. Runs for designs of 1000 nucleotides typically take 5-10 minutes; design might not even terminate in cases with particularly difficult constraints, or unfortunate random seeds. Structures up to 5000 nucleotides in length have been designed in less than 2 hours. Currently Revolvr is set to 2 hours, units of seconds with the command $timeout_length = 7200. 

########################################
Drawing RNA origami structure blueprints
########################################

In cases where you would like to specify Watson-Crick base pairing but the stem is only 1-2 base pairs long, then they should be marked with ‘!’.  This notation specifies pairs that are expected to NOT fold in the ViennaRNA energy model, but that are still required to have Watson-Crick complementarity by design. The designation will lock the sequence to be complementary while simultaneously trying to satisfy the ‘!’ unpaired designation in the 2D fold.  These are best used for 1 bp and 2 bp dovetails, as the shortest stem that can be formed in the ViennaRNA energy model is actually 3 bp long (thus there is a zero probability of forming a 1 bp or 2 bp stem in this energy model).  The ‘!’ designation is a way for us to logically design short stems irrespective of their folding in the energy model.

In the following example, the 2 bp dovetails have been designated by ‘!’ and are also specified to be G or C by designating them as ‘S’.
The 5' end has been constrained to have the sequence ‘GGAA’, which is known to improve the yield of RNA by T7 polymerase.  The sequence complementary to the 5' end is constrained to be ‘CC’ to avoid placing any GU pairs at the 5' end.
Additionally, this design includes two internal kissing loops (designated by‘*’) and two external kissing loops (marked ‘X’). The connectivity of the external loops is indicated on line 2 of the input ‘@A1’, which in this case says the first and second loops pair.
Lastly, two light-up aptamer motifs are added by putting the exact sequence and 2D structure into the diagram as corresponds to the pattern in the library file Motif_Library.txt

===============================================================================================
>Demo_Scaffold                                                                                 
@A1                                                                                                  
                             ╭──AA────╮                                                             
 ╭XXXXXXX─NNNNN───NNNNNNNN─A╮╰NNNNNN─╮╰─NNNNN───NNNNNNNN──CUG─UU─GA─GUAGAGUGUGGGCUCNNNNNNNNNGC╮  
 │        ┊┊┊┊┊   ┊┊┊┊┊┊┊┊  │ ****** │  ┊┊┊┊┊   ┊┊┊┊┊┊┊┊  ┊┊┊    ┊┊             ┊ ┊┊┊┊┊┊┊┊┊┊  │ 
 ╰────────NNNNN╮╭─NNNNNNNN─╮╰─NNNNNN╮╰A─NNNNN─╮╭NNNNNNNN──GAC────CU──GGGCUGG──GAGUGNNNNNNNNNUU╯ 
               ^╰────╮     ╰────AA──╯         ^╰────╮                        
               ╰────╮^         ╭──AA────╮     ╰────╮^                        
         ╭UUNNNN──SS╯╰NNNNNN─A╮╰NNNNNN─╮╰─NNNNN──SS╯╰NNNNGC╮                 
         │  ┊┊┊┊  !!  ┊┊┊┊┊┊  │ ****** │  ┊┊┊┊┊  !!  ┊┊┊┊  │                 
         ╰CGNNNN╮╭SS──NNNNNN─╮╰─NNNNNN╮╰A─NNNNN╮╭SS──NNNNUU╯                 
                ^╰────╮      ╰────AA──╯        ^╰────╮                       
                ╰────╮^                        ╰────╮^                       
       ╭XXXXXXX─NNNNN╯╰NNNNNNN───────CCNNNNNNNNNNNNN╯╰NNNNN─UGUAC───CC─────╮ 
       │        ┊┊┊┊┊  ┊┊┊┊┊┊┊       ┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊  ┊┊┊┊┊ ┊┊┊┊┊   ┊┊     │ 
       ╰────────NNNNN──NNNNNNN3     5GGAANNNNNNNNNNN──NNNNN─ACAUG-A-GG-AUCA╯ 
                                                                                                    
===============================================================================================

For quick reference, below are the four aptamer motifs used in this work, drawn with the 5' end entering at the bottom left. In Motif_Library.txt each of the four aptamers is drawn in all four possible orientations for convenience when using these aptamers as modules in different designs.                             

===============================================
                                               
Mango Aptamer                                  
                                                 
      3─CACG─AGA─GG─AGA─GG────╮                
        ┊┊┊┊                  │                
      5─GUGC─GAA─GG─GAC─GG─UGC╯                
                                                
iSpinach Aptamer                               
                                                
      3─CUCGGGUGUGAGAUG─AG─UU─GUC─5           
        ┊ ┊             ┊┊    ┊┊┊              
      5─GUGAG──GGUCGGG──UC────CAG─3          
                                                   
  MS2                                          
      3─UGUAC───CC─────╮                      
        ┊┊┊┊┊   ┊┊     │                     
      5─ACAUG─A─GG─AUCA╯                    
                                           
  PP7                                      
      3─CCGUG───CUUC───────╮                
        ┊┊┊┊┊   ┊┊┊┊       │               
      5─GGCAC─A─GAAG─AUAUGG╯               
                                             
===============================================

Now try to generate a different variation of the Demo_Scaffold by switching the aptamers with different ones from the above list. To Copy/Paste 2D patterns hold down the Option key to change the cursor to the cross (this works in TextEdit on Mac), and select a rectangle around the residues of the motif you wish to copy, select an insertion point and paste (or select a rectangle to replace text).


Sometimes, when working with complicated diagrams it can be helpful to be able to generate flipped copies of the diagram either horizontally or vertically. The program flip_trace.pl will read an input blueprint and produce an output with all four orientations of the pattern.  

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
perl flip_trace.pl pattern.txt > pattern_flipped.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This command will read pattern.txt and produce an output file called pattern_flipped.txt containing the four flipped variants.

Now, let's try to design a structure with more helices and see what issues arise:

Take the following 3-helix tall scaffold and, using copy/paste, extend it to be 4 helices tall.  You can do this by selecting the middle five rows and copy/pasting them.

===============================================================================================
>Demo_Scaffold
@AB21

                                              ╭──AA────╮                                              
 ╭XXXXXXX─NNNNNNNNNNNNNNNN──NNN──NNNNNNNN─NNA╮╰NNNNNN─╮╰NN─NNNNNNNNN──NNN──NNNNNNNNNNNNNNNNNN────────╮
 │        ┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊        │
 ╰────────NNNNNNNNNNNNNNNN╮╭NNN──NNNNNNNN─NN╮╰─NNNNNN╮╰ANN─NNNNNNNNN╮╭NNN──NNNNNNNNNNNNNNNNNN─XXXXXXX╯
                          ^╰────╮           ╰────AA──╯              ^╰────╮                             
                          ╰────╮^             ╭──AA────╮            ╰────╮^                           
                ╭UUNNNNNNN──NNN╯╰NNNNNNNN─NNA╮╰NNNNNN─╮╰NN─NNNNNNNNN──NNN╯╰NNNNNNNGC╮                 
                │  ┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊  │                 
                ╰CGNNNNNNN╮╭NNN──NNNNNNNN─NN╮╰─NNNNNN╮╰ANN─NNNNNNNNN╮╭NNN──NNNNNNNUU╯                 
                          ^╰────╮           ╰────AA──╯              ^╰────╮                           
                          ╰────╮^                                   ╰────╮^                           
     ╭XXXXXXX─NNNNNNNNNNNN──NNN╯╰NNNNNNNNNNNN─────NNNNNNNNNNNNNNNNNN──NNN╯╰NNNNNNNNNNN────────╮       
     │        ┊┊┊┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊┊┊┊┊     ┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊┊┊┊        │       
     ╰────────NNNNNNNNNNNN──NNN──NNNNNNNNNNNN3   5GGAANNNNNNNNNNNNNN──NNN──NNNNNNNNNNN─XXXXXXX╯       

===============================================================================================

Now, trace the new pattern and try to design a single sequence for it:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
perl trace_pattern.pl pattern.txt > target.txt
perl batch_revolvr.pl 1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For comparison, paste the following pattern.txt and see how the design time required for Revolvr changes. By pre-specifying certain base pairs to be wobbles (K-K) we reduce the amount of optimization required to remove complementary DNA sequences for the DNA genes.

Likewise, by specifying all the short 3 bp stems to be ‘S’ we can speed up the design time, since conveniently in the Vienna2 folding model the only stable 3 bp stems happen to be 100% GC.

===============================================================================================
>Demo_Scaffold_2
@AB21

                                              ╭──AA────╮                                                 
 ╭XXXXXXX─NNNNKNNNNNNKNNNN──NNN──NNNKNNNN─NNA╮╰NNNNNN─╮╰NN─NNNKNNNNN──NNN──NKNNNNNNKNNNNNKNNN────────╮           
 │        ┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊        │             
 ╰────────NNNNKNNNNNNKNNNN╮╭NNN──NNNKNNNN─NN╮╰─NNNNNN╮╰ANN─NNNKNNNNN╮╭NNN──NKNNNNNNKNNNNNKNNN─XXXXXXX╯           
                          ^╰────╮           ╰────AA──╯              ^╰────╮                                                                                                                                                                                                         
                          ╰────╮^             ╭──AA────╮            ╰────╮^                              
                ╭UUNNNKNNN──SSS╯╰NNNKNNNN─NNA╮╰NNNNNN─╮╰NN─NNNKNNNNN──SSS╯╰NNNKNNNGC╮                    
                │  ┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊  │                    
                ╰CGNNNKNNN╮╭SSS──NNNKNNNN─NN╮╰─NNNNNN╮╰ANN─NNNKNNNNN╮╭SSS──NNNKNNNUU╯                    
                          ^╰────╮           ╰────AA──╯              ^╰────╮                              
                          ╰────╮^             ╭──AA────╮            ╰────╮^                              
                ╭UUNNNKNNN──SSS╯╰NNNKNNNN─NNA╮╰NNNNNN─╮╰NN─NNNKNNNNN──SSS╯╰NNNKNNNGC╮                    
                │  ┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊  │                    
                ╰CGNNNKNNN╮╭SSS──NNNKNNNN─NN╮╰─NNNNNN╮╰ANN─NNNKNNNNN╮╭SSS──NNNKNNNUU╯                    
                          ^╰────╮           ╰────AA──╯              ^╰────╮                              
                          ╰────╮^                                   ╰────╮^                              
     ╭XXXXXXX─NNNKNNNNNNKN──NNN╯╰NNKNNNNKNNNN─────NNNNNKNNKNNNNNKNNN──NNN╯╰NNNNKNNNNNN────────╮          
     │        ┊┊┊┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊┊┊┊┊     ┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊┊┊┊        │          
     ╰────────NNNKNNNNNNKN──NNN──NNKNNNNKNNNN3   5GGAANKNNKNNNNNKNNN──NNN──NNNNKNNNNNN─XXXXXXX╯          

===============================================================================================


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
perl trace_pattern.pl pattern.txt > target.txt
perl batch_revolvr.pl 1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Try the above examples again with even more rows pasted in so they are more than 4 helices tall.  See how the design time increases as the length, number of KLs, and number of short dovetails that require stabilization increases.

Now, try to make the above example into a 2-column-wide design.  Do this by copy/pasting a box by holding down the option button as shown below.

First, draw the box around the core that will be copied, and then copy the section:

                          ╭─────────────────────────────────────────────────╮                                                       
                          |                    ╭──AA────╮                   |                              
 ╭XXXXXXX─NNNNKNNNNNNKNNNN|──NNN──NNNKNNNN─NNA╮╰NNNNNN─╮╰NN─NNNKNNNNN──NNN──|NKNNNNNNKNNNNNKNNN────────╮           
 │        ┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊|  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊  ┊┊┊  |┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊        │             
 ╰────────NNNNKNNNNNNKNNNN|╮╭NNN──NNNKNNNN─NN╮╰─NNNNNN╮╰ANN─NNNKNNNNN╮╭NNN──|NKNNNNNNKNNNNNKNNN─XXXXXXX╯           
                          |^╰────╮           ╰────AA──╯              ^╰────╮|                                                                                                                                                                                                         
                          |╰────╮^             ╭──AA────╮            ╰────╮^|                              
                ╭UUNNNKNNN|──SSS╯╰NNNKNNNN─NNA╮╰NNNNNN─╮╰NN─NNNKNNNNN──SSS╯╰|NNNKNNNGC╮                    
                │  ┊┊┊┊┊┊┊|  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊  ┊┊┊  |┊┊┊┊┊┊┊  │                    
                ╰CGNNNKNNN|╮╭SSS──NNNKNNNN─NN╮╰─NNNNNN╮╰ANN─NNNKNNNNN╮╭SSS──|NNNKNNNUU╯                    
                          |^╰────╮           ╰────AA──╯              ^╰────╮|                              
                          |╰────╮^             ╭──AA────╮            ╰────╮^|                              
                ╭UUNNNKNNN|──SSS╯╰NNNKNNNN─NNA╮╰NNNNNN─╮╰NN─NNNKNNNNN──SSS╯╰|NNNKNNNGC╮                    
                │  ┊┊┊┊┊┊┊|  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊  ┊┊┊  |┊┊┊┊┊┊┊  │                    
                ╰CGNNNKNNN|╮╭SSS──NNNKNNNN─NN╮╰─NNNNNN╮╰ANN─NNNKNNNNN╮╭SSS──|NNNKNNNUU╯                    
                          |^╰────╮           ╰────AA──╯              ^╰────╮|                              
                          |╰────╮^                                   ╰────╮^|                              
     ╭XXXXXXX─NNNKNNNNNNKN|──NNN╯╰NNKNNNNKNNNN─────NNNNNKNNKNNNNNKNNN──NNN╯╰|NNNNKNNNNNN────────╮          
     │        ┊┊┊┊┊┊┊┊┊┊┊┊|  ┊┊┊  ┊┊┊┊┊┊┊┊┊┊┊┊     ┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊  ┊┊┊  |┊┊┊┊┊┊┊┊┊┊┊        │          
     ╰────────NNNKNNNNNNKN|──NNN──NNKNNNNKNNNN3   5GGAANKNNKNNNNNKNNN──NNN──|NNNNKNNNNNN─XXXXXXX╯          
                          ╰─────────────────────────────────────────────────╯                                                        

Next, draw the box around the section that will be replaced/extended, and then paste the section over:

                                                                    ╭───────╮                                                       
                                              ╭──AA────╮            |       |                              
 ╭XXXXXXX─NNNNKNNNNNNKNNNN──NNN──NNNKNNNN─NNA╮╰NNNNNN─╮╰NN─NNNKNNNNN|──NNN──|NKNNNNNNKNNNNNKNNN────────╮           
 │        ┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊|  ┊┊┊  |┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊        │             
 ╰────────NNNNKNNNNNNKNNNN╮╭NNN──NNNKNNNN─NN╮╰─NNNNNN╮╰ANN─NNNKNNNNN|╮╭NNN──|NKNNNNNNKNNNNNKNNN─XXXXXXX╯           
                          ^╰────╮           ╰────AA──╯              |^╰────╮|                                                                                                                                                                                                         
                          ╰────╮^             ╭──AA────╮            |╰────╮^|                              
                ╭UUNNNKNNN──SSS╯╰NNNKNNNN─NNA╮╰NNNNNN─╮╰NN─NNNKNNNNN|──SSS╯╰|NNNKNNNGC╮                    
                │  ┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊|  ┊┊┊  |┊┊┊┊┊┊┊  │                    
                ╰CGNNNKNNN╮╭SSS──NNNKNNNN─NN╮╰─NNNNNN╮╰ANN─NNNKNNNNN|╮╭SSS──|NNNKNNNUU╯                    
                          ^╰────╮           ╰────AA──╯              |^╰────╮|                              
                          ╰────╮^             ╭──AA────╮            |╰────╮^|                              
                ╭UUNNNKNNN──SSS╯╰NNNKNNNN─NNA╮╰NNNNNN─╮╰NN─NNNKNNNNN|──SSS╯╰|NNNKNNNGC╮                    
                │  ┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊|  ┊┊┊  |┊┊┊┊┊┊┊  │                    
                ╰CGNNNKNNN╮╭SSS──NNNKNNNN─NN╮╰─NNNNNN╮╰ANN─NNNKNNNNN|╮╭SSS──|NNNKNNNUU╯                    
                          ^╰────╮           ╰────AA──╯              |^╰────╮|                              
                          ╰────╮^                                   |╰────╮^|                              
     ╭XXXXXXX─NNNKNNNNNNKN──NNN╯╰NNKNNNNKNNNN─────NNNNNKNNKNNNNNKNNN|──NNN╯╰|NNNNKNNNNNN────────╮          
     │        ┊┊┊┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊┊┊┊┊     ┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊|  ┊┊┊  |┊┊┊┊┊┊┊┊┊┊┊        │          
     ╰────────NNNKNNNNNNKN──NNN──NNKNNNNKNNNN3   5GGAANKNNKNNNNNKNNN|──NNN──|NNNNKNNNNNN─XXXXXXX╯          
                                                                    ╰───────╯                                                        


Now, repair the strand so that there is only one 5' and 3' end. There are two different but very similar options for where to put the ends, let’s choose the left side to repair and make continuous:


                                              ╭──AA────╮                                ╭──AA────╮                                                 
 ╭XXXXXXX─NNNNKNNNNNNKNNNN──NNN──NNNKNNNN─NNA╮╰NNNNNN─╮╰NN─NNNKNNNNN──NNN──NNNKNNNN─NNA╮╰NNNNNN─╮╰NN─NNNKNNNNN──NNN──NKNNNNNNKNNNNNKNNN────────╮           
 │        ┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊        │             
 ╰────────NNNNKNNNNNNKNNNN╮╭NNN──NNNKNNNN─NN╮╰─NNNNNN╮╰ANN─NNNKNNNNN╮╭NNN──NNNKNNNN─NN╮╰─NNNNNN╮╰ANN─NNNKNNNNN╮╭NNN──NKNNNNNNKNNNNNKNNN─XXXXXXX╯           
                          ^╰────╮           ╰────AA──╯              ^╰────╮           ╰────AA──╯              ^╰────╮                                                                                                                                                                                                         
                          ╰────╮^             ╭──AA────╮            ╰────╮^             ╭──AA────╮            ╰────╮^                              
                ╭UUNNNKNNN──SSS╯╰NNNKNNNN─NNA╮╰NNNNNN─╮╰NN─NNNKNNNNN──SSS╯╰NNNKNNNN─NNA╮╰NNNNNN─╮╰NN─NNNKNNNNN──SSS╯╰NNNKNNNGC╮                    
                │  ┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊  │                    
                ╰CGNNNKNNN╮╭SSS──NNNKNNNN─NN╮╰─NNNNNN╮╰ANN─NNNKNNNNN╮╭SSS──NNNKNNNN─NN╮╰─NNNNNN╮╰ANN─NNNKNNNNN╮╭SSS──NNNKNNNUU╯                    
                          ^╰────╮           ╰────AA──╯              ^╰────╮           ╰────AA──╯              ^╰────╮                              
                          ╰────╮^             ╭──AA────╮            ╰────╮^             ╭──AA────╮            ╰────╮^                              
                ╭UUNNNKNNN──SSS╯╰NNNKNNNN─NNA╮╰NNNNNN─╮╰NN─NNNKNNNNN──SSS╯╰NNNKNNNN─NNA╮╰NNNNNN─╮╰NN─NNNKNNNNN──SSS╯╰NNNKNNNGC╮                    
                │  ┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊  │                    
                ╰CGNNNKNNN╮╭SSS──NNNKNNNN─NN╮╰─NNNNNN╮╰ANN─NNNKNNNNN╮╭SSS──NNNKNNNN─NN╮╰─NNNNNN╮╰ANN─NNNKNNNNN╮╭SSS──NNNKNNNUU╯                    
                          ^╰────╮           ╰────AA──╯              ^╰────╮           ╰────AA──╯              ^╰────╮                              
                          ╰────╮^                                   ╰────╮^                                   ╰────╮^                              
     ╭XXXXXXX─NNNKNNNNNNKN──NNN╯╰NNKNNNNKNNNN─────NNNNNKNNKNNNNNKNNN──NNN╯╰NNKNNNNKNNNN─────NNNNNKNNKNNNNNKNNN──NNN╯╰NNNNKNNNNNN────────╮          
     │        ┊┊┊┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊┊┊┊┊     ┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊┊┊┊┊     ┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊┊┊┊        │          
     ╰────────NNNKNNNNNNKN──NNN──NNKNNNNKNNNN3   5GGAANKNNKNNNNNKNNN──NNN──NNKNNNNKNNNN3   5GGAANKNNKNNNNNKNNN──NNN──NNNNKNNNNNN─XXXXXXX╯          



Finally, paste the following pattern into pattern.txt and examine it with trace_analysis.pl:



===============================================================================================================================================================
>Demo_Scaffold                                                                                                                                  
@AB21                                                                                                                                           
                                              ╭──AA────╮                                ╭──AA────╮                                              
 ╭XXXXXXX─NNNNKNNNNNNKNNNN──NNN──NNNKNNNN─NNA╮╰NNNNNN─╮╰NN─NNNKNNNNN──NNN──NNNKNNNN─NNA╮╰NNNNNN─╮╰NN─NNNKNNNNN──NNN──NKNNNNNNKNNNNNKNNN────────╮
 │        ┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊        │
 ╰────────NNNNKNNNNNNKNNNN╮╭NNN──NNNKNNNN─NN╮╰─NNNNNN╮╰ANN─NNNKNNNNN╮╭NNN──NNNKNNNN─NN╮╰─NNNNNN╮╰ANN─NNNKNNNNN╮╭NNN──NKNNNNNNKNNNNNKNNN─XXXXXXX╯
                          ^╰────╮           ╰────AA──╯              ^╰────╮           ╰────AA──╯              ^╰────╮                           
                          ╰────╮^             ╭──AA────╮            ╰────╮^             ╭──AA────╮            ╰────╮^                           
                ╭UUNNNKNNN──SSS╯╰NNNKNNNN─NNA╮╰NNNNNN─╮╰NN─NNNKNNNNN──SSS╯╰NNNKNNNN─NNA╮╰NNNNNN─╮╰NN─NNNKNNNNN──SSS╯╰NNNKNNNGC╮                 
                │  ┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊  │                 
                ╰CGNNNKNNN╮╭SSS──NNNKNNNN─NN╮╰─NNNNNN╮╰ANN─NNNKNNNNN╮╭SSS──NNNKNNNN─NN╮╰─NNNNNN╮╰ANN─NNNKNNNNN╮╭SSS──NNNKNNNUU╯                 
                          ^╰────╮           ╰────AA──╯              ^╰────╮           ╰────AA──╯              ^╰────╮                           
                          ╰────╮^             ╭──AA────╮            ╰────╮^             ╭──AA────╮            ╰────╮^                           
                ╭UUNNNKNNN──SSS╯╰NNNKNNNN─NNA╮╰NNNNNN─╮╰NN─NNNKNNNNN──SSS╯╰NNNKNNNN─NNA╮╰NNNNNN─╮╰NN─NNNKNNNNN──SSS╯╰NNNKNNNGC╮                 
                │  ┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊  │                 
                ╰CGNNNKNNN╮╭SSS──NNNKNNNN─NN╮╰─NNNNNN╮╰ANN─NNNKNNNNN╮╭SSS──NNNKNNNN─NN╮╰─NNNNNN╮╰ANN─NNNKNNNNN╮╭SSS──NNNKNNNUU╯                 
                          ^╰────╮           ╰────AA──╯              ^╰────╮           ╰────AA──╯              ^╰────╮                           
                          ╰────╮^                                   ╰────╮^                                   ╰────╮^                           
     ╭XXXXXXX─NNNKNNNNNNKN──NNN╯╰NNKNNNNKNNNN─────NNNNNKNNKNNNNNKNNN──NNN╯╰NNKNNNNKNNNN─────NNNNNKNNKNNNNNKNNN──NNN╯╰NNNNKNNNNNN────────╮       
     │        ┊┊┊┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊┊┊┊┊     ┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊┊┊┊┊     ┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊┊┊┊        │       
     ╰────────NNNKNNNNNNKN──NNN──NNKNNNNKNNNN─────NNNNNKNNKNNNNNKNNN──NNN──NNKNNNNKNNNN3   5GGAANKNNKNNNNNKNNN──NNN──NNNNKNNNNNN─XXXXXXX╯       
                                                                                                                                                ===============================================================================================================================================================


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
perl trace_analysis.pl pattern.txt > pattern_analysis.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Notice that at the end of the file there is an additional output called Structural Barriers.  This portion contains a strand-path analysis, specifically it looks at closing ‘)’ partner of pairs that contain a closing KL loop ‘]’ between itself and its partner ‘(’. In these cases, if the KL forms very fast there is a chance the KL could form and associate before the closing pair of the stem (especially if the two parts of the stem are transcriptionally far apart). In such situations, ‘X’ marks any regions that are more than 6 bps deep (or a half-turn) long where a topological clash may occur.
 

===============================================================================================================================================================

Highlighting Structural Barriers
                                                                                                                                                         
                                              ╭──◦◦────╮                                ╭──◦◦────╮                                                       
 ╭◦◦◦◦◦◦◦─◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦──XXX──XXXX~~~~─~~◦╮╰◦◦◦◦◦◦─╮╰◦◦─◦◦◦◦◦◦◦◦◦──XXX──XXXX~~~~─~~◦╮╰◦◦◦◦◦◦─╮╰◦◦─◦◦◦◦◦◦◦◦◦──◦◦◦──◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦────────╮              
 │        ││││││││││││││││  │││  ││││││││ ││ │ ││││││ │ ││ │││││││││  │││  ││││││││ ││ │ ││││││ │ ││ │││││││││  │││  ││││││││││││││││││        │                
 ╰────────◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦╮╭~~~──~~~~~~~~─~~╮╰─◦◦◦◦◦◦╮╰◦◦◦─◦◦◦◦◦◦◦◦◦╮╭~~~──~~~~~~~~─~~╮╰─◦◦◦◦◦◦╮╰◦◦◦─◦◦◦◦◦◦◦◦◦╮╭◦◦◦──◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦─◦◦◦◦◦◦◦╯              
                          │╰────╮           ╰────◦◦──╯              │╰────╮           ╰────◦◦──╯              │╰────╮                                                                                                                                                                                                            
                          ╰────╮│             ╭──◦◦────╮            ╰────╮│             ╭──◦◦────╮            ╰────╮│                                    
                ╭◦◦◦◦◦◦◦◦◦──~~~╯╰XXXX~~~~─~~◦╮╰◦◦◦◦◦◦─╮╰◦◦─◦◦◦◦◦◦◦◦◦──~~~╯╰XXXX~~~~─~~◦╮╰◦◦◦◦◦◦─╮╰◦◦─◦◦◦◦◦◦◦◦◦──◦◦◦╯╰◦◦◦◦◦◦◦◦◦╮                          
                │  │││││││  │││  ││││││││ ││ │ ││││││ │ ││ │││││││││  │││  ││││││││ ││ │ ││││││ │ ││ │││││││││  │││  │││││││  │                          
                ╰◦◦◦◦◦◦◦◦◦╮╭~~~──~~~~~~~~─~~╮╰─◦◦◦◦◦◦╮╰◦◦◦─◦◦◦◦◦◦◦◦◦╮╭~~~──~~~~~~~~─~~╮╰─◦◦◦◦◦◦╮╰◦◦◦─◦◦◦◦◦◦◦◦◦╮╭◦◦◦──◦◦◦◦◦◦◦◦◦╯                          
                          │╰────╮           ╰────◦◦──╯              │╰────╮           ╰────◦◦──╯              │╰────╮                                    
                          ╰────╮│             ╭──◦◦────╮            ╰────╮│             ╭──◦◦────╮            ╰────╮│                                    
                ╭◦◦◦◦◦◦◦◦◦──~~~╯╰XXXX~~~~─~~◦╮╰◦◦◦◦◦◦─╮╰◦◦─◦◦◦◦◦◦◦◦◦──~~~╯╰XXXX~~~~─~~◦╮╰◦◦◦◦◦◦─╮╰◦◦─◦◦◦◦◦◦◦◦◦──◦◦◦╯╰◦◦◦◦◦◦◦◦◦╮                          
                │  │││││││  │││  ││││││││ ││ │ ││││││ │ ││ │││││││││  │││  ││││││││ ││ │ ││││││ │ ││ │││││││││  │││  │││││││  │                          
                ╰◦◦◦◦◦◦◦◦◦╮╭~~~──~~~~~~~~─~~╮╰─◦◦◦◦◦◦╮╰◦◦◦─◦◦◦◦◦◦◦◦◦╮╭~~~──~~~~~~~~─~~╮╰─◦◦◦◦◦◦╮╰◦◦◦─◦◦◦◦◦◦◦◦◦╮╭◦◦◦──◦◦◦◦◦◦◦◦◦╯                          
                          │╰────╮           ╰────◦◦──╯              │╰────╮           ╰────◦◦──╯              │╰────╮                                    
                          ╰────╮│                                   ╰────╮│                                   ╰────╮│                                    
     ╭◦◦◦◦◦◦◦─◦◦◦◦◦◦◦◦◦◦◦◦──◦◦◦╯╰~~~~~~~~~~~~─────~~~~~~~~~~~~~~~~~~──~~~╯╰~~~~~~~~~~~~─────◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦──◦◦◦╯╰◦◦◦◦◦◦◦◦◦◦◦────────╮                
     │        ││││││││││││  │││  ││││││││││││     ││││││││││││││││││  │││  ││││││││││││     ││││││││││││││││││  │││  │││││││││││        │                
     ╰────────◦◦◦◦◦◦◦◦◦◦◦◦──◦◦◦──~~~~~~XXXXXX─────XXXXXXXXXXXXXXXXXX──XXX──XXXXXXXXXXXX3   5◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦──◦◦◦──◦◦◦◦◦◦◦◦◦◦◦─◦◦◦◦◦◦◦╯                
                                                                                                                                                         
================================================================================================================================================================

In the above example, the short stretches of ‘X’ to the left of each KL are likely not long enough to be a significant problem for folding.  However, at the 3' tail there is a quite long region of ‘X’es that indicates this might be a folding issue.  Indeed, in the case where all the KLs were formed and the 3' end not folded, that strand would have to thread through a narrow space 3 times!

Now, try moving the position of the 3' and 5' ends of the above scaffold to another position and run trace_analysis again.  The trace_analysis tool is often useful when evaluating and comparing different variants of an RNA design.


Next, try using copy/paste to swap the KL at the top right with the section in the bottom row where 3' and 5' ends meet. There are many different possibilities for routing the strand. You can even move the locations of crossovers or remove them entirely, as long as everything is connected together into one strand:

=================================================================================================================================================================
>Demo_Scaffold  
@AB12
                                              ╭──AA────╮  
 ╭XXXXXXX─NNNNKNNNNNNKNNNN──NNN──NNNKNNNN─NNA╮╰NNNNNN─╮╰NN─NNNKNNNNN──NNN──NNKNNNNKNNNN─────NNNNNKNNKNNNNNKNNN──NNN──NKNNNNNNKNNNNNKNNN────────╮ 
 │        ┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊┊┊┊┊     ┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊        │ 
 ╰────────NNNNKNNNNNNKNNNN╮╭NNN──NNNKNNNN─NN╮╰─NNNNNN╮╰ANN─NNNKNNNNN╮╭NNN──NNKNNNNKNNNN3   5GGAANKNNKNNNNNKNNN╮╭NNN──NKNNNNNNKNNNNNKNNN─XXXXXXX╯ 
                          ^╰────╮           ╰────AA──╯              ^╰────╮                                   ^╰────╮          
                          ╰────╮^                                   ╰────╮^             ╭──AA────╮            ╰────╮^           
                ╭UUNNNKNNN──SSS╯╰NNKNNNNKNNNN─────NNNNNKNNKNNNNNKNNN──SSS╯╰NNNKNNNN─NNA╮╰NNNNNN─╮╰NN─NNNKNNNNN──SSS╯╰NNNKNNNGC╮ 
                │  ┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊┊┊┊┊     ┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊  │ 
                ╰CGNNNKNNN╮╭SSS──NNKNNNNKNNNN─────NNNNNKNNKNNNNNKNNN╮╭SSS──NNNKNNNN─NN╮╰─NNNNNN╮╰ANN─NNNKNNNNN──SSS──NNNKNNNUU╯ 
                          ^╰────╮                                   ^╰────╮           ╰────AA──╯                                
                          ╰────╮^             ╭──AA────╮            ╰────╮^                                                     
                ╭UUNNNKNNN──SSS╯╰NNNKNNNN─NNA╮╰NNNNNN─╮╰NN─NNNKNNNNN──SSS╯╰NNKNNNNKNNNN─────NNNNNKNNKNNNNNKNNN──SSS──NNNKNNNGC╮ 
                │  ┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊┊┊┊┊     ┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊  │ 
                ╰CGNNNKNNN╮╭SSS──NNNKNNNN─NN╮╰─NNNNNN╮╰ANN─NNNKNNNNN╮╭SSS──NNKNNNNKNNNN─────NNNNNKNNKNNNNNKNNN╮╭SSS──NNNKNNNUU╯ 
                          ^╰────╮           ╰────AA──╯              ^╰────╮                                   ^╰────╮           
                          ╰────╮^             ╭──AA────╮            ╰────╮^             ╭──AA────╮            ╰────╮^           
     ╭XXXXXXX─NNNKNNNNNNKN──NNN╯╰NNNKNNNN─NNA╮╰NNNNNN─╮╰NN─NNNKNNNNN──NNN╯╰NNNKNNNN─NNA╮╰NNNNNN─╮╰NN─NNNKNNNNN──NNN╯╰NNNNKNNNNNN────────╮  
     │        ┊┊┊┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊ ┊┊ │ ****** │ ┊┊ ┊┊┊┊┊┊┊┊┊  ┊┊┊  ┊┊┊┊┊┊┊┊┊┊┊        │  
     ╰────────NNNKNNNNNNKN──NNN──NNNKNNNN─NN╮╰─NNNNNN╮╰ANN─NNNKNNNNN──NNN──NNNKNNNN─NN╮╰─NNNNNN╮╰ANN─NNNKNNNNN──NNN──NNNNKNNNNNN─XXXXXXX╯  
                                            ╰────AA──╯                                ╰────AA──╯                                                                                                                                                            
=================================================================================================================================================================
Typically, making the strand-path more complex and convoluted will increase the number of problem regions marked ‘X’ by trace_analysis. (We discuss this more in the section on RNAvis which marks these problem regions in a 3D animation.)
Repeat the analysis on a few different designs and compare the outputs.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
perl trace_analysis.pl pattern.txt > pattern_analysis.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############################
Troubleshooting RNA Designs
############################

Below are some common problems that may cause Revolvr to fail to finish a design.

•Short helices cause Revolvr to fail to solve the initial fold. 
Solution - Lengthen short helices.

• Designs are almost completely solved, diverging from the target only near junctions. 
Solution - The structures are too highly branched. Add offset dovetails, such as +11/ 11 bp dovetails, to decrease the branching number. 

Consider another example: a weak 2 bp dovetail between two crossovers essentially behaves like a 6-way multi loop that can transiently form 2 bps in its middle. We desire that the dovetail be stable. Increasing the 2 bp dovetail to a full stem of 13 makes the 6 way loop into a much less dynamic pair of 3-way branches, which is much easier to design.

•Designs with a motif having weak or noncanonical structure (e.g. an aptamer) do not complete. 
Solution - there are not enough ‘N’ residues surrounding the motif to stabilize it. Increase the length of the helix adjacent to the motif.

•Designs with duplicated motifs having constrained sequence do not complete because Revolvr 
cannot design adjacent sequence without creating longer repeats. 
Solution  - to each instance of the motif, add an additional one or more base pairs that differentiate it from all other instances of that motif. 

#########################################################################
Drawing more complex strand diagrams for Revolvr (Incompatible with RNAbuild/RNAvis)
#########################################################################

IMPORTANT: RNAbuild and RNAvis will NOT be able to model any of the following non-RNA origami designs in 3D, and will produce strange errors and/or structures if attempted.  However, Revolvr will attempt sequence design for any input that it can parse, allowing users to create designs that are not within the limited scope of RNAbuild’s motif-library. 

So far we have described blueprints for which the strand path is simple and does not cross itself. But what about cases where we cannot draw the structure easily without some lines crossing? In such cases we can use the ‘+’ symbol to indicate places where two strands cross.

See the following example, a small nuclear RNA (snoRNA) with crossing strands, drawn in plain text for simplicity. Notice that the slanted lines ‘/’ and ‘\’ are interpreted as going 5' to 3' or 5' to 3' based on context. To accommodate crossing, some helices are drawn in the vertical direction: base pairs in vertically oriented helices are marked differently, with a ‘=’.





====================================================================================
>SnoRNA_shaped_example                                                                                                                         
                                                                                
       /N\        /N\                                                                         
       N N        N N                                                                      
       N N        N N                                                                     
       N=N        N=N                                                                      
       N=N        N=N                                                                     
       N=N        N=N                                                                      
       N=N        N=N                                                                      
     /-/ \-\    /-/ \-\                                                                        
     | /N\ |    | /N\ |                                                                           
     N*N N*N    N*N N*N                                                                        
     N*N N*N    N*N N*N                                                                        
     N*N N*N    N*N N*N                                                                       
     N*N N*N    N*N N*N                                                                        
     N*N N*N    N*N N*N                                                                  
 3-N-+-/ \-+-NN-+-/ \-+-\                                                                                
     \-\ /-/    \-\ /-/ |                                                                     
       N=N        N=N   N                                                                  
       N=N        N=N   N                                                                     
       N=N        N=N   N                                                                     
       N=N        N=N   N                                                                    
       N=N        N=N   N                                                                    
   5-NN/ \NNNNNNNN/ \NNN/                                                                           
                                                                                
====================================================================================

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
perl trace_pattern.pl pattern.txt > target.txt
perl batch_revolvr.pl 1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Example output from Revolvr design
Please note that Revolvr is not intended to design snoRNAs, this is simply an example of complex and arbitrary design with vertical helices, crossing strands, and pseudoknotted Watson-Crick base pairs. In this example the target MFE in the Vienna model is the two hairpins followed by an unfolded single-strand; the contributions of pseudoknot interactions to the folding energy are not considered (i.e., they are only logically designed based on indicated base-pairs).

====================================================================================
                                                                                                                                      
       ╭U╮        ╭A╮                                                                                                                 
       A G        A U                                                                                                                 
       A A        U U                                                                                                                 
       C=G        A=U                                                                                                                 
       G=C        G=C                                                                                                                 
       C=G        A=U                                                                                                                 
       A=U        C=G                                                                                                                 
     ╭─╯ ╰─╮    ╭─╯ ╰─╮                                                                                                               
     │ ╭A╮ │    │ ╭A╮ │                                                                                                               
     A*U U*A    A*U U*A                                                                                                               
     G*C U*A    U*A C*G                                                                                                               
     A*U G*C    A*U U*A                                                                                                               
     C*G U*A    G*C U*A                                                                                                               
     U*A A*U    G*C U*A                                                                                                               
 3─C─┼─╯ ╰─┼─AA─┼─╯ ╰─┼─╮                                                                                                             
     ╰─╮ ╭─╯    ╰─╮ ╭─╯ │                                                                                                             
       U=A        G=C   A                                                                                                             
       U=A        A=U   A                                                                                                             
       C=G        A=U   A                                                                                                             
       G=C        G=C   A                                                                                                             
       U=A        A=U   C                                                                                                             
   5─AA╯ ╰AAAUACCA╯ ╰AAA╯                                                                                                             
                                                                                                                                      
====================================================================================

Also note that in the snoRNA example above, if the 5 bp pseudoknotted sections are extended to 6 bp and 7 bp, this will cause Revolvr to treat each pseudoknot as either a 180KL or 120KL and it will attempt to energetically orthogonalize them (this is by design as Revolvr is intended to make RNA origami with primarily those 2 motifs!). Also notice that these types of designs are many times more difficult to design than standard RNA origami structures! Several attempts running the script may be required to get even a single folding result with such difficult inputs.











Below is another example that is quite difficult for Revolvr to design. It arose during as a design challenge in the EteRNA game (https://eternagame.org/home/) and it has no biological significance. It showcases the ability of blueprints to represent highly unusual and branched structures, as well as the ability of Revolvr to design them. Short stems of only 2-3 bps are very weak, and thus it is often difficult to find a sequence that will stabilize them. Keep this in mind when designing structures: the larger number of short helices, the more difficult it will be to find a sequence that folds that way in the energy model (and structures with such features are also very unlikely to fold properly in the LAB!).

====================================================================================
>Not_a_knot                                                                                      
                                                                                                                                                                                                                     
             ╭NNNNNNNNN╮     ╭NNNNNNNNN╮                                                                                                          
             N  ┊┊┊┊┊  N     N  ┊┊┊┊┊  N                                                                                                          
             N ╭NNNNN╮ N     N ╭NNNNN╮ N                                                                                                          
             N=N     N=N     N=N     N=N                                                                                                  
             N=N     N=N     N=N     N=N                                                                                                  
             N=N     N N     N=N     N=N                                                                                                  
             N=N     N N     N=N     N=N                                                                                                  
             N=N     ╰─╯     N=N     N=N                                                                                                        
     ╭NNNNNN╮N ╰NNNNNNN─NNNNN╯ │╭NNNN╯ N                                                                                               
     N  ┊┊  │N  ┊┊┊┊┊   ┊┊┊┊┊  N│  ┊┊  N                                                                                         
     N ╭NNNN╯│ ╭NNNNN╮ ╭NNNNN╮ N╰NNNNNN╯                                                                                         
     N=N     N=N     N=N     N=N                                                                                                                                                         
     N=N     N=N     N=N     N=N                                                                                                                                                         
     N=N     N N     N N     N=N                                                                                                 
     N=N     N N     N N     N=N                                                                                                 
     N=N     ╰─╯     ╰─╯     N=N                                                                                                                                                                                                                                         
     N ╰NNNNN─NNNNNNNNN─NNNNN╯ ╰NNNNNNN╮                                                                                                                                                                                             
     N  ┊┊┊┊┊   ┊┊┊┊┊   ┊┊┊┊┊   ┊┊┊┊┊  N                                                                                                                                                                                         
     ╰NNNNNNN╮ ╭NNNNN─NNNNNNNNN─NNNNN╮ N                                                                                         
             N=N             ╭─╮     N=N                                                                                         
             N=N     5 3     N N     N=N                                                                                         
             N=N     N=N     N N     N=N                                                                                         
             N=N     N=N     N=N     N=N                                                                                         
             N=N     N=N     N=N     N=N                                                                                         
     ╭NNNNNN╮N ╰NNNNN╯ ╰NNNN╮│ ╰NNNNN╯ N                                                                                         
     N  ┊┊  │N  ┊┊┊┊┊   ┊┊  │N  ┊┊┊┊┊  N                                                                                                                                                                                                                                            
     N ╭NNNN╯│ ╭NNNNNNN─NNNN╯N ╭NNNNNNN╯                                                                                                                                                                                                                                               
     N=N     N=N     ╭─╮     N=N                                                                                                                                                                                          
     N=N     N=N     N N     N=N                                                                                                 
     N=N     N=N     N N     N=N                                                                                                 
     N=N     N=N     N=N     N=N                                                                                                 
     N=N     N=N     N=N     N=N                                                                                                 
     N ╰NNNNN╯ N     N ╰NNNNN╯ N                                                                                                 
     N  ┊┊┊┊┊  N     N  ┊┊┊┊┊  N                                                                                                 
     ╰NNNNNNNNN╯     ╰NNNNNNNNN╯                                                                                                 
                                                                                                                             
====================================================================================
Note that the above blueprint is drawn using the fancier UNICODE strand path characters, because it was possible to draw a much more compact blueprint, and easier to understand strand path using UNICODE rather than ASCII. This is often the case for complex designs.








Sample Revolvr design:

====================================================================================
GCCCAGGAUCCUUGUAAUGACUAAUAAUGAAAAGCGACAUUGGCAGAAGCAAUCUACGUGUACCACAAGGAAAAUACACCGACGUAGACAAUGAAAAUCGCUCGUAUUUCAUCCGGAACUGAAGAAGGAUACCCCAAAAGGAAAAGGUAUACAACCUUCGAAAUUCAGACCGAAAGAAUGCUGAAUGUGCCCAUAUGGAAAAGCACAAACAUUCAGAAGCGCUAGCCAUUCGGAAAACCUCGGUCCGGAAUUAGAAAUUACUAAACUCGUAAGCGAACUCGAAAACUUACAAUAGAGUUAAGGAAAAGAAGUCCCUCAGGCCCCUGUGGAAAAGCCUGGUGCAGGGAAACACUUCUUCCUGAAGCAAAAGCGGC

                                                                                                  
                                                                                                                                                                                                                                                                                                                       
             ╭AAACACGAA╮     ╭AGCUUCCAA╮                                                                                                                                                                                                            
             C  ┊┊┊┊┊  A     A  ┊┊┊┊┊  C                                                                                                                                                                                                            
             A ╭UGUGC╮ A     A ╭GAAGG╮ A                                                                                                                                                                                                            
             U=A     C=G     U=A     A=U                                                                                                                                                                                                    
             U=A     C=G     U=A     U=A                                                                                                                                                                                                    
             C=G     A U     C=G     A=U                                                                                                                                                                                                    
             A=U     U A     A=U     C=G                                                                                                                                                                                                    
             G=C     ╰─╯     G=C     C=G                                                                                                                                                                                                          
     ╭AAGCUC╮A ╰GUAAGAA─AGCCA╯ │╭AACC╯ A                                                                                                                                                                                                 
     A  ┊┊  │A  ┊┊┊┊┊   ┊┊┊┊┊  A│  ┊┊  A                                                                                                                                                                                           
     A ╭CGAA╯│ ╭CAUUC╮ ╭UCGGU╮ A╰AAGGAA╯                                                                                                                                                                                           
     C=G     G=C     G=C     C=G                                                                                                                                                                                                                                                           
     U=A     C=G     G=C     C=G                                                                                                                                                                                                                                                           
     U=A     G A     A A     G=C                                                                                                                                                                                                   
     A=U     C U     A A     G=C                                                                                                                                                                                                   
     C=G     ╰─╯     ╰─╯     A=U                                                                                                                                                                                                                                                                                                                                           
     A ╰CUCAA─AUCAUUAAA─GAUUA╯ ╰ACUUUAU╮                                                                                                                                                                                                                                                                                               
     A  ┊┊┊┊┊   ┊┊┊┊┊   ┊┊┊┊┊   ┊┊┊┊┊  G                                                                                                                                                                                                                                                                                           
     ╰UAGAGUU╮ ╭GUAAU─GACUAAUAA─UGAAA╮ C                                                                                                                                                                                           
             A=U             ╭─╮     A=U                                                                                                                                                                                           
             A=U     5 3     A G     G=C                                                                                                                                                                                           
             G=C     G=C     A A     C=G                                                                                                                                                                                           
             G=C     C=G     G=C     G=C                                                                                                                                                                                           
             A=U     C=G     C=G     A=U                                                                                                                                                                                           
     ╭AAGGUG╮A ╰AGGAC╯ ╰CGAA╮│ ╰GUUAC╯ A                                                                                                                                                                                           
     A  ┊┊  │A  ┊┊┊┊┊   ┊┊  │A  ┊┊┊┊┊  A                                                                                                                                                                                                                                                                                                                                              
     A ╭CCCU╯│ ╭UCCUGAA─GCAA╯A ╭CAAUGAA╯                                                                                                                                                                                                                                                                                                                                                 
     G=C     A=U     ╭─╮     U=A                                                                                                                                                                                                                                                                                            
     C=G     G=C     A C     C=G                                                                                                                                                                                                   
     C=G     A=U     A A     U=A                                                                                                                                                                                                   
     U=A     A=U     G=C     A=U                                                                                                                                                                                                   
     G=C     G=C     G=C     C=G                                                                                                                                                                                                   
     G ╰UCCCU╯ A     A ╰AUGUG╯ C                                                                                                                                                                                                   
     U  ┊┊┊┊┊  C     A  ┊┊┊┊┊  A                                                                                                                                                                                                   
     ╰GCAGGGAAA╯     ╰AAUACACCG╯                                                                                                                                                                                                   
                                                                                                                                                                                                                                  
====================================================================================



########################################################################
How to create folding movie visualizations with RNAvis
########################################################################

The purpose of RNAvis is to help structure designers visualize the gradual build up of an RNA origami structure as it is transcribed and folds. 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
perl RNAvis.pl pattern.txt $synthesis_step_size $KL_delay
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The software generates a series of PDB models visualizing the elongation and folding of a nascent RNA strand as it forms the secondary structure in the blueprint pattern.txt. The arguments $synthesis_step_size and $KL_delay are optional. One keyframe model is produced for each step, with the step-size between keyframes ($synthesis_step_size) set to 15 nts as a default. By default, to save on file size, the variable $pmode is set to 1. With this setting, RNAvis produces PDB files that only contain the phosphate atom of each nucleotide position. To produce full-atom keyframes instead, set $pmode to 0 within the file RNAvis.pl.

Outputs are stored in a folder, which contains all of the keyframes and an analysis file containing a summary of the substructures generated. Along with these files, a Chimera script command file marked_commands.cmd is created, which allows the animation of keyframes and morphing between them to be set up and rendered in Chimera v1.10, automatically. Note that these scripts do not work correctly in later versions of Chimera.

In this work we chose to synthesize structures that were predicted to have no topological clashes or at most a couple less-serious ones, such as ZigZag-B-4X (Supplementary Movie 1, zero clashes) or Ribbon-9H (Supplementary Movie 2, two clashes, each involving no more than 1.5 turns of helix), for the particular delay ($KL_delay) between the transcription time and KL folding, which is set to 300 nt by default. 

However, it is illuminating to consider a single shape (which we did not synthesize) with the same phosphate backbone design but in two versions: one that has a significant number of predicted topological clashes, and one that does not. In general, the relative position of crossovers, dovetails and KLs all affect the order of folding events and propensity for a design to have a putative topological clash. Such effects are often nonobvious and unintuitive. Here, we will see that simply changing the position of the 5' and 3' ends of a design will have a large impact on the viability of a design for cotranscriptional folding. Below we give the trace_analysis output of a 9-helix, two-column design with a highly serpentine folding pathway, call it 9H-9topo:

=================================================================================================================================
Highlighting Structural Barriers

  Plot of delay = 300 nts before closing KL interactions. 

                                                                                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                                                                                        
     ╭◦◦◦◦◦◦◦◦◦◦◦──XX──XXXXXXXXXXXXXXXX─────XXXXXXXXXXXXXXX──XX──XXXXXXXXXXXXXXXX─────XXXXXXXXX~~~~~~───◦◦◦◦◦◦◦◦◦◦◦─◦◦◦◦◦◦◦◦◦◦◦╮                                                                                                                                                                                                                                        
     │  │││││││││  ││  ││││││││││││││││     │││││││││││││││  ││  ││││││││││││││││     │││││││││││││││   │││││││││││ │││││││││  │                                                                                                                                                                                                                                        
     ╰◦◦◦◦◦◦◦◦◦◦◦╮╭~~──~~~~~~~~~~~~~~~~─────~~~~~~~~~~~~~~~──~~──~~~~~~~~~~~~~~~~─────~~~~~~~~~~~~~~~╮╭─◦◦◦◦◦◦◦◦◦◦◦─◦◦◦◦◦◦◦◦◦◦◦╯                                                                                                                                                                                                                                        
                 ^╰───╮                                                                              ^╰───╮                                                                                                                                                                                                                                                                                                                                                                                
                 ╰───╮^                     ╭──◦◦────╮                                ╭──◦◦────╮     ╰───╮^                                                                                                                                                                                                                                                             
         ╭◦◦◦◦◦◦◦◦◦◦◦╯╰~~~~~~~~~───◦◦◦◦◦◦◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦───◦◦◦◦◦◦◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦◦╯╰XXX~~~~~~──◦◦◦◦◦◦◦◦◦◦◦╮                                                                                                                                                                                                                                      
         │  │││││││││  │││││││││   │││││││ │ ││││││ │ │ │││││││  │││││││││   │││││││ │ ││││││ │ │ │││││││  │││││││││  │││││││││  │                                                                                                                                                                                                                                      
         ╰◦◦◦◦◦◦◦◦◦◦◦──~~~~~~XXX╮╭─◦◦◦◦◦◦◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦╮╭─◦◦◦◦◦◦◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦◦──~~~~~~~~~╮╭◦◦◦◦◦◦◦◦◦◦◦╯                                                                                                                                                                                                                                      
                                ^╰───╮    ╰────◦◦──╯                      ^╰───╮    ╰────◦◦──╯                      ^╰───╮                                                                                                                                                                                                                                              
                                ╰───╮^         ╭──◦◦────╮                 ╰───╮^         ╭──◦◦────╮                 ╰───╮^                                                                                                                                                                                                                                              
           ╭◦◦◦◦◦◦◦◦◦◦◦──XXXXXXXXXXX╯╰◦◦◦◦◦─◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦◦◦╯╰◦◦◦◦◦─◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦──~~~~~~~~~~~╯╰◦◦◦◦◦◦◦◦◦◦◦╮                                                                                                                                                                                                                                  
           │  │││││││││  │││││││││││  │││││ │ │ ││││││ │ │ ││││││  │││││││││││  │││││ │ │ ││││││ │ │ ││││││  │││││││││││  │││││││││  │                                                                                                                                                                                                                                  
           ╰◦◦◦◦◦◦◦◦◦◦◦╮╭~~~~~~~~~~~──◦◦◦◦◦─◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦╮╭◦◦◦◦◦◦◦◦◦◦◦──◦◦◦◦◦─◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦╮╭XXXXXXXXXXX──◦◦◦◦◦◦◦◦◦◦◦╯                                                                                                                                                                                                                                  
                       ^╰───╮                ╰────◦◦──╯          ^╰───╮                ╰────◦◦──╯          ^╰───╮                                                                                                                                                                                                                                                       
                       ╰───╮^                     ╭──◦◦────╮     ╰───╮^                     ╭──◦◦────╮     ╰───╮^                                                                                                                                                                                                                                                       
               ╭◦◦◦◦◦◦◦◦◦◦◦╯╰~~~~~~~~~───◦◦◦◦◦◦◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦◦╯╰◦◦◦◦◦◦◦◦◦───◦◦◦◦◦◦◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦◦╯╰XXX~~~~~~──◦◦◦◦◦◦◦◦◦◦◦╮                                                                                                                                                                                                                                
               │  │││││││││  │││││││││   │││││││ │ ││││││ │ │ │││││││  │││││││││   │││││││ │ ││││││ │ │ │││││││  │││││││││  │││││││││  │                                                                                                                                                                                                                                
               ╰◦◦◦◦◦◦◦◦◦◦◦──~~~~~~XXX╮╭─◦◦◦◦◦◦◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦╮╭─◦◦◦◦◦◦◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦◦──~~~~~~~~~╮╭◦◦◦◦◦◦◦◦◦◦◦╯                                                                                                                                                                                                                                
                                      ^╰───╮    ╰────◦◦──╯                      ^╰───╮    ╰────◦◦──╯                      ^╰───╮                                                                                                                                                                                                                                        
                                      ╰───╮^         ╭──◦◦────╮                 ╰───╮^         ╭──◦◦────╮                 ╰───╮^                                                                                                                                                                                                                                        
                 ╭◦◦◦◦◦◦◦◦◦◦◦──XXXXXXXXXXX╯╰◦◦◦◦◦─◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦◦◦╯╰◦◦◦◦◦─◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦◦◦╯╰◦◦◦◦◦◦◦◦◦◦◦╮                                                                                                                                                                                                                            
                 │  │││││││││  │││││││││││  │││││ │ │ ││││││ │ │ ││││││  │││││││││││  │││││ │ │ ││││││ │ │ ││││││  │││││││││││  │││││││││  │                                                                                                                                                                                                                            
                 ╰◦◦◦◦◦◦◦◦◦◦◦╮╭~~~~~~~~~~~──◦◦◦◦◦─◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦╮╭◦◦◦◦◦◦◦◦◦◦◦──◦◦◦◦◦─◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦╮╭◦◦◦◦◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦◦◦╯                                                                                                                                                                                                                            
                             ^╰───╮                ╰────◦◦──╯          ^╰───╮                ╰────◦◦──╯          ^╰───╮                                                                                                                                                                                                                                                 
                             ╰───╮^                     ╭──◦◦────╮     ╰───╮^                     ╭──◦◦────╮     ╰───╮^                                                                                                                                                                                                                                                 
                     ╭◦◦◦◦◦◦◦◦◦◦◦╯╰~~~~~~~~~───◦◦◦◦◦◦◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦◦╯╰◦◦◦◦◦◦◦◦◦5 3◦◦◦◦◦◦◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦◦╯╰◦◦◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦◦◦╮                                                                                                                                                                                                                          
                     │  │││││││││  │││││││││   │││││││ │ ││││││ │ │ │││││││  │││││││││   │││││││ │ ││││││ │ │ │││││││  │││││││││  │││││││││  │                                                                                                                                                                                                                          
                     ╰◦◦◦◦◦◦◦◦◦◦◦──~~~~~~XXX╮╭─◦◦◦◦◦◦◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦╮╭─◦◦◦◦◦◦◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦╮╭◦◦◦◦◦◦◦◦◦◦◦╯                                                                                                                                                                                                                          
                                            ^╰───╮    ╰────◦◦──╯                      ^╰───╮    ╰────◦◦──╯                      ^╰───╮                                                                                                                                                                                                                                  
                                            ╰───╮^         ╭──◦◦────╮                 ╰───╮^         ╭──◦◦────╮                 ╰───╮^                                                                                                                                                                                                                                  
                       ╭◦◦◦◦◦◦◦◦◦◦◦──XXXXXXXXXXX╯╰◦◦◦◦◦─◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦──~~~~~~~~~~~╯╰◦◦◦◦◦─◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦◦◦╯╰◦◦◦◦◦◦◦◦◦◦◦╮                                                                                                                                                                                                                      
                       │  │││││││││  │││││││││││  │││││ │ │ ││││││ │ │ ││││││  │││││││││││  │││││ │ │ ││││││ │ │ ││││││  │││││││││││  │││││││││  │                                                                                                                                                                                                                      
                       ╰◦◦◦◦◦◦◦◦◦◦◦╮╭~~~~~~~~~~~──◦◦◦◦◦─◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦╮╭XXXXXXXXXXX──◦◦◦◦◦─◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦╮╭◦◦◦◦◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦◦◦╯                                                                                                                                                                                                                      
                                   ^╰───╮                ╰────◦◦──╯          ^╰───╮                ╰────◦◦──╯          ^╰───╮                                                                                                                                                                                                                                           
                                   ╰───╮^                     ╭──◦◦────╮     ╰───╮^                     ╭──◦◦────╮     ╰───╮^                                                                                                                                                                                                                                           
                           ╭◦◦◦◦◦◦◦◦◦◦◦╯╰~~~~~~~~~───◦◦◦◦◦◦◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦◦╯╰XXX~~~~~~───◦◦◦◦◦◦◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦◦╯╰◦◦◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦◦◦╮                                                                                                                                                                                                                    
                           │  │││││││││  │││││││││   │││││││ │ ││││││ │ │ │││││││  │││││││││   │││││││ │ ││││││ │ │ │││││││  │││││││││  │││││││││  │                                                                                                                                                                                                                    
                           ╰◦◦◦◦◦◦◦◦◦◦◦──~~~~~~XXX╮╭─◦◦◦◦◦◦◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦◦──~~~~~~~~~╮╭─◦◦◦◦◦◦◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦╮╭◦◦◦◦◦◦◦◦◦◦◦╯                                                                                                                                                                                                                    
                                                  ^╰───╮    ╰────◦◦──╯                      ^╰───╮    ╰────◦◦──╯                      ^╰───╮                                                                                                                                                                                                                            
                                                  ╰───╮^                                    ╰───╮^                ╭──◦◦────╮          ╰───╮^                                                                                                                                                                                                                            
                              ╭◦◦◦◦◦◦◦◦◦◦◦─◦◦◦◦◦◦◦◦◦◦◦╯╰~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~───────~~╯╰◦◦◦◦◦◦◦◦◦◦◦─◦◦◦╮╰◦◦◦◦◦◦─╮╰◦◦─◦◦◦◦◦◦◦──◦◦╯╰◦◦◦◦◦◦◦◦◦◦◦╮                                                                                                                                                                                                                
                              │  │││││││││ │││││││││││  │││││││││││││││││││││││││││││││       ││  │││││││││││ ││ │ ││││││ │ ││ │││││││  ││  │││││││││  │                                                                                                                                                                                                                
                              ╰◦◦◦◦◦◦◦◦◦◦◦─◦◦◦◦◦◦◦◦◦◦◦──~~~~~~XXXXXXXXXXXXXXXXXXXXXXXXX───────XX──◦◦◦◦◦◦◦◦◦◦◦─◦◦╮╰─◦◦◦◦◦◦╮╰◦◦◦─◦◦◦◦◦◦◦──◦◦──◦◦◦◦◦◦◦◦◦◦◦╯                                                                                                                                                                                                                
                                                                                                                ╰────◦◦──╯                                                                                                                                                                                                                                              
=================================================================================================================================

A pre-rendered folding video of this design is provided in Supplementary Movie 3.

The 5' and 3' ends of 9H-9topo are located roughly in the middle of the design. As the structure folds, many loops are created (orange in movie, ‘~’ in diagram) that the nascent strand (red in movie, ‘X’ diagram) has to thread through. Nine potential structural barriers due to topological clashes are highlighted on the strand path diagram, and animated in the movie as well. For the topological clash in the top helix (the bottom helix as presented in Supplementary Movie 3), the nascent strand must wind through a loop in previously formed structure six times for the affected section of duplex to form correctly. 

As a side note: trace_analysis also has the variable $KL_delay, but it must be set by changing the variable within trace_analysis rather than as an argument at the command line. For both trace_analysis and RNAvis setting $KL_delay = 0 will cause prediction of a topological clash for almost every KL.  Setting $KL_delay to be equal or greater to the length of a design will result in no topological clashes being predicted. Try this for 9H-9topo, if you desire.

Now, with $KL_delay = 300, consider the following diagram, in which we carefully relocated the 5' and 3' end to reduce the number of potential topological clashes to create 9H-0topo:

=================================================================================================================================

Highlighting Structural Barriers

  Plot of delay = 300 nts before closing KL interactions. 

                                                                                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                                                                                        
     ╭◦◦◦◦◦◦◦◦◦◦◦──◦◦──◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦─────◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦──◦◦──◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦─────◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦───◦◦◦◦◦◦◦◦◦◦◦─◦◦◦◦◦◦◦◦◦◦◦╮                                                                                                                                                                                                                                        
     │  │││││││││  ││  ││││││││││││││││     │││││││││││││││  ││  ││││││││││││││││     │││││││││││││││   │││││││││││ │││││││││  │                                                                                                                                                                                                                                        
     ╰◦◦◦◦◦◦◦◦◦◦◦╮╭◦◦──◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦─────◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦──◦◦──◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦─────◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦╮╭─◦◦◦◦◦◦◦◦◦◦◦─◦◦◦◦◦◦◦◦◦◦◦╯                                                                                                                                                                                                                                        
                 ^╰───╮                                                                              ^╰───╮                                                                                                                                                                                                                                                                                                                                                                                  
                 ╰───╮^                     ╭──◦◦────╮                                ╭──◦◦────╮     ╰───╮^                                                                                                                                                                                                                                                             
         ╭◦◦◦◦◦◦◦◦◦◦◦╯╰◦◦◦◦◦◦◦◦◦───◦◦◦◦◦◦◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦───◦◦◦◦◦◦◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦◦╯╰◦◦◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦◦◦╮                                                                                                                                                                                                                                      
         │  │││││││││  │││││││││   │││││││ │ ││││││ │ │ │││││││  │││││││││   │││││││ │ ││││││ │ │ │││││││  │││││││││  │││││││││  │                                                                                                                                                                                                                                      
         ╰◦◦◦◦◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦╮╭─◦◦◦◦◦◦◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦╮╭─◦◦◦◦◦◦◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦╮╭◦◦◦◦◦◦◦◦◦◦◦╯                                                                                                                                                                                                                                      
                                ^╰───╮    ╰────◦◦──╯                      ^╰───╮    ╰────◦◦──╯                      ^╰───╮                                                                                                                                                                                                                                              
                                ╰───╮^         ╭──◦◦────╮                 ╰───╮^         ╭──◦◦────╮                 ╰───╮^                                                                                                                                                                                                                                              
           ╭◦◦◦◦◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦◦◦╯╰◦◦◦◦◦─◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦◦◦╯╰◦◦◦◦◦─◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦◦◦╯╰─◦◦◦◦◦◦◦◦◦◦◦╮                                                                                                                                                                                                                                 
           │  │││││││││  │││││││││││  │││││ │ │ ││││││ │ │ ││││││  │││││││││││  │││││ │ │ ││││││ │ │ ││││││  │││││││││││   │││││││││  │                                                                                                                                                                                                                                 
           ╰◦◦◦◦◦◦◦◦◦◦◦╮╭◦◦◦◦◦◦◦◦◦◦◦──◦◦◦◦◦─◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦╮╭◦◦◦◦◦◦◦◦◦◦◦──◦◦◦◦◦─◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦╮╭◦◦◦◦◦◦◦◦◦◦◦───◦◦◦◦◦◦◦◦◦◦◦╯                                                                                                                                                                                                                                 
                       ^╰───╮                ╰────◦◦──╯          ^╰───╮                ╰────◦◦──╯          ^╰───╮                                                                                                                                                                                                                                                       
                       ╰───╮^                     ╭──◦◦────╮     ╰───╮^                     ╭──◦◦────╮     ╰───╮^                                                                                                                                                                                                                                                       
               ╭◦◦◦◦◦◦◦◦◦◦◦╯╰◦◦◦◦◦◦◦◦◦───◦◦◦◦◦◦◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦◦╯╰◦◦◦◦◦◦◦◦◦───◦◦◦◦◦◦◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦◦╯╰◦◦◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦◦◦╮                                                                                                                                                                                                                                
               │  │││││││││  │││││││││   │││││││ │ ││││││ │ │ │││││││  │││││││││   │││││││ │ ││││││ │ │ │││││││  │││││││││  │││││││││  │                                                                                                                                                                                                                                
               ╰◦◦◦◦◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦╮╭─◦◦◦◦◦◦◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦╮╭─◦◦◦◦◦◦◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦╮╭◦◦◦◦◦◦◦◦◦◦◦╯                                                                                                                                                                                                                                
                                      ^╰───╮    ╰────◦◦──╯                      ^╰───╮    ╰────◦◦──╯                      ^╰───╮                                                                                                                                                                                                                                        
                                      ╰───╮^         ╭──◦◦────╮                 ╰───╮^         ╭──◦◦────╮                 ╰───╮^                                                                                                                                                                                                                                        
                 ╭◦◦◦◦◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦◦◦╯╰◦◦◦◦◦─◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦◦◦╯╰◦◦◦◦◦─◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦◦◦╯╰◦◦◦◦◦◦◦◦◦◦◦╮                                                                                                                                                                                                                            
                 │  │││││││││  │││││││││││  │││││ │ │ ││││││ │ │ ││││││  │││││││││││  │││││ │ │ ││││││ │ │ ││││││  │││││││││││  │││││││││  │                                                                                                                                                                                                                            
                 ╰◦◦◦◦◦◦◦◦◦◦◦╮╭◦◦◦◦◦◦◦◦◦◦◦──◦◦◦◦◦─◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦╮╭◦◦◦◦◦◦◦◦◦◦◦──◦◦◦◦◦─◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦╮╭◦◦◦◦◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦◦◦╯                                                                                                                                                                                                                            
                             ^╰───╮                ╰────◦◦──╯          ^╰───╮                ╰────◦◦──╯          ^╰───╮                                                                                                                                                                                                                                                 
                             ╰───╮^                     ╭──◦◦────╮     ╰───╮^                     ╭──◦◦────╮     ╰───╮^                                                                                                                                                                                                                                                 
                     ╭◦◦◦◦◦◦◦◦◦◦◦╯╰◦◦◦◦◦◦◦◦◦───◦◦◦◦◦◦◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦◦╯╰◦◦◦◦◦◦◦◦◦───◦◦◦◦◦◦◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦◦╯╰◦◦◦◦◦◦◦◦◦5 3◦◦◦◦◦◦◦◦◦◦◦╮                                                                                                                                                                                                                         
                     │  │││││││││  │││││││││   │││││││ │ ││││││ │ │ │││││││  │││││││││   │││││││ │ ││││││ │ │ │││││││  │││││││││   │││││││││  │                                                                                                                                                                                                                         
                     ╰◦◦◦◦◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦╮╭─◦◦◦◦◦◦◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦╮╭─◦◦◦◦◦◦◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦╮╭─◦◦◦◦◦◦◦◦◦◦◦╯                                                                                                                                                                                                                         
                                            ^╰───╮    ╰────◦◦──╯                      ^╰───╮    ╰────◦◦──╯                      ^╰───╮                                                                                                                                                                                                                                  
                                            ╰───╮^         ╭──◦◦────╮                 ╰───╮^         ╭──◦◦────╮                 ╰───╮^                                                                                                                                                                                                                                  
                       ╭◦◦◦◦◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦◦◦╯╰◦◦◦◦◦─◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦◦◦╯╰◦◦◦◦◦─◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦◦◦╯╰◦◦◦◦◦◦◦◦◦◦◦╮                                                                                                                                                                                                                      
                       │  │││││││││  │││││││││││  │││││ │ │ ││││││ │ │ ││││││  │││││││││││  │││││ │ │ ││││││ │ │ ││││││  │││││││││││  │││││││││  │                                                                                                                                                                                                                      
                       ╰◦◦◦◦◦◦◦◦◦◦◦╮╭◦◦◦◦◦◦◦◦◦◦◦──◦◦◦◦◦─◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦╮╭◦◦◦◦◦◦◦◦◦◦◦──◦◦◦◦◦─◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦╮╭◦◦◦◦◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦◦◦╯                                                                                                                                                                                                                      
                                   ^╰───╮                ╰────◦◦──╯          ^╰───╮                ╰────◦◦──╯          ^╰───╮                                                                                                                                                                                                                                           
                                   ╰───╮^                     ╭──◦◦────╮     ╰───╮^                     ╭──◦◦────╮     ╰───╮^                                                                                                                                                                                                                                           
                           ╭◦◦◦◦◦◦◦◦◦◦◦╯╰◦◦◦◦◦◦◦◦◦───◦◦◦◦◦◦◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦◦╯╰◦◦◦◦◦◦◦◦◦───◦◦◦◦◦◦◦◦╮╰◦◦◦◦◦◦─╮╰◦─◦◦◦◦◦◦◦╯╰◦◦◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦◦◦╮                                                                                                                                                                                                                    
                           │  │││││││││  │││││││││   │││││││ │ ││││││ │ │ │││││││  │││││││││   │││││││ │ ││││││ │ │ │││││││  │││││││││  │││││││││  │                                                                                                                                                                                                                    
                           ╰◦◦◦◦◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦╮╭─◦◦◦◦◦◦◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦╮╭─◦◦◦◦◦◦◦╮╰─◦◦◦◦◦◦╮╰◦◦─◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦╮╭◦◦◦◦◦◦◦◦◦◦◦╯                                                                                                                                                                                                                    
                                                  ^╰───╮    ╰────◦◦──╯                      ^╰───╮    ╰────◦◦──╯                      ^╰───╮                                                                                                                                                                                                                            
                                                  ╰───╮^                                    ╰───╮^                ╭──◦◦────╮          ╰───╮^                                                                                                                                                                                                                            
                              ╭◦◦◦◦◦◦◦◦◦◦◦─◦◦◦◦◦◦◦◦◦◦◦╯╰◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦───────◦◦╯╰◦◦◦◦◦◦◦◦◦◦◦─◦◦◦╮╰◦◦◦◦◦◦─╮╰◦◦─◦◦◦◦◦◦◦──◦◦╯╰◦◦◦◦◦◦◦◦◦◦◦╮                                                                                                                                                                                                                
                              │  │││││││││ │││││││││││  │││││││││││││││││││││││││││││││       ││  │││││││││││ ││ │ ││││││ │ ││ │││││││  ││  │││││││││  │                                                                                                                                                                                                                
                              ╰◦◦◦◦◦◦◦◦◦◦◦─◦◦◦◦◦◦◦◦◦◦◦──◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦───────◦◦──◦◦◦◦◦◦◦◦◦◦◦─◦◦╮╰─◦◦◦◦◦◦╮╰◦◦◦─◦◦◦◦◦◦◦──◦◦──◦◦◦◦◦◦◦◦◦◦◦╯                                                                                                                                                                                                                
                                                                                                                ╰────◦◦──╯                                                                                                                                                                                                                                              
=================================================================================================================================


By positioning the 5' and 3' ends at the right edge (the fourth helix from the bottom) of the design, the 9H-0topo has no predicted topological clashes according to trace_analysis. A pre rendered folding video of this design is provided in Supplementary Movie 4.

Based only topological concerns, 9H-0topo appears to be a good candidate to attempt cotranscriptional folding. However other features of the design are of potential concern with respect to cotranscriptional folding. Compare Supplementary Movie 1 (ZigZag-B-4X) to Supplementary Movie 4 (for 9H-0topo). Both have no topological defects, but the folding of ZigZag-B-4X appears to proceed in a more orderly fashion. After the initial synthesis of a long single-stranded region (118 nt, almost 11 turns), ZigZag-B-4X begins to exhibit a single, compact, folded core that gets larger as more columns of KLs are built and sequentially double strand the long single-strand. On the other hand, 9H-0topo begins with five flexibly arranged hairpins and a 6-turn single-stranded region. The single-stranded region does not become double-stranded until more than halfway through the synthesis and the flexibly arranged hairpins do not get locked in place until the very end of synthesis. Thus one might imagine that 9H-0topo could be more prone to misfolding or aggregation than say ZigZag B 4X. This and many other features of RNA origami folding can now be examined in RNAvis, and hypotheses about them can be generated and experimentally tested in future studies.  

We note that, by design, RNAvis adds a randomly generated wobble to single-stranded regions and represents double-stranded regions as rigidified. This wobbling is in part of what makes synthesis and folding of 9H-0topo look more complicated than that of ZigZag B 4X. The extent of the wobbling can be controlled using the $eccentricity variable in the RNAvis script.


 


