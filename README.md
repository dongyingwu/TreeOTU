TreeOTU
=======

TreeOTU: Operational Taxonomic Unit Classification Based on Phylogenetic Trees  

Dongying Wu, Ladan Doroud, Jonathan A. Eisen

University of California, Davis, Davis, California 95616, USA


TreeOTU is a package of perl scripts that classify Operational Taxonomic Unit (OTU) based on phylogenetic trees. 
The package including the following components:

<b>1. Main script: TreeOTU.pl</b>

   TreeOTU.pl takes a rooted phylogenetic tree and a PN (Position of Node) cutoff (0<=PN<1) and outputs a set of OTUs
   
   Usage: TreeOTU.pl -i input_tree -c PN_cutoff -o output_file

    Input format:
     -i inputtree, the input tree must be in newick format
     -c PN cutoff, PN cutoff must >=0 and <1
    
    Output format:
     -o output file, if no output file is define, the script outputs to STDOUT
    
     The format of the output is demonstrated by the following example:

     input_tree_PN_cutoff_0.05  A,B  C D,E
     

     The first column is the OTU set description including tree file name and PN cutoff, 
     from the second column to the end of the line are the OTUs. 
     different OTUs are separated by spaces while the taxa in each OTU are separated by commas. 
    
<b>2. OTU comparison scripts</b>

   <b>A. Compare OTU sets based on Adjusted Mutual Information (AMI): OTUcompare_AMI.pl</b>
   
   OTUcompare_AMI.pl compares one of OTU set (Query) with multiple sets of OTUs (Target), 
   and outputs the AMI values (Ajusted Mutual Information) between the query and target(s). 
   AMI is a number between 0 and 1. AMI value of 1 indicates two OTU sets are identical, 
   the more different the two OTU sets, the smaller AMI gets. 
   
   Usage: OTUcompare_AMI.pl -query query_OTUset -target target_OTUset -shuffle 100
   
    Input format:
   
     -query query_OTUset, this input is a file with only ONE line of OTU set, 
                         the format of the input is exactly the same as the output of TreeOTU.pl
                        
     -target target_OTUset, this input target file can contain more than one OTU sets. 
                            One OTU set takes one line, the format of the OTU layout is  
                            exactly the same as the output of TreeOTU.pl  
   
     -shuffle an_integer, the number of shuffle times to determine the expect value in AMI calculation, 
                          the default value is 100
   
    Output format:
   
     The script outputs results to STDOUT, the output format is examplified by the following line:
   
     AMI=0.5696	Query_ID=BA00028_0.080	Target_ID=concat38_0.010	Shuffle=100	Query_taxa_count=4573	Target_taxa_count=4436	Shared_taxa_count=4424
 
     The output explains itself. Only the shared taxa in both OTU sets are included in the AMI calculation.


   <b>B. Identify the maximum F1 score of one OTU in a OTU set: OTUcompare_Fmeasure.pl</b>
   
   This script compares one OTU (query) against all the OTUs in a OTU dataset (target), and report the maximum F measure score.

   Usage: OTUcompare_Fmeasure.pl -query query_OTU -target target_OTUsets -shuffle 10 -beta 1
   
     Input format:

      -query query_OTU, the OTU group of interest, format: A,B,C,D

      -target target_OTUsets, the target OTU set that the query is compared against. The input target file can contain 
                              more than one OTU sets. One OTU set takes one line and the format of the OTU layout is
                              exactly the same as the output of TreeOTU.pl
      -shuffle integer, target random shuffling times to evaluate the expected maximum F measure score (the default is 10)

      -beta: F measure beta factor, the default is 1

     Output format:
    
      The script outputs results to STDOUT, and the output takes the format of the following line: 

       targetID=concat38_0.600	Adjusted_max_Fscore=0.0032	Max_Fscore=0.0640	Query_taxa_count=135	Query_taxa_evaluated=132	BETA=1	Shuffle=5	Expected_maxF=0.0610


<b>3. Tree rooting scripts</b>

   <b>A. Mid-point rooting: MidpointRooting.pl</b> 

   Usage: MidpointRotting.pl -i input_tree -o output_tree (input must be in newick format)

   <b>B. Rooting a tree between two nodes defined by six taxa: SixPointRooting.pl</b>
 
   The script roots the tree in the middle of two nodes. Three taxa define one node, thus the script needs 6 taxa to carry
   out the rooting process.
   
   Usage: SixPointRooting.pl -i input_tree -l taxa_list -o output_tree 
   
      input format: -i input tree, the tree must be in newick format
                    -l input taxa list, two lines, each line includes 3 taxa separated by spaces 
                    
   <b>C. Rooting the tree according to an input outgroup: OutgroupRooting.pl</b>

   The script takes a list of taxa, and tries to find the best matching monophyletic clade as the outgroup, 
   eventually roots the tree in the middle of the edge that connects the outgroup and the rest.
   
   Usage: OutgroupRooting.pl -i input_tree -l input_taxa -o output_tree
   
       input format: -i input tree, the tree must be in newick format
                     -l input list of taxa as the outgroup, can be non-monophyletic in the tree. The format of this input
                        is multiple lines with one aceesion in each line.
       output format: -o the rooted tree. There is a log file with the extension ".OutgroupRootingLog" that documents the 
                         actual outgroup the script decides to use, and its similarity with the input outgroup (F1 score, 
                         1 means identical, the F1 score should be >=0.5 for the script to carry out the rooting)
                         
                         
<b>4. Additional Data</b>                        
                  
      The file TreeOTU_ali_tree_taxa.tgz inlcudes alignments, trees and taxonomy information for 40 PhyEco marker families 
      and ssu-rRNAs from the IMG database as well as the ssu-rRNAs from the "All-Species Living Tree" project. 
                     
   

