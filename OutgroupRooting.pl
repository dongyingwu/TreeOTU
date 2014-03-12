#!/usr/bin/perl

use strict;
##usage: OutgroupRooting.pl -i tree -l taxa_list -o output 
my %opt=@ARGV;
my $opt_i=$opt{-i}; ##input tree file
my $opt_o=$opt{-o}; ##output tree file
my $opt_l=$opt{-l}; ##list of taxa as outgroup
my $logfile=$opt_o.".OutgroupRootingLog";

##take in the tree
my $in;
my $tree;
open(I,$opt_i) || die "cannot open input tree file \n";
while(<I>){$in.=$_;}
close I;

$tree=TREE_MOD->new($in);

##make sure the acceesions in the taxa_list are all in the tree
my @acc;
my $acc;
open(I,$opt_l) || die "cannot open taxa list input\n";
while(<I>){
my @t=split(/\s+/);
@acc=(@acc,@t);
}
close I;

foreach $acc(@acc){
if($tree->{tx_id}->{$acc}){next;}
die "$acc is not in the tree\n";
}


##if @acc only have 1 taxon, find it's parent and root in between
## or if the @acc leave out only one taxon

my $total_taxa_count=0;
my $in_list_count=scalar @acc;
foreach	my $tx(keys %{$tree->{tx_id}}){$total_taxa_count++;}
my $left_out_count=$total_taxa_count-$in_list_count;
my ($min_count)=sort {$a <=> $b} ($in_list_count,$left_out_count);
$min_count=int($min_count/2);
if($min_count<2){$min_count=2;}

my $acc_ref;
foreach $acc(@acc){$acc_ref->{$acc}=1;}

if($left_out_count<1){die "error:the input taxa include all the taxa in the tree\n";}
if($left_out_count==1){
foreach my $tx(keys %{$tree->{tx_id}}){
   if($acc_ref->{$tx}){next;}
  else{@acc=($tx);last;}}
}


if(@acc == 1){
 $acc=$acc[0];
 my $tID=$tree->{tx_id}->{$acc};
 my $nID=$tree->{parent}->{$tID};
 $tree->root_mid_between_node_tx($nID,$tID);
 open(OUT,">$opt_o") || die "cannot create output file $opt_o \n";
 print OUT $tree->output_newick();
 close OUT;

 open(OUT,">$logfile") || die "cannot output log file $logfile \n";
 print OUT "root mid point between $acc and its parent node \n";
 close OUT;
 exit;
}

#####split all the edges and evaluate the clusters till the best matches are met

my %cut;
my %tx;
my $acc_ref;
my @max_c1;
my @max_c2;
my $max;


foreach $acc(@acc){$acc_ref->{$acc}=1;}

foreach my $t(keys %{$tree->{id_tx}}){
   $tx{$t}=$tree->{id_tx}->{$t};
    my $this_tn=$t;
    while($tree->{parent}->{$this_tn}){
     if($tree->{branch_length}->{$this_tn} > 0){    ##edage need to be >=0 to be cut
           $cut{$this_tn}->{$t}=1;
     }
    $this_tn=$tree->{parent}->{$this_tn};
   }
}


foreach my $tn(keys %cut){
my @c1;
my @c2;
my %tmp;
   foreach my $c(keys %{$cut{$tn}}){
   push(@c1,$tx{$c});
   $tmp{$c}=1;
  }

  foreach my $c(keys %tx){
  if($tmp{$c}){next;}else{push(@c2,$tx{$c});}
  }

if((@c1 >= $min_count) && (@c2 >= $min_count)) {
  my ($p,$r)=&get_pr($acc_ref,\@c1);
  my $f1=&Fmeasure($p,$r);
  if($f1=~/\d/){if($f1>$max){$max=$f1;@max_c1=@c1; @max_c2=@c2;}}
  
 
  ($p,$r)=&get_pr($acc_ref,\@c2);
  $f1=&Fmeasure($p,$r);
  if($f1=~/\d/){if($f1>$max){$max=$f1;@max_c1=@c1; @max_c2=@c2;}}
  } 
  if($max==1){last;}
}

if($max < 0.5){
print STDERR "Cannot output to $opt_o , see $logfile for detail \n"; 
open(LOG,">$logfile") || die "cannot create log file $logfile \n";
print LOG "failed to find comparable clade to the input outgroup list: best F1 score is only $max \n";
close LOG;
exit;
}



###############root the tree between @c1 and @c2####

my %in;
$tree->root_to_txID($tree->{tx_id}->{$max_c1[0]});
foreach my $t(@max_c1){
 my $this=$tree->{tx_id}->{$t};
 $in{$this}=1;
while($tree->{parent}->{$this}){
   $in{$this}=1;
   $this=$tree->{parent}->{$this};
  }
 }

my $this=$tree->{tx_id}->{$max_c2[0]};
my $node_r;
my $node_l;

while($tree->{parent}->{$this}){
   my $p=$tree->{parent}->{$this};
   if($in{$p} && !($in{$this})){
     $node_l=$this;$node_r=$p;last;
    }
    $this=$p;
  }

if($node_l=~/^N/){
$tree->root_mid_between_nodes($node_r,$node_l);
 }
elsif($node_l=~/^T/){
$tree->root_mid_between_node_tx($node_r,$node_l);
}
##reroot

open(OUT,">$opt_o") || die"cannot create output file\n";
print OUT $tree->output_newick();
close OUT;

open(LOG,">$logfile") || die "cannot create log	file $logfile \n";
print LOG "####max F1: $max \n";
print LOG "###outgroup:\n".join("\n",@max_c1)."\n";
close LOG;


##############




sub Fmeasure{
my ($p,$r,$beta)=@_;
if(!$beta){$beta=1;}
$beta=$beta*$beta;

if($beta*$p+$r==0){return "NA";}

return sprintf("%.4f",(1+$beta)*$p*$r/($beta*$p+$r));

}


sub get_pr{
  my ($std_ref,$q_array)=@_;
  my $tp=0;
  my $fp=0;
  my $fn=0;

  my %otu;
  foreach my $t(@{$q_array}){
   $otu{$t}=1;
  }
  foreach my $t(keys %$std_ref){if($otu{$t}){$tp++;}else{$fp++;}}
  foreach my $t(keys %otu){if(!($std_ref->{$t})){$fn++;}}
  my $precision=0;
  my $recall=0;
  if(($tp+$fp)>0){$precision=sprintf("%.5f",$tp/($tp+$fp));}
  if(($tp+$fn)>0){$recall=sprintf("%.5f",$tp/($tp+$fn));}
  return ($precision,$recall);
}

#############################
package TREE_MOD;
use strict;

sub sub_tree_info{
###take a array of taxa, find the common acestor
##report: nodeID node_dist_from_root node_weighted_average_dist_from_input_taxa subtree_PD node_label
my ($class,@tx)=@_;
my $array_size=scalar @tx;
if($array_size ==0){return "##ERROR: empty_input";}
my @id;

foreach my $t(@tx){
if(length($class->{tx_id}->{$t})>=1){push(@id,$class->{tx_id}->{$t});}
else{return "##ERROR: $t is not in the tree";}
}

if($array_size ==1 ){
return $id[0]."\t".$class->{depth}->{$id[0]}."\t0\t0\tleaf";
}


my %pass_count;


foreach my $id(@id){
my $t=$id;
$pass_count{$t}++;
   while(length($class->{parent}->{$t})){
    $t=$class->{parent}->{$t};
    $pass_count{$t}++;
  }
}


my $t=$id[0];
my $common_node;
while(length($class->{parent}->{$t})){
my $p=$class->{parent}->{$t};
if($pass_count{$p}==$array_size){
$common_node=$p;last;
}
$t=$p;
}

my $common_depth_from_root=$class->float_number($class->{depth}->{$common_node});

my $label=$class->{label}->{$common_node};
my $sub_PD;
foreach my $t(keys %pass_count){
if($pass_count{$t} >= $array_size){
next;
 }
$sub_PD+=$class->{branch_length}->{$t};
}

###calculate weighted depth of the common node
###  @id, $common_node, $pass_count
my %true_depth_to_common;
my %inflated_depth_to_common;
my $weighted_node_to_tip=0;

foreach my $id(@id){
my $t=$id;
   while(length($class->{parent}->{$t})>=1){
   my $p=$class->{parent}->{$t};
   $inflated_depth_to_common{$id}+=$class->{branch_length}->{$t};
   $true_depth_to_common{$id}+=$class->{branch_length}->{$t}*(1/$pass_count{$t});
   if($p eq $common_node){last;}else{$t=$p;}
   }

  $weighted_node_to_tip+=$inflated_depth_to_common{$id}*$true_depth_to_common{$id};

}


  if($sub_PD == 0){
   $weighted_node_to_tip="NA";
  }
  else{
  $weighted_node_to_tip=$class->float_number($weighted_node_to_tip/$sub_PD);
  }

   $sub_PD=$class->float_number($sub_PD);
 
return $common_node."\t".$common_depth_from_root."\t".$weighted_node_to_tip."\t".$sub_PD."\t".$label;

}




##six_point_rooting: take in two sets of three taxa, identify the two defined nodes
##                   root the tree in the middle of the two nodes

sub six_point_rooting{
my ($class,@taxa)=@_;
my $count=scalar @taxa;
if(!($count==6)){return 0;}
my (@array1,@array2);
while($count>3){
my $tx=pop @taxa;
$count--;
push(@array1,$tx);
}
@array2=@taxa;


##root the tree to the first tx in @array1
##identify the node for array1
foreach my $t(@array1){
    $t=$class->{tx_id}->{$t};
   

    if(!(length($t)>=1)){die "taxa in input list are not in the tree\n";}
}

foreach my $t(@array2){
    $t=$class->{tx_id}->{$t};
   
    if(!(length($t)>=1)){die "taxa in input list are not in the tree\n";}
}

$class->root_to_txID($array1[0]);
my %path=();
my $thisID=$array1[1];
$path{$thisID}=1;
while($class->{parent}->{$thisID}){
$thisID=$class->{parent}->{$thisID};
$path{$thisID}=1;
}

$thisID=$array1[2];
while($class->{parent}->{$thisID}){
$thisID=$class->{parent}->{$thisID};
if($path{$thisID}){last;}
}
my $node1=$thisID;

##root the tree to the first tx in @array2
##identify the node for array2
$class->root_to_txID($array2[0]);
%path=();
$thisID=$array2[1];
$path{$thisID}=1;
while($class->{parent}->{$thisID}){
$thisID=$class->{parent}->{$thisID};
$path{$thisID}=1;
}

$thisID=$array2[2];
while($class->{parent}->{$thisID}){
$thisID=$class->{parent}->{$thisID};
if($path{$thisID}){last;}
}
my $node2=$thisID;

if ($node1 eq $node2){
$class->root_to_nodeID($node1);
return 1;
}

##$class->root_to_nodeID($node1,$node2); replaced by the following line
$class->root_mid_between_nodes($node1,$node2);

return 1;
} 



#####function added 5/7/2012 
sub relative_density{
    my $class=shift;
    &tree_clean($class);
    &all_depth($class);

    my %pass_count; 

    foreach my $txID(keys %{$class->{id_tx}}){
    my $this_id=$txID;
    $pass_count{$this_id}=1;
       while(length($class->{parent}->{$this_id})>=1){
       $this_id=$class->{parent}->{$this_id};
       $pass_count{$this_id}++;
       }    
    }

    my %PD_node_txID;
    my %inflated_node_txID;    

    foreach my $txID(keys %{$class->{id_tx}}){
    my $this_id=$txID;
    my $parent_id;
    my $inflated_PD=0;
    my $PD=0;
 
       while(length($class->{parent}->{$this_id})>=1){
       $parent_id=$class->{parent}->{$this_id};
       $PD+=$class->{branch_length}->{$this_id}/$pass_count{$this_id};       
       $inflated_PD+=$class->{branch_length}->{$this_id};
       
       $PD_node_txID{$parent_id}->{$txID}=$PD; 
       $inflated_node_txID{$parent_id}->{$txID}=$inflated_PD;

       $this_id=$parent_id;
       }
    }





   my %relative_density;

   foreach my $nodeID(keys %{$class->{depth}}){
       if($nodeID=~/^T/){next;}

      my $subPD=0;

     foreach my $txID(keys %{$PD_node_txID{$nodeID}}){
        $subPD+=$PD_node_txID{$nodeID}->{$txID};
     }

     if($subPD == 0){next;}

    my $weighted_average=0;

    foreach my $txID(keys %{$PD_node_txID{$nodeID}}){
        $weighted_average+=$inflated_node_txID{$nodeID}->{$txID} * $PD_node_txID{$nodeID}->{$txID};
     }

    $weighted_average/=$subPD;

       if($weighted_average == 0){next;}

##print $class->{depth}->{$nodeID}."\t".$weighted_average."\t".$subPD."\n";
 
       my $relative_density=(($class->{depth}->{$nodeID}+$weighted_average)*$subPD)/($weighted_average*$weighted_average);
       if($relative_density >=1000){$relative_density{$nodeID}=sprintf("%.2e",$relative_density);}
       else{$relative_density{$nodeID}=sprintf("%.2f",$relative_density);}

   } ##end of foreach


    $class->{relative_density}=\%relative_density;

} ##end of sub



############modified to the module 7/15/2012

sub OTU_collapse{ ##check 2012/03/24

##this function will work on the tree object directory, thus once called 
##the tree structure will be altered
    my ($class,$line)=@_;   
    if($class->{root}=~/^T/){die "cannot call OTU_collapse function with a tree rooted from one of the taxa\n";}

    my ($acc,@tx)=split(/\s+/,$line);
##check if the depth calculation is there, if yes, skip, if no do so
    if(!(length($class->{depth})>=1)){
	$class->all_depth();
    }
##calculate the going through time of all the nodes covered by the input tx
## if one tx is not in the tree, means there are repeats in the input OTU list, die in the circumstance
### call sub_tree_info function to get the common node, and weighted aeverage branch length  

  if(@tx==1){return;}
  my $subtree_info=$class->sub_tree_info(@tx);
  my ($OTU_node,$OTU_node_depth,$OTU_node_branch_length,$sub_PD,$label)=split(/\s+/,$subtree_info); 
  if (!($OTU_node=~/^N/)){return;}

    my $log;
    my $node_depth=$class->{depth}->{$OTU_node};
    my $node_depth_out=$class->format_number($node_depth);
    
     

    foreach my $tx(@tx){
        my $this_id=$class->{tx_id}->{$tx};
        my $t=$class->{depth}->{$this_id} - $node_depth;
	$log.=$acc."\t".$node_depth_out."\t".$tx."\t".$class->format_number($t)."\n";
    }
##record the depth of the collaping node and all the branch length from the tx to this node
 
	$log="#TAXON:".$acc."\tAttached_node_depth:".$node_depth_out."\tWeighted_average_Leave_Depth:".$OTU_node_branch_length."\n".$log;

##remove the input taxa
   foreach my $tx(@tx){
       my $txID=$class->{tx_id}->{$tx};
       delete $class->{id_tx}->{$txID};
       delete $class->{tx_id}->{$tx};
       delete $class->{parent}->{$txID};
       delete $class->{branch_length}->{$txID};
       delete $class->{depth}->{$txID};
   }

##add a leaf to the tree 
    my $OTU_number=0;
    if($class->{max_tx_number}>=1){
	$class->{max_tx_number}++;
        $OTU_number=$class->{max_tx_number};
    }
    else{
	$class->taxa_ID_count();
	$class->{max_tx_number}++;
        $OTU_number=$class->{max_tx_number};
    }

    my $txID="T".$OTU_number."DYW";
    
##    print $OTU_node."\t".$tx."\t".$ave_depth."\n";
    $class->{id_tx}->{$txID}=$acc;
    $class->{tx_id}->{$acc}=$txID;
    $class->{parent}->{$txID}=$OTU_node;
    $class->{branch_length}->{$txID}=$OTU_node_branch_length;
##collapsing the tx in the list, and output the stat
    return $log;
}


#############end tree collapes#######


sub taxa_ID_count{ ##check 2012/03/24
    my $class=shift;
    my  $max=0;
    foreach my $t(keys %{$class->{id_tx}}){
        if($t=~/(\d+)/){if($1>=$max){$max=$1;}}
    }
    $class->{max_tx_number}=$max;
    return $max;
}



sub OTU_by_normalized_depth_from_tip_weighted{ ##check 2012/03/23
    my ($class,$cutoff)=@_;
    if(!($cutoff>=0 && $cutoff<1)){
        die "the normalized cutoff from root must be larger than 0 (root) and <= 1 (break into individual taxa)\n";
    }
    $cutoff=1-$cutoff;
    my $cluster=&OTU_by_normalized_depth_from_root_weighted($class,$cutoff);
    return $cluster;
}

############
sub OTU_by_normalized_depth_from_root_weighted{
my ($class,$cutoff)=@_;

if($class->{root}=~/^T/){die "OTU by normalized depth need a tree that is NOT rooted at one of its taxa\n";}
    foreach my $t(keys %{$class->{depth}}){
        if($class->{depth}->{$t} == 0 && $t=~/^T/){
            die "there are leaves with 0 PD from the root, please reroot the tree \n";
        }
    }

   if(!($cutoff<=1 && $cutoff>0)){
        print $cutoff."  &&& \n";
        die "the normalized cutoff from root must be larger than 0 (root) and <= 1 (break into individual 
taxa)\n";
    }

    &tree_clean($class);
    &all_depth($class);

################

my (%child,@proc,$cluster);
my %leaf_to_node_PD_contribution;
my %leaf_to_node_depth;

    foreach my $t(keys %{$class->{parent}}){
        my $p=$class->{parent}->{$t};
        $child{$p}->{$t}=1;
        if($t=~/^T/){push (@proc,$t);
          $leaf_to_node_PD_contribution{$t}=$class->{branch_length}->{$t};
          $leaf_to_node_depth{$t}=$class->{branch_length}->{$t}; 
          } ##end of if 
    } ##end of foreach

while(@proc >=1){
my $proc;
while(@proc){
$proc=pop @proc;
###if $proc if a leaf
if($proc=~/^T/){
  my ($verdict,$depth);
  my $p=$class->{parent}->{$proc};
  my $depth=$class->{depth}->{$p}/$class->{depth}->{$proc};
  $depth=sprintf($class->{in_float},$depth);
 ##verdict 0 no cut verdict 1 cut
  if($depth >= $cutoff){$verdict=0;}else{$verdict=1;}
  if($verdict ==1){
  $cluster.=$class->{id_tx}->{$proc}."\n";
  delete $child{$p}->{$proc};
  }
###do nothing if the verdict is no cut
 } ##end if $proc is a leaf
elsif($proc=~/^N/){
   my $verdict;
   my $depth;
   my @child=();
   my $p=$class->{parent}->{$proc};
   foreach my $tmp(keys %{$child{$proc}}){push(@child,$tmp);}

   if(!(scalar @child) >=1){
        if($child{$p}->{$proc}>=1){delete $child{$p}->{$proc};}
        next;
       }
###############
             my  $child_ave=0;
             my  $subPD=0;  

      foreach my $child(@child){
      $subPD+=$leaf_to_node_PD_contribution{$child};                                          
      $child_ave+=$leaf_to_node_PD_contribution{$child}*$leaf_to_node_depth{$child}; 
              }        ###end_of_foreach
   if($subPD == 0){$child_ave=0;}else{$child_ave=$child_ave/$subPD;} 
##############
$depth=($class->{depth}->{$p})/($class->{depth}->{$proc}+$child_ave); 
$depth=sprintf($class->{in_float},$depth);


                if($depth >= $cutoff){$verdict=0;}
                else{$verdict=1;}




############
               if($verdict == 0){
                    my $PDadd=0;                                 
                    my $depthadd=$class->{branch_length}->{$proc};
                    my $child_count=scalar @child;                          
                       $PDadd=$depthadd/$child_count;

                  foreach my $tmp(@child){
                    $child{$p}->{$tmp}=1;   ##add the leaves to the parent node
                       $leaf_to_node_PD_contribution{$tmp}+=$PDadd;  
                       $leaf_to_node_depth{$tmp}+=$depthadd; 				
                      } ##end of foreacch

                    delete $child{$p}->{$proc};

                } ##end if verdict ==0
              else{   ##if the verdict is cut
                        my @t=();
                        foreach my $tmp(keys %{$child{$proc}}){
                            push(@t,$class->{id_tx}->{$tmp});
                        } ##end of foreach
                        $cluster.=join("\t",@t)."\n";
                    delete $child{$p}->{$proc};
              } ##end of else if verdict ==1

} ##end of elsif $proc is N
else{
die "ERROR non Node and leave in tree structure\n";
 } ##end of else

delete $child{$proc};

} ##end of while loop

##rebuild @proc 

################
##decide the next round of nodes

###identify the node that have only leaves as children
@proc=();
foreach my $t(keys %child){
    my $n_count=0;
    my $t_count=0;
    foreach my $child(keys %{$child{$t}}){
      
      if($child=~/^N/){$n_count++;}
      if($child=~/^T/){$t_count++;}
   } ##end $child
   if($t_count>=1 && $n_count==0 ){push(@proc,$t);}
   if($t_count==0 && $n_count==0) {delete $child{$t}; 
    my $this_p=$class->{parent}->{$t};
    delete $child{$this_p}->{$t};
   }
}  ##end $t

   

        } ### end of if(@proc>=1)

return $cluster;
 
    }

		                




sub average{ ##check 2012/03/24
    my @t=@_;
    if(!(@t>=1)){
	return 0;
    }
    my $total=0;
    my $count=scalar @t;
    my $t;

    foreach $t(@t){
	$total+=$t;
    }
   $total=$total/$count;
   return $total;   
}



sub all_depth{ ##check 2012/03/24

##this function calculate the PD from node to the root
##and store the results in $class->{node_depth}
    my $class=shift;
    undef $class->{depth};
    my %depth;
    $depth{$class->{root}}=0;    

    my ($tx,$this,@t,$depth,$pa);
    foreach $tx(keys %{$class->{id_tx}}){
	$this=$tx;
        $depth=0;
        $pa="";   
        @t=();
        while(length($this)>=1){
            push(@t,$this);
	    $pa=$class->{parent}->{$this};
            if(length($depth{$pa})>=1){last;}
            $this=$pa;
	}
        $depth=$depth{$pa};
        while(@t >=1){
	    $this=pop @t;
            $depth+=$class->{branch_length}->{$this};
            $depth{$this}=$depth;
       }	
    }

    $class->{depth}=\%depth;

    return $class;


}



sub subtree{ ##check 2012/03/24
    my ($new_sub,$class,$list_ref)=@_;
    my @list=@$list_ref;
    my ($tx,$txID,$this);
    my (%parent,%branch_length,%tx_id,%id_tx,%label,$root);
   
    foreach $tx(@list){
	if(length($class->{tx_id}->{$tx})>=1){
            $tx_id{$tx}=$class->{tx_id}->{$tx};
            $id_tx{$class->{tx_id}->{$tx}}=$tx;
	}
	else{
            die "$tx is not in the tree \n";
	}
    }

    foreach $txID(keys %id_tx){
	$this=$txID;
       
	while(length($class->{parent}->{$this})>=1){
	    $parent{$this}=$class->{parent}->{$this};
	    $branch_length{$this}=$class->{branch_length}->{$this};
            if(length($class->{label}->{$this})>=1){$label{$this}=$class->{label}->{$this};}
          
            $this=$class->{parent}->{$this};
           
	}
    }

############

    $root=$class->{root};
    if(length($class->{label}->{$root})>=1){$label{$root}=$class->{label}->{$root}; }



my $self={};

$self->{label}=\%label;
$self->{parent}=\%parent;
$self->{branch_length}=\%branch_length;
$self->{root}=$root;
$self->{tx_id}=\%tx_id;
$self->{id_tx}=\%id_tx;
$self->{in_float}=$class->{in_float};
$self->{out_float}=$class->{out_float};


$self=&tree_clean($self);

   bless $self, $new_sub;
   
   return $self;

}


sub tree_clean{  ##check 2012/03/23
##get rid of zero edges that connect two nodes
##get rid of passing though nodes
##input is tree objects
    my $class=shift;
    my ($t,$p,$c,%child,@zero_node);

##delete all the nodes that are not touched by any taxa
    my %in;
    foreach $t(keys %{$class->{id_tx}}){
	my $this=$t;
        $in{$this}=1;
        while(length($class->{parent}->{$this}) >=1){
	    $p=$class->{parent}->{$this};
            $in{$p}=1;
            $this=$p;
	}
    }

    foreach $t(keys %{$class->{parent}}){
        if($in{$t} ==1){next;}
	delete $class->{parent}->{$t};
	delete $class->{branch_length}->{$t};
    }



###go thought the data structure once, record nodes that are 0 dist from their parents
## keep all the children of a node
    foreach $t(keys %{$class->{parent}}){
	$p=$class->{parent}->{$t};
        $child{$p}->{$t}=1;
        if(($t=~/^N/) && ($class->{branch_length}->{$t} == 0)){
            push(@zero_node,$t);
	}
    }

    foreach $t(@zero_node){
	$p=$class->{parent}->{$t};
        foreach $c(keys %{$child{$t}}){
	    $class->{parent}->{$c}=$p;
            $child{$p}->{$c}=1;   
	}
    
	delete $child{$t};
        delete $child{$p}->{$t};
        delete $class->{parent}->{$t};
    }

###end of zero edges


   
###get rid of passing though nodes   

    my ($child_count,@pass);

    foreach $t(keys %child){
	$child_count=0;
        foreach $c(keys %{$child{$t}}){
           $child_count++;     
	}
	if(($child_count==1)&&(length($class->{parent}->{$t})>=1)){
	    push(@pass,$t);
	}
    }

    while(@pass >=1){
	$t=shift @pass;
        $p=$class->{parent}->{$t};
        foreach my $tmp(keys %{$child{$t}}){$c=$tmp;}
	    $class->{parent}->{$c}=$p;
            $class->{branch_length}->{$c}+=$class->{branch_length}->{$t};
            delete $child{$p}->{$t};	   
      	 $child{$p}->{$c}=1;  
          delete $class->{parent}->{$t};
          delete $class->{branch_length}->{$t};
    }


##deal with the stick out root 
##if the root is a taxon, leave it alone,if the root is a node, cut it out
    $child_count=0;
    my $root=$class->{root};
    if(($root=~/^N/)||( $root=~/^T/ && !(length($class->{id_tx}->{$root}) >=1) )){
     foreach my $tmp(keys %{$child{$root}}){$c=$tmp;$child_count++;}
     if($child_count == 1){
     ##    print $root."\t".$c."\t".$class->{branch_length}->{$c}."\n";
	 $class->{root}=$c;
         delete $class->{parent}->{$c};
         delete $class->{branch_length}->{$c};
     }
    }

    return $class;
}



sub mid_point_root{ ##check 2012/03/24
    my $class=shift;
    my ($to,$from,$max,$mid,$node_number,$id,$i,$t,$dist_from,$dist_to);
    my @path;
 
   ############
    ($from,$max)=$class->max_path_by_reroot();

   ############


    $i=0;
    foreach $t(keys %{$class->{parent}}){
	if($t=~/N(\d+)DYW/){
	    if($1 > $i){$i=$1;}
        }
    }
    $i++;
    $id="N".$i."DYW";  ## this is the Node of midpoint
##  create a mid point between tip and root
    $mid=$class->float_number($max/2);
    
    $dist_from=0;
    $dist_to=0;

    while(length($class->{parent}->{$from})>=1){
	$to=$class->{parent}->{$from};
        $dist_to+=$class->{branch_length}->{$from};
        if(($dist_to >= $mid) && ($dist_from < $mid)){last;}     
        $from=$to;
        $dist_from=$dist_to;
    }
    

## 
    $dist_to=$class->float_number($dist_to-$mid);
    $dist_from=$class->float_number($mid-$dist_from);

    if($dist_to == 0){
	$class->root_to_nodeID($to);
    }
    else{
##create mid point ID 
	$class->{parent}->{$from}=$id;
        $class->{branch_length}->{$from}=$dist_from;
        
	$class->{parent}->{$id}=$to;
        $class->{branch_length}->{$id}=$dist_to;
	$class->root_to_nodeID($id);
        $class->{label}->{$id}="";
    }
  
    return 1;


}




sub format_number{ ##check 2012/03/23
    my ($class,$number)=@_;
    $number=sprintf($class->{out_float},$number);
    $number=~s/([^0])0+$/$1/;
    $number=~s/\.$//;
    return $number;
}


sub float_number{ ##check 2012/03/23
    my ($class,$number)=@_;
    return sprintf($class->{in_float},$number);
}


sub max_path_by_reroot{ ##check 2012/03/24  
##this step include reroot, the root will be distroyed
    my $class=shift;
    my $start_id=$class->pick_one_txID();
    $class->root_to_txID($start_id);
    my ($tx1,$tx2,$max);
    ($tx1,$max)=$class->deepest_txID();
    $class->root_to_txID($tx1);
    ($tx2,$max)=$class->deepest_txID();
    return ($tx2,$max);
}



sub deepest_txID{ ##check 2012/03/24
    my $class=shift;
    my $farest_tx;
    my $max=-999999;
    my $tx;

    $class->all_depth();

    foreach $tx(keys %{$class->{id_tx}}){ ##modified 05/01/2012
	if($class->{depth}->{$tx} >=$max){
	    $max=$class->{depth}->{$tx};
            $farest_tx=$tx;
	}
    }

    return ($farest_tx,$max);
}


sub output_newick{ ## check 2012/03/24
    my $class=shift;
    my (%in,%remain,%child);
    my ($id,$upper,$t,$tmp);
    my (@proc,@t,%tmp);
    my $skip=1;
    my $float=$class->{float};
    my $newick; 

    my $mod_node=""; ## if the tree is rooted at a taxon
    my $mod_tx="";

 
    if($class->{root}=~/^T/){
        $mod_tx=$class->{root};
	$class->{parent}->{$mod_tx}='FAKEROOT';
	$class->{branch_length}->{$mod_tx}=0;

        my $tmp_tx=$class->pick_one_txID();
        while(length($class->{parent}->{$tmp_tx})>=1){
            if($class->{parent}->{$tmp_tx} eq $class->{root}){$mod_node=$tmp_tx; last;}
	    $tmp_tx=$class->{parent}->{$tmp_tx};
	}
       ###tmp_tx is the one node connect to root

	if(!($mod_node=~/^N/)){die "cannot output newick, tree struct problem\n";}
        $class->{parent}->{$mod_node}='FAKEROOT';
        $class->{root}='FAKEROOT'; 
} 
    
 ##set a hash for all the children

####if the root is tx, make a fake node ###########


    foreach $t(keys %{$class->{parent}}){
	$upper=$class->{parent}->{$t};
        $remain{$upper}=1;
        push (@{$child{$upper}},$t);
    }

    foreach $id(keys %{$class->{id_tx}}){
        $in{$id}=$class->{id_tx}->{$id};
        $tmp{$class->{parent}->{$id}}=1;
    }
################

    foreach $t(keys %tmp){
	my $good=0;
        my $bad=0;
        foreach $tmp(@{$child{$t}}){
	    if(length($in{$tmp})>=1){$good++;}
            if(length($remain{$tmp})>=1){$bad++;}
	}
	if($bad>=1){next;}
        if($good>=1){push(@proc,$t);}
    }

while(@proc>=1){
##merge all the nodes in @proc

    foreach $id(@proc){
	@t=();
    foreach $t(@{$child{$id}}){
	push(@t, $in{$t}.":".$class->format_number($class->{branch_length}->{$t}));
        delete $in{$t};
    }
	$in{$id}="(".join(",",@t).")";
 
        if(length($class->{label}->{$id})>=1){
	    $in{$id}.=$class->{label}->{$id};
	}
    delete $remain{$id}; 
    }  ###end of proc id

##figure out the next bunch
    @proc=();
    %tmp={};
    foreach $t(keys %in){$tmp{$class->{parent}->{$t}}=1;}
  
    foreach $t(keys %tmp){
	my $good=0;
        my $bad=0;
        foreach $tmp(@{$child{$t}}){
	    if(length($in{$tmp})>=1){$good++;}
            if(length($remain{$tmp})>=1){$bad++;}
	}
	if($bad>=1){next;}
        if($good>=1){push(@proc,$t);}
    }
}    



    $newick=$in{$class->{root}}.";";
 
    if($class->{root} eq 'FAKEROOT'){
	$class->{parent}->{$mod_node}=$mod_tx;
        $class->{root}=$mod_tx;
        delete 	$class->{parent}->{$mod_tx};
        delete $class->{branch_length}->{$mod_tx};
    }
   

   
    return $newick;

}


sub root_mid_between_node_tx{
my ($class,$node,$tx)=@_;
if($class->{parent}->{$tx} ne $node){return 0;}


my   $i=0;
    foreach my $t(keys %{$class->{parent}}){
        if($t=~/N(\d+)DYW/){
            if($1 > $i){$i=$1;}
        }
    }
    $i++;
    my $newID="N".$i."DYW";  ## this is the new node for rooting

###############
$class->root_to_nodeID($node);
my $mid=$class->float_number(($class->{branch_length}->{$tx})/2);


$class->{parent}->{$tx}=$newID;
$class->{branch_length}->{$tx}=$mid;
$class->{parent}->{$newID}=$node;
$class->{branch_length}->{$newID}=$mid;

$class->root_to_nodeID($newID);
return 1;
}

sub root_mid_between_nodes{
my ($class,$node1,$node2)=@_;


if ($node1 eq $node2){
$class->root_to_nodeID($node1);
return 1;
}


my   $i=0;
    foreach my $t(keys %{$class->{parent}}){
        if($t=~/N(\d+)DYW/){
            if($1 > $i){$i=$1;}
        }
    }
    $i++;
    my $newID="N".$i."DYW";  ## this is the new node for rooting

##insert $newID between $node1 and $node2
##root to node1
$class->root_to_nodeID($node1);
##identify the distance between node2 to root
my $dist=0;
my $thisID=$node2;
while($class->{parent}->{$thisID}){
$dist+=$class->{branch_length}->{$thisID};
$thisID=$class->{parent}->{$thisID};
}
my   $mid=$class->float_number($dist/2);

my $dist_from=0;
my $from=$node2;
my $to;
my $dist_to=0;


while($class->{parent}->{$from}){
$to=$class->{parent}->{$from};
$dist_to+=$class->{branch_length}->{$from};
if($dist_to>=$mid){
last;
}
$from=$to;
$dist_from=$dist_to;
}

if($dist_to > $mid){
$class->{parent}->{$from}=$newID;
$class->{branch_length}->{$from}=$mid-$dist_from;
$class->{parent}->{$newID}=$to;
$class->{branch_length}->{$newID}=$dist_to-$mid;
}

$class->root_to_nodeID($newID);
return 1;

}


sub root_to_nodeID{ ##check 2012/03/24
    my ($class,$nodeID)=@_;
    my $old_root=$class->{root};
    if($nodeID eq $old_root){return;}

    my $thisID=$nodeID;
    my @path=($nodeID);

    my (%tmp_parent,%tmp_length);
    my $id;
    my $upper;

    while(length($class->{parent}->{$thisID})>=1){
	$thisID=$class->{parent}->{$thisID};
        push(@path,$thisID);
   } 


    while(@path>=1){
	$id=shift(@path);
    if(@path>=1){
        $upper=$path[0];
        $tmp_parent{$upper}=$id;
        $tmp_length{$upper}=$class->{branch_length}->{$id};          
     }    
    }

    foreach $id(keys %tmp_parent){
      $class->{parent}->{$id}=$tmp_parent{$id};
      $class->{branch_length}->{$id}=$tmp_length{$id};
    }

    delete $class->{parent}->{$nodeID};
    delete $class->{branch_length}->{$nodeID};

  
    my @down;
    foreach $id (keys %{$class->{parent}}){
    if($class->{parent}->{$id} eq $old_root){
	push(@down,$id);
    }
    }

   if(@down == 1){
       my $up=$old_root;
       my $down=$down[0];
 
       $class->{parent}->{$down}=$class->{parent}->{$up};
       $class->{branch_length}->{$down}=$class->float_number($class->{branch_length}->{$down}+$class->{branch_length}->{$old_root});
       delete  $class->{parent}->{$old_root};
       delete  $class->{branch_length}->{$old_root};
   } 

    $class->{root}=$nodeID;   


}


sub root_to_txID{ ##check 2012/03/24
    my ($class,$txID)=@_;
    my $old_root=$class->{root};
       
    if($old_root eq $txID){return;} 
    my $thisID=$txID;
    my @path=($txID);
    my (%tmp_parent,%tmp_length);
    my $id;
    my $upper;

    while(length($class->{parent}->{$thisID})>=1){
	$thisID=$class->{parent}->{$thisID};
        push(@path,$thisID);
   } 

#####

    while(@path>=1){
	$id=shift(@path);
    if(@path>=1){
        $upper=$path[0];
        $tmp_parent{$upper}=$id;
        $tmp_length{$upper}=$class->{branch_length}->{$id};          
     }    
    }

    foreach $id(keys %tmp_parent){
      $class->{parent}->{$id}=$tmp_parent{$id};
      $class->{branch_length}->{$id}=$tmp_length{$id};
    }

    delete $class->{parent}->{$txID};
    delete $class->{branch_length}->{$txID};

  

##if original root node connect two nodes, delete it

    my @down;
    foreach $id (keys %{$class->{parent}}){
    if($class->{parent}->{$id} eq $old_root){
	push(@down,$id);
    }
    }

   if(@down == 1){
       my $up=$old_root;
       my $down=$down[0]; 
       $class->{parent}->{$down}=$class->{parent}->{$up};
       $class->{branch_length}->{$down}=$class->float_number($class->{branch_length}->{$down}+$class->{branch_length}->{$old_root});
       delete  $class->{parent}->{$old_root};
       delete  $class->{branch_length}->{$old_root};
   } 
    $class->{root}=$txID;   
   
}

sub pick_one_txID{ ##check 2012/03/24
    my $class=shift;
    my $id;
    my $return_id;
    foreach $id(keys %{$class->{id_tx}}){
	if($class->{branch_length}->{$id} > 0){
        $return_id=$id;
        last; 
        }
    }
    return $return_id;
}

sub str4pattern{  ##check 2012/03/24
##this function take in a tring and make it ofr patern match
##it backslash \ | ( ) [ { ^ $ * + ? .
    my $str=shift;
    $str=~s/([\\\|\(\)\[\{\^\$\*\+\?\.])/\\$1/g;
    return $str;
}





sub new{  ##check 2012/03/23
my ($class,$tree)=@_;
my @node;
my $node_i=0;
my $tx_i=0;
my (%parent,%label,%edge_length);
my (%tx_id,%id_tx);
my ($tx,$value,$label, $tmp_str,$track,$root);
my ($txID,$nodeID,$i,@t,$len);
my $float=0;


@t=split(//,$tree);
$len=scalar @t;

for $i(0..$len-1){
if($t[$i] eq "\n"){next;}
if($t[$i] eq "\t"){next;}

if($t[$i] eq '('){
##start a node: create a nodeID,reset $tx,update track
$track=$t[$i];
$node_i++;
$nodeID="N".$node_i."DYW";
push(@node,$nodeID);
$tmp_str="";
next;}

if($t[$i] eq ':'){
##start of a value for node or taxon, 
##if it is for a node
if($track eq ')'){
   $label=$tmp_str;  
  }
else{
##if it is for a tx
  $tx=$tmp_str; 
  }
$tmp_str='';
next;
}

if($t[$i] eq ','){
## a internal tx or node
$value=$tmp_str;

#####setup float point
$value=~s/^\s+//;
$value=~s/\s+$//;
if($value=~/\.(\d+)$/){if(length($1) > $float){$float=length($1);}}
####end of float point setting

$tmp_str='';

if($track eq ")"){
##setup the internal node
$nodeID=pop @node;
$parent{$nodeID}=$node[-1];
$edge_length{$nodeID}=$value;
if(length($label)>=1){
$label{$nodeID}=$label;
$label='';
  }
$track=',';
next;
}
else{
###end of internal tx, track should be either , or (
##update relation branch and relationship 
##reset tx and value
$tx_i++;
$txID="T".$tx_i."DYW";

$parent{$txID}=$node[-1];
$edge_length{$txID}=$value;
$tx_id{$tx}=$txID;
$id_tx{$txID}=$tx;
$tx='';
$value='';
$track=',';
next;
}
}

if($t[$i] eq ')'){
###end of a node
$value=$tmp_str;
$tmp_str='';


#####setup float point
$value=~s/^\s+//;
$value=~s/\s+$//;
if($value=~/\.(\d+)$/){if(length($1) > $float){$float=length($1);}}
####end of float point setting


###if the end is a node
if($track eq ")"){
##setup the internal node
$nodeID=pop @node;
$parent{$nodeID}=$node[-1];
$edge_length{$nodeID}=$value;
if(length($label)>=1){
$label{$nodeID}=$label;
$label='';
   }
  }
else{
###if the end is a tx
###end of track should be either , or (
##update relation branch and relationship 
##reset tx and value
$tx_i++;
$txID="T".$tx_i."DYW";

$parent{$txID}=$node[-1];
$edge_length{$txID}=$value;
$tx_id{$tx}=$txID;
$id_tx{$txID}=$tx;
$tx='';
$value='';
}

##start of label and update track
$label='';
$tmp_str='';
$track=')';
next;
}

if(($t[$i] eq ';') && ($track eq ')')){
  $nodeID=pop(@node);
  if(length($tmp_str)>=1){$label{$nodeID}=$tmp_str;}
  $root=$nodeID;
  last;
}

  $tmp_str.=$t[$i];

}



if(@node == 1 && ($track eq ')')){
  $nodeID=pop(@node);
  if(length($tmp_str) >=1 ){$label{$nodeID}=$tmp_str;}
  $root=$nodeID;
}

if(@node >=1){
    die "tree format problem\n";
}


###turn all the edge into a float point number 
###setup internal float point
##setup external float point


my $self={};

$self->{label}=\%label;
$self->{parent}=\%parent;
$self->{root}=$root;
$self->{tx_id}=\%tx_id;
$self->{id_tx}=\%id_tx;
$float+=2;

$self->{out_float}='%.'.$float.'f';
$float+=2;
$self->{in_float}='%.'.$float.'f';

foreach my $t(keys %edge_length){
    $edge_length{$t}=sprintf($self->{in_float},$edge_length{$t});
} 
$self->{branch_length}=\%edge_length;

bless $self, $class;

$self->tree_clean();

   return $self;
}

sub output_struct{
    my $class=shift;
    my $struct;
##$self->{label}=\%label;
    $struct.="#label\n";
    foreach my $t(keys %{$class->{label}}){
	$struct.=$t."\t".$class->{label}->{$t}."\n";
    }
    $struct.="#end\n";
##$self->{parent}=\%parent;
    $struct.="#parent\n";
    foreach my $t(keys %{$class->{parent}}){
	$struct.=$t."\t".$class->{parent}->{$t}."\n";
    }
    $struct.="#end\n";
##$self->{root}=$root;
    $struct.="#root\n";
    $struct.="root\t".$class->{root}."\n";
    $struct.="#end\n";
##$self->{tx_id}=\%tx_id;
    $struct.="#tx_id\n";
    foreach my $t(keys %{$class->{tx_id}}){
	$struct.=$t."\t".$class->{tx_id}->{$t}."\n";
    }
    $struct.="#end\n";
##$self->{id_tx}=\%id_tx;
    $struct.="#id_tx\n";
    foreach my $t(keys %{$class->{id_tx}}){
	$struct.=$t."\t".$class->{id_tx}->{$t}."\n";
    }
    $struct.="#end\n";
##$self->{out_float}='%.'.$float.'f';
    $struct.="#out_float\n";
    $struct.="out_float\t".$class->{out_float}."\n";
    $struct.="#end\n";
##$self->{in_float}='%.'.$float.'f';
    $struct.="#in_float\n";
    $struct.="in_float\t".$class->{in_float}."\n";
    $struct.="#end\n";
##$self->{branch_length}=\%edge_length;
    $struct.="#branch_length\n";
    foreach my $t(keys %{$class->{branch_length}}){
	$struct.=$t."\t".$class->{branch_length}->{$t}."\n";
    }
    $struct.="#end\n";
return $struct;  
}

sub new_from_struct{
    my ($class,$struct)=@_;
    my $self={};
    my @t;
    my @line=split(/\n/,$struct);
    my $line;
    my $current="";
    while(@line){
	$line=shift @line;
        if($line=~/^#/){
	    $line=~s/^\#+//;
            $line=~s/\s+$//;
            if($line eq 'end'){$current="";}
            elsif(length($line)>=1){$current=$line;}
	    next;
	}
        @t=split(/\s+/,$line);
       if(!(@t>=2)){next;} 

	if($current eq 'label'){
	    $self->{label}->{$t[0]}=$t[1];
	}
        elsif($current eq 'parent'){
   	    $self->{parent}->{$t[0]}=$t[1];         
	}
        elsif($current eq 'root'){
	    $self->{root}=$t[1];
	}
        elsif($current eq 'tx_id'){
            $self->{tx_id}->{$t[0]}=$t[1];
	}
        elsif($current eq 'id_tx'){
            $self->{id_tx}->{$t[0]}=$t[1];
	}
        elsif($current eq 'out_float'){
	    $self->{out_float}=$t[1];
	}
        elsif($current eq 'in_float'){
	    $self->{in_float}=$t[1];
	}
        elsif($current eq 'branch_length'){
   	    $self->{branch_length}->{$t[0]}=$t[1];         
	}
    }

bless $self, $class;
return $self;

}


sub proc_seq{
## this one lineup node sequences for processing
## the results is a string of nodes under the key "proc_seq"
## also under key child the hash key=string of child will be setup
    my $class=shift;
##get all the children
    my (%tmp,%child);
    foreach my $c(keys %{$class->{parent}}){
	my $p=$class->{parent}->{$c};
        $tmp{$p}->{$c}=1;
    }
    foreach my $p(keys %tmp){
    my @c;
       foreach my $c(keys %{$tmp{$p}}){
	   push(@c,$c);
      } 
    @c=sort @c;
    $child{$p}=join("\t",@c);
    }

##set current children
    my %current;
    my @proc;
    %tmp=();
    foreach my $id(keys %{$class->{id_tx}}){
     $current{$id}=1;
     $tmp{$class->{parent}->{$id}}=1;   
    }

##get their parent
    
    foreach my $p(keys %tmp){
    my @c=split(/\t/,$child{$p});
    my $ready=1;
    foreach my $c(@c){
	if($tmp{$c}>=1){
	    next;
	}
        else{
	    $ready=0;
	}
    }
    if($ready>=1){push(@proc,$p);}
    } 

    my @final=();
    
    while(@proc){
##sort @proc
	@proc=sort @proc;
##output @proc;
        @final=(@final,@proc);
##register @proc;
##clean @proc and all its children
        foreach my $p(@proc){
        $tmp{$p}=1;
        my @c=split(/\t/,$child{$p});
        foreach my $c(@c){
	    delete $tmp{$c};
          }         
	}
	@proc=();
##set up next round
   
    foreach my $p(keys %tmp){
    my @c=split(/\t/,$child{$p});
    my $ready=1;
    foreach my $c(@c){
	if($tmp{$c}>=1){
	    next;
	}
        else{
	    $ready=0;
	}
    }
    if($ready>=1){push(@proc,$p);}
    }
    } ###loop until no @proc can be processed

    $class->{child}=\%child;
    $class->{proc_seq}=join("\t",@proc);    

    return 1;
}


1;

