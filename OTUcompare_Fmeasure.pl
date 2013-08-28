#!/usr/bin/perl

use strict;

#Usage: OTUcompare_Fmeasure.pl -query query_OTU -target target_OTUsets -shuffle 10 -beta 1
### -query query_OTU: the OTU group of interest, format: A,B,C,D
### -target target_OTUsets: clusters of OTU sets
### -shuffle integer: target random shuffle times to evaluate the expected maximum F measure score
### -beta: F measure beta factor, default 1
 


my %opt=@ARGV;

my $master=$opt{-query}; ##input OTU, format: A,B,C,D
my $clusterfile=$opt{-target}; ##clusterID  A,B,C  D  E,F 

my $random=10;

if(length($opt{-shuffle})>=1){$random=int($opt{-shuffle});}

 ## randomize cluster times, default 10, if target has only one OTU, this value is set to 0



#############
my $prune=$opt{-prune};         ## 1 means yes, 0 means no, default yes
                                ## for one cluster in put, $prune has to be set to 0
                                ## because it exclude members from the query that are not
                                ## in target clusters, thus one group cluster always end up with precision of 1
                                ## currently, if the target have only one OTU, $prune will be set to 0 automatically
 
my $beta=$opt{-beta}; ##beta for F measure, default 1 (0.5, emphasis precision, 2 emphasis recall)




my $random_hash;

##if($random>=1){$random_hash=&setup_rand_cluster_hash($clusterfile,$random);}

if(!$beta){$beta=1;}
##output the best F1 measure value for one cluster 
my %master_cluster;
my $master_acc_count=0;

if($prune=~/^([Nn]|0)/){$prune=0;}else{$prune=1;}

open(I,$master) || die;
while(<I>){
##my ($acc,$id)=split(/\s+/);
my @acc=split(/\s+|,/,$_);
foreach my $acc(@acc){if(length($acc)>=1){$master_cluster{$acc}=1;}}
}
close I;
foreach my $acc(keys %master_cluster){
$master_acc_count++;
}

open(I,$clusterfile) || die;
while(<I>){
if($_=~/\#/){next;}
if($_=~/^\s+$/){next;}
my $line=$_;
my @test=split(/\s+/,$line);
my $this_prune=$prune;
my $this_random=$random;

if(@test<=2){$this_prune=0;$this_random=0;} ##if a OTU set contains only one OTU, prune is turned off
my ($id,$maxF,$sub_master_acc_count)=&proc_cluster($line,\%master_cluster,$this_prune);
my $rand_ave=0;


if($this_random>1) {
my $ave=0;
for my $this_rand(1..$random){

##  my $rand_line=&get_rand_line_from_hash($line,$random_hash,$this_rand);

  my $rand_line=&get_rand_line($line);
  my ($id,$randF)=&proc_cluster($rand_line,\%master_cluster,$this_prune);
  $ave+=$randF;
 }

 $ave=sprintf("%.4f",$ave/$random);
$rand_ave=$ave;
}



my $adjust_maxF="NA";

if(1-$rand_ave > 0){$adjust_maxF=sprintf("%.4f",($maxF-$rand_ave)/(1-$rand_ave));
if($adjust_maxF <= 0){$adjust_maxF=0;}}
elsif($rand_ave == 1){$adjust_maxF=1;} ##a singleton against OTU set of singltons

if($this_prune==0){$sub_master_acc_count=$master_acc_count;}

print "targetID=".$id."\tAdjusted_max_Fscore=".$adjust_maxF;
print "\tMax_Fscore=".$maxF."\t"."Query_taxa_count=".$master_acc_count."\tQuery_taxa_evaluated=".$sub_master_acc_count."\tBETA=".$beta;
print "\tShuffle=".$this_random."\tExpected_maxF=".$rand_ave."\n";

}
close I;


sub get_rand_line{
my $line=shift;
my $rand_line;
my ($OTUid,@acc)=split(/\s+|,/,$line);
my ($rest,@t)=split(/\s+/,$line);

my $i=scalar @acc;
  while($i--){
     my $j=int rand($i+1);
     @acc[$i,$j]=@acc[$j,$i];
   } ##end while

my @rand_t;
foreach my $t(@t){
   my @a=split(/,/,$t);
   my $t_count=scalar @a;
   my $i=0;
   my @tmp;
   for $i(1..$t_count){my $tmp=pop @acc; push(@tmp,$tmp);}
   push(@rand_t,join(",",@tmp));
 }
 return "RAND".$OTUid."\t".join("\t",@rand_t)."\n";
}


sub proc_cluster{
my ($line,$master_cluster,$prune)=@_;

my ($id,@t)=split(/\s+/,$line);
my ($rest,@acc)=split(/\s+|,/,$line);

  ##make sure all the members in the master cluster are present in the OTUs

  my %sub_master_cluster;
  my %acc_this_group;
  my $master_acc_count=0;
  my $sub_master_acc_count=0;

  foreach my $acc(@acc){
    if($acc){
    if($master_cluster->{$acc}){$sub_master_cluster{$acc}=1; $sub_master_acc_count++;}
    $master_acc_count++;
    } 
  }


  my $maxF=0.0000;
  my $return_acc_count=0;



  foreach my $otu(@t){
  my ($precision,$recall);
   if($prune){($precision,$recall)=&get_pr(\%sub_master_cluster,$otu);$return_acc_count=$sub_master_acc_count;}
   else{($precision,$recall)=&get_pr(\%master_cluster,$otu);$return_acc_count=$master_acc_count; }
     my $f=&Fmeasure($precision,$recall,$beta);
     if($f=~/\d+/ && $f>$maxF){$maxF=$f;}
  }


return ($id,$maxF,$return_acc_count);
 }



sub Fmeasure{
my ($p,$r,$beta)=@_;
if(!$beta){$beta=1;}
$beta=$beta*$beta;

if($beta*$p+$r==0){return "NA";}

return sprintf("%.4f",(1+$beta)*$p*$r/($beta*$p+$r));

}



sub get_pr{
  my ($std_ref,$otu)=@_;
  my $tp=0;
  my $fp=0;
  my $fn=0;

  my %otu;
  my @t=split(/,/,$otu);
  foreach my $t(@t){
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
