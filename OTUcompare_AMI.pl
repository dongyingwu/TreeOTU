#!/usr/bin/perl

use strict;

##Usage: OTUcompare_AMI.pl -query query_OTUset -target target_OTUset -shuffle 100

my %opt=@ARGV;
my $master=$opt{-query}; ##standard cluster: format: ClusterID  A,B C,D (only one line is taken)
my $clusterfile=$opt{-target}; ##clusterID  A,B,C  D  E,F (can be multiple OTU sets, one per line)
my $rand_sample_time=$opt{-shuffle}; #default 100, if you don't want randomization input 0



my $report_cutoff=$opt{-cutoff}; ##only report ami >= this cutoff, default -9999 (everything)
if(!length($report_cutoff)){$report_cutoff=-9999;}
if(!(length($rand_sample_time)>=1)){$rand_sample_time=100;}

my $master_cluster;
my $qID;

open(I,$master) || die "cannot open query file\n";
while(<I>){
if(/^\s+$/){next;}
else{
my @t;
($qID,@t)=split(/\s+/);
my $q_i=0;

  foreach my $t(@t){
  my @acc=split(/,/,$t);
  $q_i++;
  my $qotu_id="MASTERQDWU".$q_i."ID";
   foreach my $acc(@acc){
    $master_cluster->{$acc}=$qotu_id;
   }
  }
last;
}

##my ($acc,$id)=split(/\s+/);
##if(length($acc) && length($id)){$master_cluster->{$acc}="MASTER".$id;}
}
close I;

open(I,$clusterfile) || die "cannot open target file\n";
while(<I>){
my ($id,@t)=split(/\s+/);
my $i=0;
if(@t){
my $this_cluster;
foreach my $t(@t){
$i++;
 my @acc=split(/,/,$t);
 foreach my $acc(@acc){$this_cluster->{$acc}="TMP".$i;}
 }
##now I have $this_cluster and $master_cluster

my ($clean1,$clean2,$clean_log)=&double_check_two_clusters($master_cluster,$this_cluster);
my $h1=&entropy($clean1);
my $h2=&entropy($clean2);
my $cluster3=&cross_cluster($clean1,$clean2);
my $h3=&entropy($cluster3);

my $I=$h1+$h2-$h3;
my $ami;

my $E;

if($rand_sample_time >=1){
$E=&get_E($clean1,$clean2,$rand_sample_time);
($E)=split(/\s+/,$E);
}
else{
$E=0;
}


if($h1+$h2-2*$E==0){$ami="NA";}
else{
$ami=sprintf("%.4f",(2*$I-2*$E)/($h1+$h2-2*$E));
}


if($ami >=$report_cutoff){
print "AMI=".$ami."\t";
print "Query_ID=".$qID."\t";
print "Target_ID=".$id."\t";
print "Shuffle=".$rand_sample_time."\t";

my ($q_count,$t_count,$c_count)=split(/,|:/,$clean_log);
print "Query_taxa_count=".$q_count."\t";
print "Target_taxa_count=".$t_count."\t";
print "Shared_taxa_count=".$c_count."\n";
}

##############
undef $this_cluster;
}
else{next;}

}
close I;
###########



sub get_E{
my ($cluster1,$cluster2,$r)=@_;
my $this_c1;
my $this_c2;
if(!$r){$r=100;}


foreach my $t(keys %$cluster1){
 $this_c1->{$t}=$cluster1->{$t};
 }

foreach	my $t(keys %$cluster2){
 $this_c2->{$t}=$cluster2->{$t};
 }

my $i;
my @I;
my $sum;
for $i(1..$r){
$this_c1=&shuffle_hash($this_c1);
$this_c2=&shuffle_hash($this_c2);
my $cross=&cross_cluster($this_c1,$this_c2);  
my $h1=&entropy($this_c1);
my $h2=&entropy($this_c2);
my $h3=&entropy($cross);
my $I=$h1+$h2-$h3;
push(@I,$I);
$sum+=$I;
undef $cross;
}

my $ave=$sum/$r;

return $ave."\t"."NA"."\t".$r;
#### I am not using standard deviation

my $std;
foreach my $I(@I){
$std+=($I-$ave)*($I-$ave);
}
if($r>1){$std=$std/($r-1);}
$std=sqrt($std);

return $ave."\t".$std."\t".$r;

}

sub shuffle_hash{
my $hash_ref=shift;
my @k;
my @v;
foreach my $t(keys %$hash_ref){
push(@k,$t);
push(@v,$hash_ref->{$t});
}

##shuffle the @v
my $i=scalar @v;
 while($i--){
     my $j=int rand($i+1);
     @v[$i,$j]=@v[$j,$i];
  }
while(@k){
my $k=pop @k;
my $v=pop @v;
$hash_ref->{$k}=$v;
}
return $hash_ref;
}

sub entropy{
my $c=shift;
my %hash;
my $total=0;
my $t;
my $h;
foreach $t(keys %$c){
  my $cID=$c->{$t};
  $hash{$cID}++;
 $total++;
 }

##foreach $t(keys %hash){$h+=($hash{$t}/$total)*(&log2($total)-&log2($hash{$t}));}
foreach $t(keys %hash){$h+=($hash{$t}/$total)*(log($total)-log($hash{$t}));} 

##for Mutual information bit and nat makes no difference

return sprintf("%.10f",$h);

}

sub log2{
        my $n = shift;
        return log($n)/log(2);
    }


sub double_check_two_clusters{
my ($c1,$c2)=@_;
my (%clean1,%clean2,%check);
my ($c1_count,$c2_count,$share_count);
foreach my $t(keys %$c1){$check{$t}="1"; $c1_count++;}
foreach my $t(keys %$c2){$check{$t}.="2"; $c2_count++;}

foreach my $t(keys %$c1){if($check{$t} eq "12"){$clean1{$t}=$c1->{$t}; $clean2{$t}=$c2->{$t};$share_count++;}}




return (\%clean1,\%clean2, $c1_count.",".$c2_count.":".$share_count);
}


sub cross_cluster{
my ($c1,$c2)=@_;
my %hash;
foreach my $t(keys %$c1){
$hash{$t}="FIRST".$c1->{$t};
}
foreach my $t(keys %$c2){
$hash{$t}.=":SECOND".$c2->{$t};
}
return \%hash;
}


sub get_cluster{
my ($n,$s,$track)=@_;
if($n<$s){die "n must be larger than s\n"}

my @cluster;
my %hash;
my $start=0;
for my $i(1..$s){
 my $id="C".$track."_".$i;
 push(@cluster,$id);
  $hash{$i}=$id;
  $start=$i+1;
 }

##every cluster must have one member (assign 1 to $s cluster 1 to $s)

if($start<$n && $start >=1){
for my $i($start..$n){
  my $cluster=$cluster[rand($s)];
  $hash{$i}=$cluster;
 }
}

###check the hash intergrity
my $check_count=0;
foreach my $i(keys %hash){
$check_count++;
}
if($check_count != $n) {die "CHECK ERROR\n;"}
####finish check

return \%hash;

}



