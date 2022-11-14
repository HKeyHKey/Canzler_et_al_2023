#!/usr/bin/perl
#use warnings;

$sample=$ARGV[0];

print "dans le module de comptage\n";

#Extract, for each chromosome, miRNAs coordinates (and more infos) + create reference tables (or hashes) + initialize miRNAs counts
#$miRNA_gff = 'rno.gff3'; #GFF from miRBase
$miRNA_gff = 'Coordinates_of_rno_matureOct18_on_Extended_Improved_rat_pre-miRNA.gff3'; 

open $miRNA_line, $miRNA_gff or die "Could not open $miRNA_gff: $!";
 while(<$miRNA_line>) {
 next if /^#/;
#next if /miRNA_primary_transcript/;
#$_ =~ /^(chr[^\s]*)\s.*\smiRNA\s([0-9]*)\s([0-9]*)\s\.\s([+-])\s\.\sID=(.*);Alias.*Name=(.*);.*$/ ; #GFF from miRBase
 $_ =~ /^miRNA_primary_transcript_Name=(.*)\([+-]\)\s.*miRNA\s([0-9]*)\s([0-9]*)\s\.\s([+-])\s\.\srno-([^\s]*)\s([^\s]*)\s.*$/ ; 
($pre_miRNA_name,$miRNA_start, $miRNA_stop, $miRNA_strand, $miRNA_ID, $miRNA_name) = ($1,$2,$3,$4,$6,$5);
 $miRNA_length=$miRNA_stop - $miRNA_start + 1;
# $miRNA_seq=`./Module_get_fasta.sh $miRNA_chr $miRNA_start $miRNA_stop`; #With GFF from miRBase
 $miRNA_seq=`./Module_get_fasta_from_home_made_gff.sh $pre_miRNA_name $miRNA_start $miRNA_stop` ; 
 chomp $miRNA_seq;
 print "$miRNA_name $miRNA_seq $pre_miRNA_name $miRNA_start $miRNA_stop $miRNA_length $miRNA_strand $miRNA_ID\n";
 @{$miRNA{$miRNA_ID}}=($miRNA_start,$miRNA_stop,$miRNA_length,$miRNA_strand,$miRNA_seq);
 @{$count{$miRNA_ID}}=($miRNA_name,0,0,0,0); #name - native - trimmed - tailed - trimmed then tailed
 push @{$pre_miRNA{$pre_miRNA_name}}, $miRNA_ID;
 }
close $miRNA_gff;

#Count native, trimmed and tailed miRNAs
$indice=0;
$sam='Mapping_'.$sample.'_trimmed_0.sam';
open $read, $sam or die "Could not open $sam: $!";
 while(<$read>) {
 next if /^@/;
 $indice++;
 print "sam file 0 read $indice\n";
# $_ =~ /^.*(chr[^\s]*)\s([0-9^\s]*)\s[0-9]*\s.*\s([ATGCN]*)\s.*\sAS\:i:([^\s]*)\s.*MD:Z:([0-9]*)\s.*$/ ;
 $_ =~ /^.*miRNA_primary_transcript_Name=(.*)\([+-]\)\s([0-9^\s]*)\s[0-9]*\s.*\s([ATGCN]*)\s.*\sAS\:i:([^\s]*)\s.*MD:Z:([0-9]*)\s.*$/ ;
 ($read_pre_miRNA, $read_start, $read_seq, $read_AS, $read_MD) = ($1,$2,$3,$4,$5) ;
 $read_end = $read_start + length($read_seq) - 1 + $read_AS;
  if (exists($pre_miRNA{$read_pre_miRNA})) { 
   foreach (@{%pre_miRNA{$read_pre_miRNA}}) {   
    if (@{%miRNA{$_}}[3] eq "+" && $read_start == @{%miRNA{$_}}[0]) {    
     if ($read_AS == 0) {
      if($read_end == @{%miRNA{$_}}[1] || $read_end == (@{%miRNA{$_}}[1] + 1) || $read_end == (@{%miRNA{$_}}[1] + 2) || $read_end == (@{%miRNA{$_}}[1] + 3) || $read_end == (@{%miRNA{$_}}[1] - 1)) {#5' 0 & 3' -1 0 +1 +2 +3 nt = native   
      @{%count{$_}}[1]++;  
      } elsif ($read_end == (@{%miRNA{$_}}[1] - 2) || $read_end == (@{%miRNA{$_}}[1] - 3) || $read_end == (@{%miRNA{$_}}[1] - 4) || $read_end == (@{%miRNA{$_}}[1] - 5)) {  #5' 0 & 3' -2 -3 -4 -5 nt = trimmed
      @{%count{$_}}[2]++;
      }    
     } elsif ($read_AS != 0) {
      if ((substr $read_seq,0,(@{%miRNA{$_}}[2] -1)) eq (substr @{%miRNA{$_}}[4],0,(@{%miRNA{$_}}[2] -1))) {
      @{%count{$_}}[3]++;
      } elsif ((substr $read_seq,0,(@{%miRNA{$_}}[2] -5)) eq (substr @{%miRNA{$_}}[4],0,(@{%miRNA{$_}}[2] -5))) {
      @{%count{$_}}[4]++;
      }
     }
    } elsif (@{%miRNA{$_}}[3] eq "-" && $read_end == @{%miRNA{$_}}[1]) {
     if ($read_AS == 0) {
      if ($read_start == @{%miRNA{$_}}[0] || $read_start == (@{%miRNA{$_}}[0] - 1) || $read_start == (@{%miRNA{$_}}[0] - 2) || $read_start == (@{%miRNA{$_}}[0] - 3) || $read_start == (@{%miRNA{$_}}[0] + 1)) {  #5' 0 & 3' -1 0 +1 +2 +3 nt = native     
      @{%count{$_}}[1]++;
      } elsif ($read_start == (@{%miRNA{$_}}[0] + 2) || $read_start == (@{%miRNA{$_}}[0] + 3) || $read_start == (@{%miRNA{$_}}[0] + 4) || $read_start == (@{%miRNA{$_}}[0] + 5)) {  #5' 0 & 3' -2 -3 -4 -5 nt = trimmed    
      @{%count{$_}}[2]++;
      }
     } elsif ($read_AS != 0) { 
      if ((substr $read_seq,-(@{%miRNA{$_}}[2] -1)) eq (substr @{%miRNA{$_}}[4],-(@{%miRNA{$_}}[2] -1))) {
      @{%count{$_}}[3]++;
      } elsif ((substr $read_seq,-(@{%miRNA{$_}}[2] -5)) eq (substr @{%miRNA{$_}}[4],-(@{%miRNA{$_}}[2] -5))) {
      @{%count{$_}}[4]++;
      }     
     }
    }
   }
  }
 } 
close $sam;

#Count tailed and trimmed + tailed miRNAs left
for $i (1..5) {
$indice=0;
$sam='Mapping_'.$sample.'_trimmed_'.$i.'.sam';
 open $read, $sam or die "Could not open $sam: $!";
  while(<$read>) {
  next if /^@/;
  $indice++;
  print "sam file $i read $indice\n";
#  $_ =~ /^.*(chr[^\s]*)\s([0-9^\s]*)\s[0-9]*\s.*\s([ATGCN]*)\s.*\sAS\:i:([^\s]*)\s.*MD:Z:([0-9]*)\s.*$/ ;
  $_ =~ /^.*miRNA_primary_transcript_Name=(.*)\([+-]\)\s([0-9^\s]*)\s[0-9]*\s.*\s([ATGCN]*)\s.*\sAS\:i:([^\s]*)\s.*MD:Z:([0-9]*)\s.*$/;
  ($read_pre_miRNA, $read_start, $read_seq, $read_AS, $read_MD) = ($1,$2,$3,$4,$5) ;
  $read_end = $read_start + length($read_seq) - 1 + $read_AS;
  if (exists($pre_miRNA{$read_pre_miRNA})) { 
   foreach (@{%pre_miRNA{$read_pre_miRNA}}) {   
    if (@{%miRNA{$_}}[3] eq "+" && $read_start == @{%miRNA{$_}}[0]) {    
     if ($read_AS == 0) {
      if($read_end == @{%miRNA{$_}}[1] || $read_end == (@{%miRNA{$_}}[1] + 1) || $read_end == (@{%miRNA{$_}}[1] + 2) || $read_end == (@{%miRNA{$_}}[1] + 3) || $read_end == (@{%miRNA{$_}}[1] - 1)) {#5' 0 & 3' -1 0 +1 +2 +3 nt = native   
      @{%count{$_}}[3]++;  
      } elsif ($read_end == (@{%miRNA{$_}}[1] - 2) || $read_end == (@{%miRNA{$_}}[1] - 3) || $read_end == (@{%miRNA{$_}}[1] - 4) || $read_end == (@{%miRNA{$_}}[1] - 5)) {  #5' 0 & 3' -2 -3 -4 -5 nt = trimmed
      @{%count{$_}}[4]++;
      }    
     } elsif ($read_AS != 0) {
      if ((substr $read_seq,0,(@{%miRNA{$_}}[2] -1)) eq (substr @{%miRNA{$_}}[4],0,(@{%miRNA{$_}}[2] -1))) {
      @{%count{$_}}[3]++;
      } elsif ((substr $read_seq,0,(@{%miRNA{$_}}[2] -5)) eq (substr @{%miRNA{$_}}[4],0,(@{%miRNA{$_}}[2] -5))) {
      @{%count{$_}}[4]++;
      }
     }
    } elsif (@{%miRNA{$_}}[3] eq "-" && $read_end == @{%miRNA{$_}}[1]) {
     if ($read_AS == 0) {
      if ($read_start == @{%miRNA{$_}}[0] || $read_start == (@{%miRNA{$_}}[0] - 1) || $read_start == (@{%miRNA{$_}}[0] - 2) || $read_start == (@{%miRNA{$_}}[0] - 3) || $read_start == (@{%miRNA{$_}}[0] + 1)) {  #5' 0 & 3' -1 0 +1 +2 +3 nt = native     
      @{%count{$_}}[3]++;
      } elsif ($read_start == (@{%miRNA{$_}}[0] + 2) || $read_start == (@{%miRNA{$_}}[0] + 3) || $read_start == (@{%miRNA{$_}}[0] + 4) || $read_start == (@{%miRNA{$_}}[0] + 5)) {  #5' 0 & 3' -2 -3 -4 -5 nt = trimmed    
      @{%count{$_}}[4]++;
      }
     } elsif ($read_AS != 0) { 
      if ((substr $read_seq,-(@{%miRNA{$_}}[2] -1)) eq (substr @{%miRNA{$_}}[4],-(@{%miRNA{$_}}[2] -1))) {
      @{%count{$_}}[3]++;
      } elsif ((substr $read_seq,-(@{%miRNA{$_}}[2] -5)) eq (substr @{%miRNA{$_}}[4],-(@{%miRNA{$_}}[2] -5))) {
      @{%count{$_}}[4]++;
      }     
     }
    }
   }
  }
 }
 close $sam;
}

#Write output file
open($output,'>', 'count_report_'.$sample.'.csv');
#print $output "miRNA_ID\tmiRNA_name\tnative_count\ttrimmed_count\ttailed_count\ttrimmed_and_tailed_count\n"; 
foreach $key (keys %count)
{
  @value = @{$count{$key}};
  print $output "$key\t",join("\t",@value),"\n";
}
close $output;

