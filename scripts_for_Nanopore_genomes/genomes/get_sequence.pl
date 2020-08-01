#! /usr/bin/perl
$|++;
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage; 
use Bio::Perl;
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::LocatableSeq;
use Bio::SearchIO;
use List::Util qw[min max];
sub loopthroughhspNmerge();
my $file='';
#########################################
my $genomefasta='Brazil.fasta';
my $contig="CurContig6_quiver_pilon_pilon_pilon_pilon";
my $hspstart=653179	;
my $hspend=655955;
my $tctslength=$hspend-$hspstart+1;	
my $strand='-';
###########################################
my $gfffile='';
my $familyfasta='';
my $mintctslengthcutoff=300;
my $maxgaplenintcts=500;  #merge if distance <=$maxgaplenintcts

##################################			   
##get all the contig seqs
##################################
my %contigs;
my %contiglength;
my $seqio3=Bio::SeqIO->new(-format => 'Fasta', -file=>"$genomefasta");
while(my $seqobj=$seqio3->next_seq)
{
	my $id=$seqobj->id;
	my $seq=$seqobj->seq;
    $contigs{$id}=$seq;
	$contiglength{$id}=length($seq);
}


my $contigseq=$contigs{$contig};
my $tctspiece_seq=substr($contigseq,$hspstart-1,$tctslength);

if ($strand eq '-')
{
	$tctspiece_seq=revcom( $tctspiece_seq )->seq();
}

print "chr:\n";
print "start: $hspstart end: $hspend strand $strand\n";
print "$tctspiece_seq\n";