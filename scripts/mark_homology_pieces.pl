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
my $genomefasta='';
my $gfffile='';
my $familyfasta='';
my $mintctslengthcutoff=300;
my $maxgaplenintcts=500;  #merge if distance <=$maxgaplenintcts
my $message_text  = "please specify option parameters. -i <blastout file> -g <genome fasta> -maxhspgap <num> -minTcTSLength <num> -gff <gff file> -familyfasta <fasta file>";
  my $exit_status   = 2;          ## The exit status to use
  my $verbose_level = 0;          ## The verbose level to use
  my $filehandle    = \*STDERR;   ## The filehandle to write to
GetOptions ('i=s' => \$file,'g=s' => \$genomefasta,'maxhspgap=s' => \$maxgaplenintcts,'minTcTSLength=s' => \$mintctslengthcutoff, 'gff=s' => \$gfffile, 'familyfasta=s' => \$familyfasta);
pod2usage( { -message => $message_text ,
               -exitval => $exit_status  ,
               -verbose => $verbose_level,
               -output  => $filehandle } ) if ($file eq '' or $genomefasta eq '');

open STATOUT, ">stats.final.txt";

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
##################################
##get ts coordinates from gff file
##################################
open GFF, "$gfffile";
my %coords;
	while (my $line=<GFF>) # get all ID' location
{
	if ($line=~/gene/)
	{
		if ($line=~/ID=(.+?)\;/)
		{
		##get loc
		my @line=split("\t",$line);

		my $id=$1;
		$line=~/Name=(.+?)\;/;
		my $desc=$1;
	    #print "$line[0] $line[3] $line[4] $id $desc\n";
		  $coords{$id}{"chr"}=$line[0];
          $coords{$id}{"start"}=$line[3];
		  $coords{$id}{"end"}=$line[4];
		  $coords{$id}{"strand"}=$line[6];
		  $coords{$id}{"desc"}=$desc;
		}
	}
}

##################################
##get seq from fasta file
##################################
my %familyfasta;
my $seqio4=Bio::SeqIO->new(-format => 'Fasta', -file=>"$familyfasta");
my $inputcount=0;
while(my $seqobj=$seqio4->next_seq)
{
	my $id=$seqobj->id;
	#print($id."\n");
	#$id=~/Parent=(.+?)\|/; #extract |Parent=a2dd0951-a251-49ca-a2de-e65db7e314ac-CDS:gene|
	#print($1."\n");
	#my $extracted_id=$1;
	my $extracted_id=$id;
	my $seq=$seqobj->seq;
    $familyfasta{$extracted_id}=$seq;
	$inputcount+=1
}
print STATOUT "$inputcount\n";
####################################################################
######parsing the blastout file and get the best match#######
####################################################################
my %markedcontigs;
my $in = new Bio::SearchIO(-format => 'blast',
                           -file   => "$file");
my $maxmatchlen=0;
open DIAG, ">diag.txt";

while( my $result = $in->next_result ) {

my $q_len = $result->query_length;
my $queryname=$result->query_name;
#print "result: ".$queryname."\n";

while(my $hit = $result->next_hit) {


		my $hitname=$hit->name;
        my $hitdesc=$hit->description;
		#print "hit: ".$hitname."\n";

	while( my $hsp = $hit->next_hsp) {

			my $hspstrand=$hsp->strand('hit');
			my $hspstart=$hsp->start('hit');
			my $hspend=$hsp->end('hit');
			my $hsplength=$hspstart-$hspend;
			print DIAG "$queryname\t$hitname\t$hitdesc\t$hspstart\t$hspend\t$hspstrand\n";
		my $largestkeynum=0;
		if (exists $markedcontigs{$hitname})
		{
				my $v1=$markedcontigs{$hitname};
				my %v1=%$v1;
				my $newhspflag=1;

			LOOPOVERHSP: foreach my $keynum (sort keys %v1)  # go through all hash-stored hsp groups in current contig (hitname)
			{
				my $keystart=$v1{$keynum}{'start'};
				my $keyend=$v1{$keynum}{'end'};
				my $keyhsplength=$keyend-$keystart;

				if ($hspend<=$keyend and $hspstart>=$keystart) # hsp is bracketed within marked region
				{
					$newhspflag=0;
					#print "bracketed by $keynum $keystart $keyend\n";
					last LOOPOVERHSP;
				}
				if ($hspend>=$keyend and $hspstart<=$keystart and $hsplength>=(100+$keyhsplength))   #  marked region is bracketed within hsp, and the difference is >=100bp
				{
					$v1{$keynum}{'start'}=$hspstart;
					$v1{$keynum}{'end'}=$hspend;
					$markedcontigs{$hitname}{$keynum}{'start'}=$hspstart;
					$markedcontigs{$hitname}{$keynum}{'end'}=$hspend;
					$markedcontigs{$hitname}{$keynum}{'strand'}=$hspstrand;
					$newhspflag=0;
					#print "bracketing $keynum $keystart $keyend\n";
					last LOOPOVERHSP;
				}
				if ($hspend>=$keyend and $hspstart<=$keyend and ($hspend-$keyend)>=100) # hsp overlaps marked region on regions's right side, and overlap is >=100bp
				{
					$v1{$keynum}{'end'}=$hspend;
					$markedcontigs{$hitname}{$keynum}{'end'}=$hspend;
					$newhspflag=0;
					##determine strand of the modified markedcontigs
					if (($hspend-$hspstart)>=($keyend-$keystart) ) ####update strand if hsp is longer#####
					{
						$markedcontigs{$hitname}{$keynum}{'strand'}=$hspstrand;
					}

					#print "hsp is overlapping the right of $keynum $keystart $keyend\n";
					last LOOPOVERHSP;
				}
				if ($hspend>=$keystart and $hspstart<=$keystart and ($keystart-$hspstart)>=100) # hsp overlaps marked region on region's left side, and overlap is >=100bp
				{
					$v1{$keynum}{'start'}=$hspstart;
					$markedcontigs{$hitname}{$keynum}{'start'}=$hspstart;
					$newhspflag=0;
					##determine strand of the modified markedcontigs
					if (($hspend-$hspstart)>=($keyend-$keystart) ) ####update strand if hsp is longer#####
					{
						$markedcontigs{$hitname}{$keynum}{'strand'}=$hspstrand;
					}
					#print "hsp is overlapping the left of $keynum $keystart $keyend\n";
					last LOOPOVERHSP;
				}
				if ($hspend<=$keystart and $hspend>=($keystart-$maxgaplenintcts) and $hsplength>=100) # hsp is on region's left side, merge if distance <=maxgaplenintcts, and hsp>=100bp
				{
					$v1{$keynum}{'start'}=$hspstart;
					$markedcontigs{$hitname}{$keynum}{'start'}=$hspstart;
					$newhspflag=0;
					##determine strand of the modified markedcontigs
					if (($hspend-$hspstart)>=($keyend-$keystart) ) ####update strand if hsp is longer#####
					{
						$markedcontigs{$hitname}{$keynum}{'strand'}=$hspstrand;
					}
					#print "hsp is on the left of $keynum $keystart $keyend\n";
					last LOOPOVERHSP;
				}
				if ($hspstart>=$keyend and $hspstart<=($keyend+$maxgaplenintcts) and $hsplength>=100) # hsp is on region's right side, merge if distance <=maxgaplenintcts, and hsp>=100bp
				{
					$v1{$keynum}{'end'}=$hspend;
					$markedcontigs{$hitname}{$keynum}{'end'}=$hspend;
					$newhspflag=0;
					##determine strand of the modified markedcontigs
					if (($hspend-$hspstart)>=($keyend-$keystart) ) ####update strand if hsp is longer#####
					{
						$markedcontigs{$hitname}{$keynum}{'strand'}=$hspstrand;
					}
					#print "hsp is on the right of $keynum $keystart $keyend\n";
					last LOOPOVERHSP;
				}

			}
			my $keynummax=0;
			foreach my $keynum (sort keys %v1)# get the largest key #
			{
				$keynummax=$keynum if $keynum>=$keynummax;
			}
			#print "largest key num: $keynummax\n";
			if ($newhspflag==1) # mark new region on hit
			{
					my $newkeynum=$keynummax+1;
					$markedcontigs{$hitname}{$newkeynum}{'start'}=$hspstart;
					$markedcontigs{$hitname}{$newkeynum}{'end'}=$hspend;
					$markedcontigs{$hitname}{$newkeynum}{'strand'}=$hspstrand;
					#print "this hsp $hspstart $hspend is a new hsp on $hitname\n";
			}
		}
		else{  #initialize hash, mark a new region on hit
			$markedcontigs{$hitname}{1}{'start'}=$hspstart;
			$markedcontigs{$hitname}{1}{'end'}=$hspend;
			$markedcontigs{$hitname}{1}{'strand'}=$hspstrand;
			#print "initializing the first hsp $hspstart $hspend on $hitname\n";
		}

		}
	}
}

my $markedregioncount= 0;
foreach my $hitname (sort keys %markedcontigs)
	{
		my $v1=$markedcontigs{$hitname};
		my %v1=%$v1;
		foreach my $keynum (sort keys %v1)  # go through all hash-stored hsp groups in current contig (hitname)
			{
				$markedregioncount+=1;
			}
	}
print "The count of marked regions is $markedregioncount\n";
print STATOUT "$markedregioncount\n";

##finishing parsing the blastout results
##merging hsps again
#loopthroughhspNmerge();
#loopthroughhspNmerge();

####################################
#merge tcts with same start or end
# 20190305 ATTENTION: NEED TO MERGE NESTED GENES, these could be due to different query matching to nested regions, these will get past the hsp merging step above
####################################
my %markedcontigs_merged;
print "start to merge genes\n";
foreach my $contig (sort keys %markedcontigs)
{
	my $v1=$markedcontigs{$contig};
	my %v1=%$v1;
	#print $contig."\n";
	foreach my $keynum1 (sort keys %v1) #going through all hsps in current contig
	{
				#print $keynum1."\n";
				my $hspstart=$v1{$keynum1}{'start'};
				my $hspend=$v1{$keynum1}{'end'};
				my $hspstrand=$v1{$keynum1}{'strand'};
				my $hsplength=($hspend-$hspstart+1);
				my $merged_flag=0;


				for my $contig_merged (sort keys %markedcontigs_merged) #go through merged ts
				{
					my $z1=$markedcontigs_merged{$contig_merged};
					my %z1=%$z1;
					#print $contig."\n";
					next if $contig ne $contig_merged;
					#print "try merging gene on same chr\n";
					#print $contig."\n";
					#print "$hspstart $hspend\n";
					foreach my $keynum2 (sort keys %z1) #going through all hsps in current contig
					{

						my $hspstart_merged=$markedcontigs_merged{$contig_merged}{$keynum2}{'start'};
						my $hspend_merged=$markedcontigs_merged{$contig_merged}{$keynum2}{'end'};
						my $hspstrand_merged=$markedcontigs_merged{$contig_merged}{$keynum2}{'strand'};
						my $hsplength_merged=$hspend_merged-$hspstart_merged;
						###near identical start and end###
						if ($hspstart<=$hspstart_merged and $hspstart>=$hspstart_merged-10
						   or $hspend>=$hspend_merged and $hspend<=$hspstart_merged+10 )
						{
							my $new_start=min($hspstart,$hspstart_merged);
							my $new_end=max($hspend,$hspend_merged);
							#update start and end
							$markedcontigs_merged{$contig_merged}{$keynum2}{'start'}=$new_start;
							$markedcontigs_merged{$contig_merged}{$keynum2}{'end'}=$new_end;
							#update strand, take the strand of the longer gene
							if (($hspend-$hspstart)>=($hspend_merged-$hspstart_merged) ) ####update strand if hsp is longer#####
							{
								$markedcontigs_merged{$contig_merged}{$keynum2}{'strand'}=$hspstrand;
							}

							$merged_flag=1;
							#print "merging $hspstart $hspend , $hspstart_merged $hspend_merged, result $new_start $new_end\n";
						}
						### hsp is bracketed within merged region ###
						if ($hspstart_merged<=$hspstart and $hspend<=$hspend_merged )
						{
							$merged_flag=1;
							#print "bracketed by $keynum $keystart $keyend\n";
						}
						## merged region is bracketed within hsp, and the difference is >=50bp
						if (($hspstart<=$hspstart_merged and $hspend_merged<=$hspend))
						{
							$markedcontigs_merged{$contig_merged}{$keynum2}{'start'}=$hspstart;
							$markedcontigs_merged{$contig_merged}{$keynum2}{'end'}=$hspend;
							$markedcontigs_merged{$contig_merged}{$keynum2}{'strand'}=$hspstrand;
							$merged_flag=1;
						}
						### hsp overlaps merged region on regions's right side, and overlap is >=50bp
						if ($hspend>=$hspend_merged and $hspstart<=$hspend_merged and ($hspend_merged-$hspstart)>=50)
						{
							$markedcontigs_merged{$contig_merged}{$keynum2}{'end'}=$hspend;
							##determine strand of the modified markedcontigs
							if ($hsplength>=$hsplength_merged ) ####update strand if hsp is longer#####
							{
								$markedcontigs_merged{$contig_merged}{$keynum2}{'strand'}=$hspstrand;
							}
							$merged_flag=1;
						}
						### hsp overlaps marked region on region's left side, and overlap is >=50bp
						if ($hspend>=$hspstart_merged and $hspstart<=$hspstart_merged and ($hspend-$hspstart_merged)>=50)
						{
							$markedcontigs_merged{$contig_merged}{$keynum2}{'start'}=$hspstart;
							##determine strand of the modified markedcontigs
							if ($hsplength>=$hsplength_merged ) ####update strand if hsp is longer#####
							{
								$markedcontigs_merged{$contig_merged}{$keynum2}{'strand'}=$hspstrand;
							}
							#print "MERGING, hsp is overlapping the left of $keynum $keystart $keyend\n";
							$merged_flag=1;
						}
						### hsp is on region's left side, merge if distance <=maxgaplenintcts, and hsp>=100bp
						if ($hspend<=$hspstart_merged and $hspend>=($hspstart_merged-$maxgaplenintcts) and $hsplength>=100)
						{
							$markedcontigs_merged{$contig_merged}{$keynum2}{'start'}=$hspstart;
							##determine strand of the modified markedcontigs
							if ($hsplength>=$hsplength_merged ) ####update strand if hsp is longer#####
							{
								$markedcontigs_merged{$contig_merged}{$keynum2}{'strand'}=$hspstrand;
							}
							#print "MERGING, hsp is on the left of $keynum $keystart $keyend\n";
							$merged_flag=1;
						}
						### hsp is on region's right side, merge if distance <=maxgaplenintcts, and hsp>=100bp
						if ($hspstart>=$hspend_merged and $hspstart<=($hspend_merged+$maxgaplenintcts) and $hsplength>=100)
						{
							$markedcontigs_merged{$contig_merged}{$keynum2}{'end'}=$hspend;
							##determine strand of the modified markedcontigs
							if ($hsplength>=$hsplength_merged ) ####update strand if hsp is longer#####
							{
								$markedcontigs_merged{$contig_merged}{$keynum2}{'strand'}=$hspstrand;
							}
							#print "MERGING, hsp is on the right of $keynum $keystart $keyend\n";
							$merged_flag=1;
						}

					}
				}

				if ($merged_flag==0) # add new gene
				{
					$markedcontigs_merged{$contig}{$keynum1}{'start'}=$hspstart;
					$markedcontigs_merged{$contig}{$keynum1}{'end'}=$hspend;
					$markedcontigs_merged{$contig}{$keynum1}{'strand'}=$hspstrand;

				}

				##...

	}

}
my $mergedregioncount= 0;
foreach my $hitname (sort keys %markedcontigs_merged)
	{
		my $v1=$markedcontigs_merged{$hitname};
		my %v1=%$v1;
		foreach my $keynum (sort keys %v1)  # go through all hash-stored hsp groups in current contig (hitname)
			{
				$mergedregioncount+=1;
			}
	}
print "The count of marked regions is $mergedregioncount\n";

####################################
#merge tcts with same start or end
# second merge
####################################
my %markedcontigs_merged2;
print "start to merge genes\n";
foreach my $contig (sort keys %markedcontigs_merged)
{
	my $v1=$markedcontigs_merged{$contig};
	my %v1=%$v1;
	#print $contig."\n";
	foreach my $keynum1 (sort keys %v1) #going through all hsps in current contig
	{
				#print $keynum1."\n";
				my $hspstart=$v1{$keynum1}{'start'};
				my $hspend=$v1{$keynum1}{'end'};
				my $hspstrand=$v1{$keynum1}{'strand'};
				my $hsplength=($hspend-$hspstart+1);
				my $merged_flag=0;


				for my $contig_merged (sort keys %markedcontigs_merged2) #go through merged ts
				{
					my $z1=$markedcontigs_merged2{$contig_merged};
					my %z1=%$z1;
					#print $contig."\n";
					next if $contig ne $contig_merged;
					#print "try merging gene on same chr\n";
					#print $contig."\n";
					#print "$hspstart $hspend\n";
					foreach my $keynum2 (sort keys %z1) #going through all hsps in current contig
					{

						my $hspstart_merged=$markedcontigs_merged2{$contig_merged}{$keynum2}{'start'};
						my $hspend_merged=$markedcontigs_merged2{$contig_merged}{$keynum2}{'end'};
						my $hspstrand_merged=$markedcontigs_merged2{$contig_merged}{$keynum2}{'strand'};
						my $hsplength_merged=$hspend_merged-$hspstart_merged;
						###near identical start and end###
						if ($hspstart<=$hspstart_merged and $hspstart>=$hspstart_merged-10
						   or $hspend>=$hspend_merged and $hspend<=$hspstart_merged+10 )
						{
							my $new_start=min($hspstart,$hspstart_merged);
							my $new_end=max($hspend,$hspend_merged);
							#update start and end
							$markedcontigs_merged2{$contig_merged}{$keynum2}{'start'}=$new_start;
							$markedcontigs_merged2{$contig_merged}{$keynum2}{'end'}=$new_end;
							#update strand, take the strand of the longer gene
							if (($hspend-$hspstart)>=($hspend_merged-$hspstart_merged) ) ####update strand if hsp is longer#####
							{
								$markedcontigs_merged2{$contig_merged}{$keynum2}{'strand'}=$hspstrand;
							}

							$merged_flag=1;
							#print "merging $hspstart $hspend , $hspstart_merged $hspend_merged, result $new_start $new_end\n";
						}
						### hsp is bracketed within merged region ###
						if ($hspstart_merged<=$hspstart and $hspend<=$hspend_merged )
						{
							$merged_flag=1;
							#print "bracketed by $keynum $keystart $keyend\n";
						}
						## merged region is bracketed within hsp, and the difference is >=50bp
						if (($hspstart<=$hspstart_merged and $hspend_merged<=$hspend))
						{
							$markedcontigs_merged2{$contig_merged}{$keynum2}{'start'}=$hspstart;
							$markedcontigs_merged2{$contig_merged}{$keynum2}{'end'}=$hspend;
							$markedcontigs_merged2{$contig_merged}{$keynum2}{'strand'}=$hspstrand;
							$merged_flag=1;
						}
						### hsp overlaps merged region on regions's right side, and overlap is >=50bp
						if ($hspend>=$hspend_merged and $hspstart<=$hspend_merged and ($hspend_merged-$hspstart)>=50)
						{
							$markedcontigs_merged2{$contig_merged}{$keynum2}{'end'}=$hspend;
							##determine strand of the modified markedcontigs
							if ($hsplength>=$hsplength_merged ) ####update strand if hsp is longer#####
							{
								$markedcontigs_merged2{$contig_merged}{$keynum2}{'strand'}=$hspstrand;
							}
							$merged_flag=1;
						}
						### hsp overlaps marked region on region's left side, and overlap is >=50bp
						if ($hspend>=$hspstart_merged and $hspstart<=$hspstart_merged and ($hspend-$hspstart_merged)>=50)
						{
							$markedcontigs_merged2{$contig_merged}{$keynum2}{'start'}=$hspstart;
							##determine strand of the modified markedcontigs
							if ($hsplength>=$hsplength_merged ) ####update strand if hsp is longer#####
							{
								$markedcontigs_merged2{$contig_merged}{$keynum2}{'strand'}=$hspstrand;
							}
							#print "hsp is overlapping the left of $keynum $keystart $keyend\n";
							$merged_flag=1;
						}
						### hsp is on region's left side, merge if distance <=maxgaplenintcts, and hsp>=100bp
						if ($hspend<=$hspstart_merged and $hspend>=($hspstart_merged-$maxgaplenintcts) and $hsplength>=100)
						{
							$markedcontigs_merged2{$contig_merged}{$keynum2}{'start'}=$hspstart;
							##determine strand of the modified markedcontigs
							if ($hsplength>=$hsplength_merged ) ####update strand if hsp is longer#####
							{
								$markedcontigs_merged2{$contig_merged}{$keynum2}{'strand'}=$hspstrand;
							}
							#print "hsp is on the left of $keynum $keystart $keyend\n";
							$merged_flag=1;
						}
						### hsp is on region's right side, merge if distance <=maxgaplenintcts, and hsp>=100bp
						if ($hspstart>=$hspend_merged and $hspstart<=($hspend_merged+$maxgaplenintcts) and $hsplength>=100)
						{
							$markedcontigs_merged2{$contig_merged}{$keynum2}{'end'}=$hspend;
							##determine strand of the modified markedcontigs
							if ($hsplength>=$hsplength_merged ) ####update strand if hsp is longer#####
							{
								$markedcontigs_merged2{$contig_merged}{$keynum2}{'strand'}=$hspstrand;
							}
							#print "hsp is on the right of $keynum $keystart $keyend\n";
							$merged_flag=1;
						}

					}
				}

				if ($merged_flag==0) # add new gene
				{
					$markedcontigs_merged2{$contig}{$keynum1}{'start'}=$hspstart;
					$markedcontigs_merged2{$contig}{$keynum1}{'end'}=$hspend;
					$markedcontigs_merged2{$contig}{$keynum1}{'strand'}=$hspstrand;

				}

				##...

	}

}
my $mergedregioncount2= 0;
foreach my $hitname (sort keys %markedcontigs_merged2)
	{
		my $v1=$markedcontigs_merged2{$hitname};
		my %v1=%$v1;
		foreach my $keynum (sort keys %v1)  # go through all hash-stored hsp groups in current contig (hitname)
			{
				$mergedregioncount2+=1;
			}
	}
print "The count of marked regions is $mergedregioncount2\n";
print STATOUT "$mergedregioncount2\n";


######################
## out put TcTS pieces
######################
my %printed_genes;
my @genomefasta=split (/\//,$genomefasta);
my $genomename =	$genomefasta[$#genomefasta];
#print $genomename."\n";
#print "start to output genes\n";
open OUT, ">$file.mintctslengthcutoff$mintctslengthcutoff.maxgaplenintcts$maxgaplenintcts.fasta" or die; # open file to store result
foreach my $contig (sort keys %markedcontigs_merged2)
{
	my $v1=$markedcontigs_merged2{$contig};
	my %v1=%$v1;
	my $contigseq=$contigs{$contig};
	#print $contig."\n";
	foreach my $keynum1 (sort keys %v1) #going through all hsps in current contig
	{
				#print $keynum1."\n";
				my $hspstart=$v1{$keynum1}{'start'};
				my $hspend=$v1{$keynum1}{'end'};
				my $hsplength=($hspend-$hspstart+1);
				next if $hsplength<=$mintctslengthcutoff; # length cutoff check
				##check if current ts is a known ts
				##insert code here
				##...
				#print "${contig}_${hspstart}_${hspend}\n";
				#go through all ts to check for overlapping with known ts
				my $overlap_flag=0;
				my $existing_ID=""; # store ID if current ts overlaps with a ID
				my $overlap_size=0; # overlap size, for comparing each overlap with
				for my $known_ID (sort keys %coords)
				{
					next if $contig ne $coords{$known_ID}{"chr"};
					#print "on same chr\n";
					my $ts_start=$coords{$known_ID}{"start"};
					my $ts_end=$coords{$known_ID}{"end"};

					next if $hspstart<=$ts_start and $hspend>=$ts_end ;             #hsp brackets ts, use hsp's coord, and mask the known ts gene

					if (   ($hspstart>=($ts_start+10) and $hspstart<=($ts_end-10) and (($ts_end-$hspstart)/$hsplength)>=0.9)    #hsp start is in ts, and overlap is larger than 90% length of hsp
						or ($hspend>=($ts_start+10) and $hspend<=($ts_end-10)  and (($hspend-$ts_start)/$hsplength)>=0.9)       #hsp end is in ts, and overlap is larger than 90% length of hsp
						or ($ts_start<=$hspstart and $ts_end>=$hspend )             #hsp is bracked by ts

						)
					{
						#print "overlap found\nts:${ts_start}-$ts_end\thsp:${hspstart}-$hspend $known_ID\n";
						$overlap_flag=1;
						##calculate current_overlap_size
						my $current_overlap_size=0;
						$current_overlap_size= $ts_end-$hspstart if ($hspstart>=($ts_start) and $hspstart<=($ts_end));
						$current_overlap_size= $hspend-$ts_start if ($hspend>=($ts_start) and $hspend<=($ts_end));
						#....................
						if ($current_overlap_size>=$overlap_size)
						{
							$existing_ID=$known_ID;
						}
					}

				}
				##...
				if ($overlap_flag==1) #existing family member
				{
					if (exists $familyfasta{$existing_ID}) #current ID is in familyfasta, thus directly output sequence
					{
						print OUT ">".$existing_ID."\n".$familyfasta{$existing_ID}."\n" if !(exists $printed_genes{$existing_ID}); # print out if not printed current ID yet
						$printed_genes{$existing_ID}=1
					}
					else #no sequence, this maybe a pseudogene or not a true family member , need to fetch sequence
					{
						#print "NO SEQUENCE FOR $existing_ID\n" ;

						##code to fetch sequence here
						my $contig_2_fetch=$coords{$existing_ID}{"chr"};
						my $start_2_fetch=$coords{$existing_ID}{"start"};
						my $end_2_fetch=$coords{$existing_ID}{"end"};
						my $length=$end_2_fetch-$start_2_fetch+1;
						my $strand_2_fetch=$coords{$existing_ID}{"strand"};
						my $contig_seq=$contigs{$contig_2_fetch};
						my $gene_seq=substr($contig_seq,$start_2_fetch-1,$length);

						if ($strand_2_fetch eq "+")
						{
							print OUT ">".$existing_ID."\n".$gene_seq."\n" if !(exists $printed_genes{$existing_ID});
						}
						else{ # - strand, revcom
							my $gene_seq_revcom = revcom( $gene_seq )->seq();
							print OUT ">".$existing_ID."\n".$gene_seq_revcom."\n" if !(exists $printed_genes{$existing_ID});
						}

						#if !(exists $printed_genes{$existing_ID}); # print out if not printed current ID yet
						$printed_genes{$existing_ID}=1
					}
				}
				else #not existing family member
				{
					my $tctspiece_seq=substr($contigseq,$hspstart-1,$hsplength);
					my $tctspiece_strand=$v1{$keynum1}{'strand'};
					$tctspiece_seq = revcom( $tctspiece_seq )->seq() if $tctspiece_strand==-1; #revcom if on the minus strand
					my $tctspiece_id="${contig}_${hspstart}_${hspend}_${tctspiece_strand}";
					print OUT ">$tctspiece_id\n$tctspiece_seq\n" if !(exists $printed_genes{$tctspiece_id});
					$printed_genes{$tctspiece_id}=1;
				}
	}

}

###########################end of main program##########################################

##function not used
sub loopthroughhspNmerge()
{
foreach my $contig (sort keys %markedcontigs)
{
	my $v1=$markedcontigs{$contig};
	my %v1=%$v1;
	foreach my $keynum1 (sort keys %v1) #going through all hsps in current contig
	{
				my $hspstart=$v1{$keynum1}{'start'};
				my $hspend=$v1{$keynum1}{'end'};

		LOOPOVERHSP: foreach my $keynum2 (sort keys %v1) #going through all hsps in current contig
		{
			next if $keynum2 eq $keynum1;
			my $keystart=$v1{$keynum2}{'start'};
			my $keyend=$v1{$keynum2}{'end'};

				if ($hspend<=$keyend and $hspstart>=$keystart) # hsp is bracketed within marked region
				{

					#print "bracketed by $keynum2 $keystart $keyend\n";
					delete $markedcontigs{$contig}{$keynum1};
					delete $v1{$keynum1};
					last LOOPOVERHSP;
				}
				if ($hspend>=$keyend and $hspstart<=$keystart)   #  marked region is brackted within hsp
				{
					delete $markedcontigs{$contig}{$keynum2};
					delete $v1{$keynum2};
					#print "bracketing $keynum2 $keystart $keyend\n";
					last LOOPOVERHSP;
				}
				if ($hspend>=$keyend and $hspstart<=$keyend) # hsp overlaps marked region on regions's right side
				{
					$v1{$keynum2}{'end'}=$hspend;
					delete $markedcontigs{$contig}{$keynum1};
					delete $v1{$keynum1};
					#print "hsp is overlapping the right of $keynum2 $keystart $keyend\n";
					last LOOPOVERHSP;
				}
				if ($hspend>=$keystart and $hspstart<=$keystart) # hsp overlaps marked region on region's left side
				{
					$v1{$keynum2}{'start'}=$hspstart;
					delete $markedcontigs{$contig}{$keynum1};
					delete $v1{$keynum1};
					#print "hsp is overlapping the left of $keynum2 $keystart $keyend\n";
					last LOOPOVERHSP;
				}
				if ($hspend<=$keystart and $hspend>=($keystart-300)) # hsp is on region's left side, merge if distance <=300
				{
					$v1{$keynum2}{'start'}=$hspstart;
					delete $markedcontigs{$contig}{$keynum1};
					delete $v1{$keynum1};
					#print "hsp is on the left of $keynum2 $keystart $keyend\n";
					last LOOPOVERHSP;
				}
				if ($hspstart>=$keyend and $hspstart<=($keyend+300)) # hsp is on region's right side, merge if distance <=300
				{
					$v1{$keynum2}{'end'}=$hspend;
					delete $markedcontigs{$contig}{$keynum1};
					delete $v1{$keynum1};
					#print "hsp is on the right of $keynum2 $keystart $keyend\n";
					last LOOPOVERHSP;
				}
	}
}

}
}
