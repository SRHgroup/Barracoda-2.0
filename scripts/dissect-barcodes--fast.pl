#!/usr/bin/perl

# Andrea Marquard
# Mar 16, 2015
#
#
# Usage: <script.pl> <merged-sam>
#
#
# Aim: take a merged sam file, and locate all the sub-sequences within each read (barcodes, primers etc).
#      Also map the sample sequence to a sample ID.
#
#
# Output: prints the following tab-delimited fields to STDOUT (and a header):
#
# 1)   read_id         id of the read   
# 2)   read_quality    mean PHRED quality of the read        
# 3)   sample_id       id of the sample barcode     
# 4)   nw_score        global alignment score from Needleman-Wunsch when identifying the sample ID
# 5)   pepA            id of the A peptide barcode
# 6)   pepB            id of the B peptide barcode
# 7)   sample_seq      the actual sequence of the found sample barcode      
# 8)   primA_seq       the actual sequence of the found A primer     
# 9)   N6A             the actual sequence of the found A N6
# 10)  pepA_seq        the actual sequence of the found A peptide barcode    
# 11)  anneal_seq      the actual sequence of the found annealing region     
# 12)  pepB_seq        the actual sequence of the found B peptide barcode    
# 13)  N6B             the actual sequence of the found B N6
# 14)  primB_seq       the actual sequence of the found B primer     
#
#
# Input description:
#
# A "merged sam file" has the following tab-delimited fields:
# 1) read id
# 2) id of the peptide A barcode that the read aligned with
# 3) id of the peptide B barcode that the read aligned with
# 4) the read sequence
# 5) the read quality string
# 6) mapping quality for pepA (not used for anything currently!)
# 7) mapping quality for pepB (not used for anything currently!)


############### COMMAND LINE ARGUMENTS ###################
my ($file, $fastaDir) = @ARGV;



###############     SET PARAMETERS     ###################
my $buffer = 5;



##############   GET BARCODE INFORMATION   ###############
# my $fastaDir = "./tags/fasta";
# Read tag sequences:
my $pepA   = &ReadBarcodesWithPosFromFasta("$fastaDir/a_epitope_tag.fasta",    0);  # do not rev comp
my $pepB   = &ReadBarcodesWithPosFromFasta("$fastaDir/b_epitope_tag.fasta",    1);  # rev comp each sequence
my $anneal = &ReadBarcodesWithPosFromFasta("$fastaDir/annealing.fasta",    0);  
my $primA  = &ReadBarcodesWithPosFromFasta("$fastaDir/a_forward_primer.fasta", 0);  
# my $primB  = &ReadBarcodesWithPosFromFasta("$fastaDir/primB.fa",     0); 
my $primB  = &ReadBarcodesWithPosFromFasta("$fastaDir/b_forward_primer.fasta",     1); 
# my $Tags   = &ReadBarcodesWithPosFromFasta("$fastaDir/constant.tags.fa", 0); 

### Get tag sequences from fasta file. Store in hash with key = sequence, value = identifier
my $samTags = &ReadBarcodesFromFasta2("$fastaDir/sample_id_tags.fasta", 0);  


############### MAIN PROGRAM ###################

open(IN, "<", $file) or die "Cannot open file $file\n";

# Header
print join("\t", ("read_id", "read_quality", "sample_id", "nw_score", "pepA", "pepB", "sample_seq", "primA_seq", "N6A", "pepA_seq", "anneal_seq", "pepB_seq", "N6B", "primB_seq", "pepA_MAPQ", "pepB_MAPQ")), "\n";

while(defined(my $line = <IN>)) {
	chomp $line;
	# print STDERR "Read line $.\n" if $. =~ m/000$/;
	my ($read_id, $a, $b, $read, $qual, $mapq_a, $mapq_b) = split("\t", $line);
	my ($qual_score) = MeanQuality($qual);

  my %Pos = ();
 
  # # constant sequences:
  $Pos{primA}  = FastFind($read, $primA,  $buffer);
  $Pos{anneal} = FastFind($read, $anneal, $buffer);
  $Pos{primB}  = FastFind($read, $primB,  $buffer);

  # variables sequences:
  $Pos{pepA} = FastFind($read, ${$pepA}{$a}, $buffer);
  $Pos{pepB} = FastFind($read, ${$pepB}{$b}, $buffer);

	
  my $sample = &Cut($read, 0,                    $Pos{primA}{start} - 1); 
  my $fwdpA  = &Cut($read, $Pos{primA}{start},   $Pos{primA}{end});
  my $N1     = &Cut($read, $Pos{primA}{end} + 1, $Pos{pepA}{start} - 1);
  my $codA   = &Cut($read, $Pos{pepA}{start},    $Pos{pepA}{end});
  my $ann    = &Cut($read, $Pos{anneal}{start},  $Pos{anneal}{end});
  my $codB   = &Cut($read, $Pos{pepB}{start},    $Pos{pepB}{end});
  my $N2     = &Cut($read, $Pos{pepB}{end} + 1,  $Pos{primB}{start} - 1);
  my $pB     = &Cut($read, $Pos{primB}{start},   $Pos{primB}{end});

  my ($sample_id, $nw_score);
  # Find sample ID:
  if (exists(${$samTags}{$sample})) {
    # print STDERR "fast $sample\n";
    $sample_id = ${$samTags}{$sample};
    $nw_score = 6;
  } elsif ( scalar(keys %{$samTags}) >= 1000 ) {
    # will be too slow to align if there are many sample ids (eg. >1000)
    $sample_id = "no_match";
    $nw_score = 0;
  } else {
    # print STDERR "slow $sample\n";
    my $max_score = 0;
    my $winner;
    foreach my $s (keys %{$samTags}) {
      my $score = NWalign($sample, $s);
      if ($score > $max_score) {
        $max_score = $score; 
        $winner = ${$samTags}{$s};
      } elsif ($score == $max_score && $score > 0) {
        $winner = "tie";
      }
    }
    $sample_id = $winner;
    $nw_score = $max_score;
  }
  
  print join("\t", ($read_id, $qual_score, $sample_id, $nw_score, $a, $b, $sample, $fwdpA, $N1, $codA, $ann, $codB, $N2, $pB, $mapq_a, $mapq_b)), "\n";


}

close IN;




############### SUB-ROUTINES ###################

sub RevComp {
   my($dna) = @_;
   $dna =~ tr/ATCGatcg/TAGCTAGC/;
   $dna = reverse $dna;
   return $dna;
}


# Like substr, but you can give the start and end positions, instead of start and length.
sub Cut {
  my ($string, $from, $to) = @_;
  my $length = $to - $from + 1;
  return substr($string, $from, $length);
}


# Stores ids as the keys, and the values are a hash with sequence, expected start and expected end position.
# If there is only one entry in the fasta file, the inner hash is returned instead.
sub ReadBarcodesWithPosFromFasta {
	my $file = $_[0];
	my $revcomp = scalar @_ > 1 ? $_[1] : 0;
	my %Hash = ();

	open(IN, "<", $file) or die "Cannot open file $file\n";
	while (defined(my $line = <IN>)) {
		chomp $line;
		my ($id, $start, $end, $seq);
		if ($line =~ m/^>(.*)/) {
			($id, $start, $end) = split(" ", $1);
			if (defined($line = <IN>)) {
				chomp $line;
				$seq = $line;
			} else {
				die "Corrupt FASTA format at line $.\n";
			}
		} else {
			die "Corrupt FASTA format at line $.\n";
		}
		$seq = &RevComp($seq) if $revcomp;
		$Hash{$id} = {seq => $seq, start => $start, end => $end};
	}
  if (scalar (keys %Hash) == 1) {
    my @id = keys %Hash; 
    return($Hash{$id[0]});
  }
  close IN;
	return(\%Hash);
}


# Stores the sequence as the keys, and the ids as the values
sub ReadBarcodesFromFasta2 {
  my $file = $_[0];
  my $revcomp = scalar @_ > 1 ? $_[1] : 0;
  my %Hash = ();

  open(IN, "<", $file) or die "Cannot open file $file\n";
  while (defined(my $line = <IN>)) {
    chomp $line;
    my ($id, $seq);
    if ($line =~ m/^>(\S+)/) {
      $id = $1;
      if (defined($line = <IN>)) {
        chomp $line;
        $seq = $line;
      } else {
        die "Corrupt FASTA format at line $.\n";
      }
    } else {
      die "Corrupt FASTA format at line $.\n";
    }
    $seq = &RevComp($seq) if $revcomp;
    $Hash{$seq} = $id;
  }
  close IN;
  return(\%Hash);
}


# Calculates mean PHRED quality for a read. Also returns the read length, but you can ignore this output if you want.
sub MeanQuality {
  my ($qread) = @_;
  my $base = 33;
  my $n = 0;
  my $sum = 0;
  for (my $i = 0; $i < length $qread; $i++) {
    my $char = substr($qread, $i, 1);
    my $quality = ord($char)-$base;
    $sum += $quality;
    $n++;
  }
  my $mean_qual = sprintf("%.2f", $sum/$n);
  return($mean_qual, $n)
}


sub FastFind {
  my ($seq1, $tag, $buffer) = @_;
  my $len = ${$tag}{end} - ${$tag}{start} + 1; 
  my $hit = -1;
  
  for (my $i = 0; $hit == -1 && $i <= 2; $i++) { # try all values from 0 to $buffer
    
    # try to offset the start, by substracting the value
    my $start = ${$tag}{start} - $i;
    unless ($start < 0) {
      if (${$tag}{seq} eq substr($seq1, $start, $len)) {
        $hit = $start;
      }
    }

    # try to offset the start, by adding the value
    $start = ${$tag}{start} + $i;
    unless ($i == 0 || $start > (length($seq1) - $len)) {
      if (${$tag}{seq} eq substr($seq1, $start, $len)) {
        $hit = $start;
      }
    }
    
  }
  if ($hit >= 0) {
    # print STDERR "${$tag}{start}\t$hit\n";
    # print STDERR "fast\n";
    return ( {start => $hit, end => $hit + $len - 1 });
  }
    # print STDERR "slow\n";
  return ( SWalignTarget($seq1, $tag,  $buffer) );

}


# Trims a read to the region where we expect to find the tag (plus minus buffer size)
# Aligns the tag to the trimmed read
# Adjusts the resulting alignment positions to be relative to full-length read
# Returns the start and end of the alignment in the read, in 0-indexed positions.
sub SWalignTarget {
	my ($seq1, $tag, $buffer) = @_;
    
    # where to cut from:
    my $start = ${$tag}{start} - $buffer;
    # $start--;                   # 1- to 0-indexed
    $start = 0 if $start < 0;

    # where to cut to:
    my $end = ${$tag}{end} + $buffer;
    # $end--;                     # 1- to 0-indexed
    my $l = length $seq1;
    $end = $l if $end > $l-1;

    # truncate read    
    $truncated = substr($seq1, $start, $end-$start+1);
    my $pos = SWalign($truncated, ${$tag}{seq});

    # adjust pos to be to entire read
    ${$pos}[0] += $start;
    ${$pos}[1] += $start;
	
    # output of SW is human-readable positions (1-indexed)
    ${$pos}[0]--;
    ${$pos}[1]--;

    return ( {start => ${$pos}[0], end => ${$pos}[1]} );
    # return $pos;
}



####################################
# Smith-Waterman local alignment
####################################

## Modified from: http://genomics.cribi.unipd.it/~telatin/go/teaching-smith-waterman

sub SWalign {
  my ($seq1, $seq2) = @_;

  # scoring scheme
  my $MATCH      =  3; # +1 for letters that match
  my $MISMATCH   = -3; # -1 for letters that mismatch
  # my $GAP        = -1; # -1 for any gap
  my $GAP_OPEN   = -3; # -1 for a gap opening, including the first gap
  my $GAP_EXTEND = -1; # -1 for a gap extension

  # initialization
  my @matrix;
  $matrix[0][0]{score}   = 0;
  $matrix[0][0]{pointer} = "none";
  $matrix[0][0]{gap} = 0;
  for(my $j = 1; $j <= length($seq1); $j++) {
       $matrix[0][$j]{score}   = 0;
       $matrix[0][$j]{pointer} = "none";
       $matrix[0][$j]{gap} = 0;
  }
  for (my $i = 1; $i <= length($seq2); $i++) {
       $matrix[$i][0]{score}   = 0;
       $matrix[$i][0]{pointer} = "none";
       $matrix[$i][0]{gap} = 0;
  }

  # fill
   my $max_i     = 0;
   my $max_j     = 0;
   my $max_score = 0;

   my @i_start = ();
   my @j_start = ();

   for(my $i = 1; $i <= length($seq2); $i++) {
       for(my $j = 1; $j <= length($seq1); $j++) {
           my ($diagonal_score, $left_score, $up_score);
           
           # calculate match score
           my $letter1 = substr($seq1, $j-1, 1);
           my $letter2 = substr($seq2, $i-1, 1);      
           if ($letter1 eq $letter2) {
               $diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
            }
           else {
               $diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
            }

           # gap status
           $up_GAP   = $matrix[$i-1][$j]{gap} ? $GAP_EXTEND : $GAP_OPEN;
           $left_GAP = $matrix[$i][$j-1]{gap} ? $GAP_EXTEND : $GAP_OPEN;

           
           # calculate gap scores
           $up_score   = $matrix[$i-1][$j]{score} + $up_GAP;
           $left_score = $matrix[$i][$j-1]{score} + $left_GAP;
           
           if ($diagonal_score <= 0 and $up_score <= 0 and $left_score <= 0) {
               $matrix[$i][$j]{score}   = 0;
               $matrix[$i][$j]{pointer} = "none";
               next; # terminate this iteration of the loop
            }

           push(@i_start, $i);
           push(@j_start, $j);

           # choose best score
           if ($diagonal_score >= $up_score) {
               if ($diagonal_score >= $left_score) {
                   $matrix[$i][$j]{score}   = $diagonal_score;
                   $matrix[$i][$j]{pointer} = "diagonal";
                   $matrix[$i][$j]{gap} = 0;
                }
               else {
                   $matrix[$i][$j]{score}   = $left_score;
                   $matrix[$i][$j]{pointer} = "left";
                   $matrix[$i][$j]{gap} = 1;
                }
            } else {
               if ($up_score >= $left_score) {
                   $matrix[$i][$j]{score}   = $up_score;
                   $matrix[$i][$j]{pointer} = "up";
                   $matrix[$i][$j]{gap} = 1;
                }
               else {
                   $matrix[$i][$j]{score}   = $left_score;
                   $matrix[$i][$j]{pointer} = "left";
                   $matrix[$i][$j]{gap} = 1;
                }
            }
           
         # set maximum score
           if ($matrix[$i][$j]{score} > $max_score) {
               $max_i     = $i;
               $max_j     = $j;
               $max_score = $matrix[$i][$j]{score};
            }
        }
   }

   # trace-back

   my $align1 = "";
   my $align2 = "";

   my $j = $max_j;
   my $i = $max_i;

   while (1) {
       last if $matrix[$i][$j]{pointer} eq "none";
       
       if ($matrix[$i][$j]{pointer} eq "diagonal") {
           $align1 .= substr($seq1, $j-1, 1);
           $align2 .= substr($seq2, $i-1, 1);
           $i--; $j--;
        }
       elsif ($matrix[$i][$j]{pointer} eq "left") {
           $align1 .= substr($seq1, $j-1, 1);
           $align2 .= "-";
           $j--;
        }
       elsif ($matrix[$i][$j]{pointer} eq "up") {
           $align1 .= "-";
           $align2 .= substr($seq2, $i-1, 1);
           $i--;
        }  
   }

   $align1 = reverse $align1;
   $align2 = reverse $align2;
   # print "$align1\tfrom $j to $max_j\n";
   # print "$align2\tfrom $i to $max_i\n";
   
   # Added these lines: the index
   $j++;
   $i++;
  return([$j, $max_j, $i, $max_i])
}



#####################################
# Needleman-Wunsch global alignment #
#####################################

## Modified from: http://etutorials.org/Misc/blast/Part+II+Theory/Chapter+3.+Sequence+Alignment/3.1+Global+Alignment+Needleman-Wunsch/

sub NWalign {
  my ($seq1, $seq2) = @_;
  
  # Change both to lower case
  $seq1 = lc($seq1);
  $seq2 = lc($seq2);

  # scoring scheme
  my $MATCH    =  1; # +1 for letters that match
  my $MISMATCH = -1; # -1 for letters that mismatch
  my $GAP      = -1; # -1 for any gap

  # initialization
  my @matrix;
  $matrix[0][0]{score}   = 0;
  $matrix[0][0]{pointer} = "none";
  for(my $j = 1; $j <= length($seq1); $j++) {
      $matrix[0][$j]{score}   = $GAP * $j;
      $matrix[0][$j]{pointer} = "left";
  }
  for (my $i = 1; $i <= length($seq2); $i++) {
      $matrix[$i][0]{score}   = $GAP * $i;
      $matrix[$i][0]{pointer} = "up";
  }

  # fill
  for(my $i = 1; $i <= length($seq2); $i++) {
      for(my $j = 1; $j <= length($seq1); $j++) {
          my ($diagonal_score, $left_score, $up_score);

          # calculate match score
          my $letter1 = substr($seq1, $j-1, 1);
          my $letter2 = substr($seq2, $i-1, 1);                            
          if ($letter1 eq $letter2) {
              $diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
          }
          else {
              $diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
          }

          # calculate gap scores
          $up_score   = $matrix[$i-1][$j]{score} + $GAP;
          $left_score = $matrix[$i][$j-1]{score} + $GAP;

          # choose best score
          if ($diagonal_score >= $up_score) {
              if ($diagonal_score >= $left_score) {
                  $matrix[$i][$j]{score}   = $diagonal_score;
                  $matrix[$i][$j]{pointer} = "diagonal";
              }
          else {
                  $matrix[$i][$j]{score}   = $left_score;
                  $matrix[$i][$j]{pointer} = "left";
              }
          } else {
              if ($up_score >= $left_score) {
                  $matrix[$i][$j]{score}   = $up_score;
                  $matrix[$i][$j]{pointer} = "up";
              }
              else {
                  $matrix[$i][$j]{score}   = $left_score;
                  $matrix[$i][$j]{pointer} = "left";
              }
          }
      }
  }

  # trace-back

  # start at last cell of matrix
  my $j = length($seq1);
  my $i = length($seq2);

  while (1) {
      last if $matrix[$i][$j]{pointer} eq "none"; # ends at first cell of matrix

      if ($matrix[$i][$j]{pointer} eq "diagonal") {
          $i--;
          $j--;
      }
      elsif ($matrix[$i][$j]{pointer} eq "left") {
          $j--;
      }
      elsif ($matrix[$i][$j]{pointer} eq "up") {
          $i--;
      }    
  }
  return $matrix[length($seq2)][length($seq1)]{score};
}