use List::Util qw(sum); # Enables make average in numeric arrays.

#Basic open file variables. #use of shift command to push data one by one to variable
 $fasta_file=shift;
 $fh;
open($fh, $fasta_file) or die "can't open $fasta_file: $!\n";

 %sequence_data;

	### array to define variables enabling printing everything in columns, for sub routine 1 to X.
@sequence_header = @totalAA = @positive_charged = @negative_charged = @polar_uncharged = ();
@special_class = @hydrofobic = @pH = @totalmw = @gravy = @AliphaticIndex = @amountC = ();
@totalH = @amountN = @totalO = @amountS = @totalAtoms = @countR = @statR = @countH = ();
@statH = @countK = @statK = @countD = @statD = @countE = @statE = @countS = @statS = ();
@countT = @statT = @countN = @statN = @countQ = @statQ = @countC = @statC = @countG= ();
@statG = @countP = @statP = @countA = @statA = @countV = @statV = @countI = ();
@statI = @countL = @statL = @countM = @statM = @countF = @statF = @countY= ();
@statY = @countW = @statW = @pH = ();
	
	### Algorithm for all calculations from file using a single while loop.
while (read_fasta_sequence($fh, \%sequence_data)) {

	## Expasy values to variables, count and storage of variables
	$countR = $sequence_data{seq} =~ tr/[R]//;
	$countH = $sequence_data{seq} =~ tr/[H]//;
	$countK = $sequence_data{seq} =~ tr/[K]//;
	$countD = $sequence_data{seq} =~ tr/[D]//;
	$countE = $sequence_data{seq} =~ tr/[E]//;
	$countS = $sequence_data{seq} =~ tr/[S]//;
	$countT = $sequence_data{seq} =~ tr/[T]//;
	$countN = $sequence_data{seq} =~ tr/[N]//;
	$countQ = $sequence_data{seq} =~ tr/[Q]//;
	$countC = $sequence_data{seq} =~ tr/[C]//;
	$countG = $sequence_data{seq} =~ tr/[G]//;
	$countP = $sequence_data{seq} =~ tr/[P]//;
	$countA = $sequence_data{seq} =~ tr/[A]//;
	$countV = $sequence_data{seq} =~ tr/[V]//;
	$countI = $sequence_data{seq} =~ tr/[I]//;
	$countL = $sequence_data{seq} =~ tr/[L]//;
	$countM = $sequence_data{seq} =~ tr/[M]//;
	$countF = $sequence_data{seq} =~ tr/[F]//;
	$countY = $sequence_data{seq} =~ tr/[Y]//;
	$countW = $sequence_data{seq} =~ tr/[W]//;
#print "Header:>$sequence_data{header}\nwith Sequence:$sequence_data{seq}\n"; 
#printing of the header and sequence
		push (@sequence_header, $sequence_data{header});
		push (@sequence_residues, $sequence_data{seq});

	#specific amino acids groups in percentage ARRAYS
	$totalaa = $sequence_data{seq}=~ tr/[RHKDESTNQCGPAVILMFYW]//;
		push (@totalAA, $totalaa);

	$count_of_positive_charged = $sequence_data{seq} =~ tr/[RHK]//;
	$pos_char=(($count_of_positive_charged/$totalaa)*100);
		push (@positive_charged, $pos_char);

	$count_of_negative_charged = $sequence_data{seq}=~ tr/[DE]//;
	$neg_char=(($count_of_negative_charged/$totalaa)*100);
		push (@negative_charged, $neg_char);

	$count_of_polar_uncharged = $sequence_data{seq}=~ tr/[STNQ]//;
	$per_polar=(($count_of_polar_uncharged/$totalaa)*100);
		push (@polar_uncharged, $per_polar);

	$count_of_special_class = $sequence_data{seq}=~ tr/[CGP]//;
	$per_spe=(($count_of_special_class/$totalaa)*100);
		push (@special_class, $per_spe);

	$count_of_hydrophobic = $sequence_data{seq}=~ tr/[AVILMFYW]//;
	$per_hydro=(($count_of_hydrophobic/$totalaa)*100);
		push (@hydrofobic, $per_hydro);

		#####isoelectric point, naive algorithm.
	for ($pH = 0.00; $pH < 14.00; $pH+=0.1) {
#		$pH=0.0;
#		while ($pH<=14) {
#		print "\n";
	$chargedAAcount = ($countR + $countH + $countK + $countD + $countE + $countC + $countY);
	$QN1=(-1/(1+(10**(3.65-$pH))));		#C-terminal charge
	$QN2=(-$countD/(1+(10**(3.9-$pH))));	#D charge
	$QN3=(-$countE/(1+(10**(4.07-$pH))));	#E charge
	$QN4=(-$countC/(1+(10**(8.18-$pH))));	#C charge
	$QN5=(-$countY/(1+(10**(10.46-$pH))));	#Y charge
	$QP1=($countH/(1+(10**($pH-6.04))));	#H charge
	$QP2=(1/(1+(10**($pH-8.2))));		#NH2charge
	$QP3=($countK/(1+(10**($pH-10.54))));	#K charge
	$QP4=($countR/(1+(10**($pH-12.48))));	#R charge
	$NQ=($QN1+$QN2+$QN3+$QN4+$QN5+$QP1+$QP2+$QP3+$QP4);
	$NQ = sprintf "%.7f", $NQ;
	$pH = sprintf "%.1f", $pH;			
			#PI is found at NQ = 0
#	print "pH: $pH\n";
#	print "NQ: $NQ\n";
#			if (($NQ == 0) || ($NQ < 0.01) && ($NQ >! 0.01)){			
			if ($NQ < 1) { # not working
				#print "$pH\n";
				push (@pH, $pH);				
				last;			}
				#print "\n";
				
				#@pH=$pH;	
				#print "@pHt[0]\n";	
#				$pH = sprintf "%.2f", $pH;
#					print "protein PI:".$pH."\n";
#					push (@pH, $pH);				
				#push (@pH, $pHt[0]);
#				print "@pH\n";

			} 
		##### end of isoelectric point algorithm
				
		#molecular weight algorithm
$totalmw=(($countR*(174.2)+$countH*(155.15)+$countK*(146.19)+$countD*(133.10)+$countE*(147.13)+$countS*(105.09)+$countT*(119.12)+$countN*(132.12)+$countQ*(146.14)+$countC*(121.16)+$countG*(75.07)+$countP*(115.13)+$countA*(89.09)+$countV*(117.15)+$countI*(131.17)+$countL*(131.17)+$countM*(149.21)+$countF*(165.19)+$countY*(181.19)+$countW*(204.23))-($totalaa*(18.01)));
			$totalmw = sprintf "%.3f", $totalmw; #transform to 3 decimals
#			print "protein molecular weight is:".$totalmw."\n";
			push (@totalmw, $totalmw);

		#Hydrophobicity index algorithm
$gravy = ((0.05 + ($countR*(-4.500) + $countH*(-3.200) + $countK*(-3.900) + $countD*(-3.500) + $countE*(-3.500) + $countS*(-0.800) + $countT*(-0.700) + $countN*(-4.500) + $countQ*(-3.500) + $countC*(2.500) + $countG*(-0.400) + $countP*(-1.600) + $countA*(1.800) + $countV*(4.200) + $countI*(4.500) + $countL*(3.800) + $countM*(1.900) + $countF*(2.800) + $countY*(-1.300) + $countW*(-0.900)))/($totalaa + 0.001));
			$gravy = sprintf "%.3f", $gravy; # 3 decimals transformation
#			print "GRAVY index:".$gravy."\n";
			push (@gravy, $gravy);


	## atomic composition and count
	$amountC = ($countR*(6) + $countH*(6) + $countK*(6) + $countD*(4) + $countE*(5) + $countS*(3) + $countT*(4) + $countN*(4) + $countQ*(5) + $countC*(3) + $countG*(2) + 		$countP*(5) + $countA*(3) + $countV*(5) + $countI*(6) + $countL*(6) + $countM*(5) + $countF*(9) + $countY*(9) + $countW*(11));
	
	$amountH = ($countR*(14) + $countH*(9) + $countK*(14) + $countD*(7) + $countE*(9) + $countS*(7) + $countT*(9) + $countN*(8) + $countQ*(10) + $countC*(7) + $countG*(5) 		+ $countP*(9) + $countA*(7) + $countV*(11) + $countI*(13) + $countL*(14) + $countM*(11) + $countF*(11) + $countY*(11) + $countW*(12));
	
	$amountN = ($countR*(4) + $countH*(3) + $countK*(2) + $countD*(1) + $countE*(1) + $countS*(1) + $countT*(1) + $countN*(2) + $countQ*(2) + $countC*(1) + $countG*(1) + 		$countP*(1) + $countA*(1) + $countV*(1) + $countI*(1) + $countL*(1) + $countM*(1) + $countF*(1) + $countY*(1) + $countW*(2));
	
	$amountO = ($countR*(2) + $countH*(2) + $countK*(2) + $countD*(4) + $countE*(4) + $countS*(3) + 
	$countT*(3) + $countN*(3) + $countQ*(3) + $countC*(2) + $countG*(2) + $countP*(2) + $countA*(2) + $countV*(2) + $countI*(2) + $countL*(2) + $countM*(2) + $countF*(2) 		+ $countY*(3) + $countW*(2));

			$amountS = ($countM+$countC);
			$totalH = ($amountH-(($totalaa-1)*2));
			$totalO = ($amountO-($totalaa-1));
#				print "Atomic Composition: C$amountC";
				#if (S > 0) { #loop to display S only for aa with S
					#print "S$amountS";
					#	}	#and of loop

				#to print number ot atoms
		#		print "\n";
				$totalAtoms = ($totalH + $amountN + $totalO + $amountS + $amountC);
		#		print "Total atoms: $totalAtoms\n";
				push (@totalAtoms, $totalAtoms);

					$perC=(($amountC/$totalAtoms)*100);
					$perH=(($totalH/$totalAtoms)*100);
					$perN=(($amountN/$totalAtoms)*100);
					$perO=(($totalO/$totalAtoms)*100);
					$perS=(($amountS/$totalAtoms)*100);

					push (@amountC, $perC);
	#				print "H$totalH";	
					push (@totalH, $perH);
	#				print "N$amountN";
					push (@amountN, $perN);
	#				print "O$totalO";
					push (@totalO, $perO);
	#				print "S$amountS";
					push (@amountS, $perS);		

		
	## end of atomic composition and count

		#aliphatic index algorithm
		$molarFrac = ($countR + $countH + $countK + $countD + $countE + $countS + $countT + $countN + $countQ + $countC + $countG + $countP + $countA + $countV + 			$countI + $countL + $countM + $countF + $countY + $countW);

		$molarA = ( ($countA)  /($molarFrac+0.0015));
		$molarV = ( ($countV)  /($molarFrac+0.0015));
		$molarI = ( ($countI)  /($molarFrac+0.0015));
		$molarL = ( ($countL)  /($molarFrac+0.0015));
				$AliphaticIndex = ( (100*($molarA))+2.9*(100*($molarV))+3.9*( (100*($molarI))+(100*($molarL)) ) );
				$AliphaticIndex = sprintf "%.3f", $AliphaticIndex; # 3 decimals 
		#		print "Aliphatic index:".$AliphaticIndex."\n";
				push (@AliphaticIndex, $AliphaticIndex);

#calculations of statistics and printing
#		print "Amino acid statistics of [>$sequence_data{header}]:\n";
$statR = (($countR/($totalaa))*100);
$statR = sprintf "%.2f", $statR;
#print "R($countR):$statR%\n";
push (@countR, $countR);
push (@statR, $statR);
$statH = (($countH /$totalaa)*100);
$statH = sprintf "%.2f", $statH;
#print "H($countH):$statH%\n";
push (@countH, $countH);
push (@statH, $statH);
$statK = (($countK/$totalaa)*100);
$statK = sprintf "%.2f", $statK;
#print "K($countK):$statK%\n";
push (@countK, $countK);
push (@statK, $statK);
$statD = (($countD/$totalaa)*100);
$statD = sprintf "%.2f", $statD;
#print "D($countD):$statD%\n";
push (@countD, $countD);
push (@statD, $statD);
$statE = (($countE/$totalaa)*100);
$statE = sprintf "%.2f", $statE;
##print "E($countE):$statE%\n";
push (@countE, $countE);
push (@statE, $statE);
$statS = (($countS/$totalaa)*100);
$statS = sprintf "%.2f", $statS;
##print "S($countS):$statS%\n";
push (@countS, $countS);
push (@statS, $statS);
$statT = (($countT/$totalaa)*100);
$statT = sprintf "%.2f", $statT;
#print "T($countT):$statT%\n";
push (@countT, $countT);
push (@statT, $statT);
$statN = (($countN/$totalaa)*100);
$statN = sprintf "%.2f", $statN;
#print "N($countN):$statN%\n";
push (@countN, $countN);
push (@statN, $statN);
$statQ = (($countQ/$totalaa)*100);
$statQ = sprintf "%.2f", $statQ;
#print "Q($countQ):$statQ%\n";
push (@countQ, $countQ);
push (@statQ, $statQ);
$statC = (($countC/$totalaa)*100);
$statC = sprintf "%.2f", $statC;
#print "C($countC):$statC%\n";
push (@countC, $countC);
push (@statC, $statC);
$statG = (($countG/$totalaa)*100);
$statG = sprintf "%.2f", $statG;
#print "G($countG):$statG%\n";
push (@countG, $countG);
push (@statG, $statG);
$statP = (($countP/$totalaa)*100);
$statP = sprintf "%.2f", $statP;
#print "P($countP):$statP%\n";
push (@countP, $countP);
push (@statP, $statP);
$statA = (($countA/$totalaa)*100);
$statA = sprintf "%.2f", $statA;
#print "A($countA):$statA%\n";
push (@countA, $countA);
push (@statA, $statA);
$statV = (($countV/$totalaa)*100);
$statV = sprintf "%.2f", $statV;
#print "V($countV):$statV%\n";
push (@countV, $countV);
push (@statV, $statV);
$statI = (($countI/$totalaa)*100);
$statI = sprintf "%.2f", $statI;
#print "I($countI):$statI%\n";
push (@countI, $countI);
push (@statI, $statI);
$statL = (($countL/$totalaa)*100);
$statL = sprintf "%.2f", $statL;
#print "L($countL):$statL%\n";
push (@countL, $countL);
push (@statL, $statL);
$statM = (($countM/$totalaa)*100);
$statM = sprintf "%.2f", $statM;
#print "M($countM):$statM%\n";
push (@countM, $countM);
push (@statM, $statM);
$statF = (($countF/$totalaa)*100);
$statF = sprintf "%.2f", $statF;
#print "F($countF):$statF%\n";
push (@countF, $countF);
push (@statF, $statF);
$statY = (($countY/$totalaa)*100);
$statY = sprintf "%.2f", $statY;
#print "Y($countY):$statY%\n";
push (@countY, $countY);
push (@statY, $statY);
$statW = (($countW/$totalaa)*100);
$statW = sprintf "%.2f", $statW;
#print "W($countW):$statW%\n\n";
push (@countW, $countW);
push (@statW, $statW);
}

#loop to read each fasta sequence individualy
sub read_fasta_sequence {
    ($fh, $seq_info) = @_;

   $seq_info->{seq} = undef; # clear out previous sequence

   # put the header into place
   $seq_info->{header} = $seq_info->{next_header} if $seq_info->{next_header};

    $file_not_empty = 0; 
   while (<$fh>) {
      $file_not_empty = 1;
      next if /^\s*$/;  # skip blank lines
      chomp;    

      if (/^>/) { # fasta header line
          $h = $_;    
         $h =~ s/^>//;  
         if ($seq_info->{header}) {
            $seq_info->{next_header} = $h;
            return $seq_info;   
         }              
         else { # first time through only
            $seq_info->{header} = $h;
         }              
      }         
      else {    
         s/\s+//;  # remove any white space
         $seq_info->{seq} .= $_;
      }         
   }    

   if ($file_not_empty) {
      return $seq_info;
   }    
   else {
      # clean everything up
      $seq_info->{header} = $seq_info->{seq} = $seq_info->{next_header} = undef;
      return;   
   }    
}

sub zup {
  join "\n", map { my $i = $_; join ' ', map $_->[ $i ], @_ } 0 .. $#{
+ $_[0] }
}


### printing the arrays

print "Header Total.amino.acids Positive.charged.RHK Negative.charged.DE Uncharged.STNQ Special.CGP HydrophobicAVILMFYW Molecular.Weight gravy Aliphatic.index Percentage.C Percentage.H Percentage.N Percentage.O Percentage.S Stat.R Stat.H Stat.K Stat.D Stat.E Stat.S Stat.T Stat.N Stat.Q Stat.C Stat.G Stat.P Stat.A Stat.V Stat.I Stat.L Stat.M Stat.F Stat.Y Stat.W Sequence Isoelectric.point";
print "\n";
print zup \(@sequence_header, @totalAA, @positive_charged, @negative_charged, @polar_uncharged, @special_class, @hydrofobic, @totalmw, @gravy, @AliphaticIndex, @amountC, @totalH, @amountN, @totalO, @amountS, @statR, @statH, @statK, @statD, @statE, @statS, @statT, @statN, @statQ, @statC, @statG, @statP, @statA, @statV, @statI, @statL, @statM, @statF, @statY, @statW, @sequence_residues, @pH);
#print "\n";
#print " ";
#foreach (@pH) {  print $_; print "\n";}



