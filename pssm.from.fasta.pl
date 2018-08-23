#Basic open file variables. #use of shift command to push data one by one to variable
 $fasta_file=shift;
 $fh;
open($fh, $fasta_file) or die "can't open $fasta_file: $!\n";

 %sequence_data;

while (read_fasta_sequence($fh, \%sequence_data)) {


# start of SETUP FOR PSSM for each sequence
%pwm;
    ++$pwm{ substr $sequence_data{seq}, $_, 1 }[ $_ ] for 0 .. length( $sequence_data{seq} ) -1;
$n = $.;
@$_ = map{ $_ ? $_ / $n : 0 } @$_ for values %pwm;
# end of SETUP FOR PSSM for each sequence
}

### fasta reader
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
############ end fasta reader
sub zup {
  join "\n", map { my $i = $_; join ' ', map $_->[ $i ], @_ } 0 .. $#{
+ $_[0] }
}

print "amino.acid position.one position.two position.three \n";
#@data = sprintf "%9.f", @{ $pwm{ $_ } };
#print "$_ @data } \n" for keys %pwm;
print "$_ @{ $pwm{ $_ } } \n" for keys %pwm;





