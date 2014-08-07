use strict;
use warnings;

if ( not defined $ARGV[1] ) {
  print STDERR "USAGE: $0 path_to/mapped_sequences path_to/genome_mappings\n\n";
  exit;
}

my %copies = ();
open( F, "$ARGV[0]" );
while (<F>) {
  chomp;
  next if ( $_ =~ m/^id/ );
  my @F = split(/\t/);
  next if ( $F[$#F] == 1 );
  next if ( $F[ $#F - 1 ] eq 'bacterial' );
  next if ( $F[ $#F - 1 ] eq 'fungus' );
  next if ( $F[ $#F - 1 ] eq 'vector' );
  next if ( $F[ $#F - 1 ] eq 'viral' );

  # only take unique-mappers
  next if ( $F[11] eq 'NULL' );
  next if ( $F[11] != 1 );
  $copies{ $F[0] } = $F[2];
}
close(F);

open( G, "$ARGV[1]" );
while (<G>) {
  chomp;
  next if ( $_ =~ m/^id/ );
  my @F = split(/\t/);
  next if ( not defined $copies{ $F[1] } );

  # start coordinate in BED is zero-based
  $F[6]--;
  print "$F[2]\t$F[6]\t$F[7]\tseq$F[1]\t$copies{$F[1]}\t$F[5]\n";
}
close(G);


