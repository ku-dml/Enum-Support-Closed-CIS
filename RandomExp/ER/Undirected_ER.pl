use strict;
use warnings;

open(DATAFILE, "> output.grh") or die("error :$!");

my $line = <>;
$line = <>;
while( $line = <>){
  chomp($line);
  my @edge = split(/ /, $line);
  print DATAFILE $edge[0]."\t".$edge[1]."\t".$edge[2]."\n";
  print DATAFILE $edge[2]."\t".$edge[1]."\t".$edge[0]."\n";
} 

