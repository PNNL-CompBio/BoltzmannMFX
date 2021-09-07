#!/usr/bin/env perl
$t = 0.0;
while(<STDIN>) {
  if (/old_time.*new time\s+(\S+)/) {
    $t = $1;
  }
  if (/particle: 0 position:\s+\S+\s+\S+\s+(\S+)/) {
    $z = $1-0.001848;
    print "$t $z\n";
  }
}
