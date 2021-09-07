#!/usr/bin/env perl
$t = 0.0;
$y1 = 0.0;
$y2 = 0.0;
while(<STDIN>) {
  if (/old_time.*new time\s+(\S+)/) {
    $t = $1;
  }
  if (/particle: 0 position:\s+\S+\s+(\S+)/) {
    $y1 = $1;
  }
  if (/particle: 1 position:\s+\S+\s+(\S+)/) {
    $y2 = $1;
    $r = abs($y2 - $y1);
    print "$t $r\n";
  }
}
