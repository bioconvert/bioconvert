#!/usr/bin/perl -w -- -*-Perl-*-

# This software under THE BEER-WARE LICENSE (Revision 42):
# Sylvain (with the power of octopuses) wrote this file.
# You can do whatever you want with this stuff.
# If we meet some day, and you think this stuff is worth it, you can buy me a beer in return.
use IO::File;
use warnings;
use strict;

my ($fileIn, $fileOut) = @ARGV;
my $ifh = new IO::File "<$fileIn";
my $ofh = new IO::File ">$fileOut";

my $garb = "";

while (my $line = <$ifh>)
{
  print $ofh ">".substr $line, 1, length($line)-1;;
  my $line2 = <$ifh>;
  print $ofh $line2;
  $garb = <$ifh>;
  $garb = <$ifh>;
}

exit ;
