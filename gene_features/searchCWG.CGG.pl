#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut,$fa,$gff);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"i:s"=>\$fIn,
				"g:s"=>\$gff,
				"f:s"=>\$fa,
				) or &USAGE;
&USAGE unless ($fIn and $fOut and $fa and $gff);

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	zhangyw <yz46606\@uga.edu> 
Program Date:   2016/09/11
Description:	this program is used to convert dmr_original list to standard gabit trait form
Usage:
  Options:
  -i <file>  input file,forced  
  -g <file>  input file,gff,forced  
  -f <file>  input file,fa,gene sequence,forced  
  -o <file>  output file,forced  
  -h         Help

USAGE
	print $usage;
	exit;
}
my %context=("CAG"=>"CWG","CTG"=>"CWG","CGG"=>"CCG","CCG"=>"CCG");
my %geneid;
my %len;
open (IN, $gff) or die $!;
while (<IN>) {
	chomp;
	next if (/^$/||/^\#/);
	my @lines=split/\s+/,$_;
	if ($lines[2]=~/gene/&&$lines[8]=~/ID\=(.*?)\.g/) {
		my $len=abs($lines[4]-$lines[3]);
		$geneid{$1}="$lines[0]\t$lines[3]\t$lines[4]\t$lines[6]\t$len";
		$len{$1}=$len;
	}
}
close IN;

$/="\>";
my %CHG;
my %total;
open (IN, $fa) or die $!;
while (<IN>) {
	chomp;
	next if (/^$/);
	my @lines=split/\n/,$_,2;
	my @head=split/\s+/,$lines[0];
	#print $head[0];die;
	if ($geneid{$head[0]}) {
		$lines[1]=~s/\n//;
		my @base=split//,$lines[1];
		for (my $i=0;$i<@base-2 ;$i++) {
			my $txt="$base[$i]$base[$i+1]$base[$i+2]";
			if ($context{$txt}) {
				$CHG{$head[0]}{$context{$txt}}++;
				$total{$head[0]}++;
			}
		}
		#print Dumper %CHG;die;
	}
}
close IN;

$/="\n";
my %hchg;
open (IN, $fIn) or die $!;
while (<IN>) {
	chomp;
	next if (/^$/);
	my @lines=split/\s+/,$_;
	$hchg{$lines[0]}++;
}
close IN;

open (OUT, ">$fOut") or die $!;
print OUT "Gene\tstate\tscaff\tst\tend\tstrand\tlength\tCWGn\tCCGn\tCHGn\tr1\tr2\tr3\n";
my $state;
foreach my $gene (keys %geneid) {
	if ($hchg{$gene}) {
		$state="hchg";
	}
	else{$state="not"}
	if (!$CHG{$gene}{"CWG"}) {
		$CHG{$gene}{"CWG"}=0.01;
	}
	if (!$CHG{$gene}{"CCG"}) {
		$CHG{$gene}{"CCG"}=0.01;
	}
	if (!$total{$gene}) {
		$total{$gene}=0.01;
	}
	my $r1=$CHG{$gene}{"CWG"}/$len{$gene};
	my $r2=$CHG{$gene}{"CCG"}/$len{$gene};
	my $r3=$total{$gene}/$len{$gene};
	print OUT "$gene\t$state\t$geneid{$gene}\t$CHG{$gene}{\"CWG\"}\t$CHG{$gene}{\"CCG\"}\t$total{$gene}\t$r1\t$r2\t$r3\n";
}
close OUT;
#####################################################

