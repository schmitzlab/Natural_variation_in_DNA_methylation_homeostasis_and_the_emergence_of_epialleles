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
my ($fIn,$fOut,$gff);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"i:s"=>\$fIn,
				"gff:s"=>\$gff,
				) or &USAGE;
&USAGE unless ($fIn and $fOut and $gff);

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	zhangyw <yz46606\@uga.edu> 
Program Date:   2016/09/11
Description:	this program is used to add gene feature of exon number and transcripts number to genes 
Usage:
  Options:
  -i <file>  input file,gbM LIST,forced  
  -gff <file>  input file,gff file,forced  
  -o <file>  output file,forced  
  -h         Help

USAGE
	print $usage;
	exit;
}

#mkdir $fOut if (! -d $fOut);
#####pre defined hash
my %chrlen=(
"1"=>"30427671","2"=>"19698289","3"=>"23459830","4"=>"18585056","5"=>"26975502",
);
my %context=(
"CGA"=>"CG","CGC"=>"CG","CGG"=>"CG","CGT"=>"CG","CG"=>"CG","CGN"=>"CG","CHG"=>"CHG","CHH"=>"CHH",
"CAG"=>"CHG","CCG"=>"CHG","CTG"=>"CHG","CAA"=>"CHH","CAC"=>"CHH","CAT"=>"CHH",
"CCA"=>"CHH","CCC"=>"CHH","CCT"=>"CHH","CTA"=>"CHH","CTC"=>"CHH","CTT"=>"CHH",
);
my @code=("CG","CHG","CHH");
#$/="\n";

## get the location of gbm genes


my %exonnum;
my %transcripts;
&read_gff(\%exonnum,\%transcripts);

&output(\%exonnum,\%transcripts);
sub output {
	open (IN, $fIn) or die $!;
	open (OUT, ">$fOut") or die $!;
	print OUT "Gene\tclass\tchr\tstart\tend\texpression\tgenelen\taltnum\texonum\n";
	while (<IN>) {
		chomp;
		my @lines=split/\s+/,$_;
		next if (/^$/||/^\#/);
		if (!$transcripts{$lines[0]}||!$exonnum{$lines[0]}) {
			print "$_\n"
		}
		else{
			print OUT "$_\t$transcripts{$lines[0]}\t$exonnum{$lines[0]}\n";
		}
	}
	close OUT;
}
############################################################# SUB


sub read_gff {

	open (IN, $gff) or die $!;
	while (<IN>) {
		chomp;
		my @lines=split/\s+/,$_;
		next if (/^$/||/^\#/);
		if ($lines[0]=~/^(\d+)/) {
			my $chr=$1;
			if ($lines[2]=~/mRNA/&&$lines[8]=~/ID\=(.*?)\.(\d+)\.Araport11\.447/) {
				$transcripts{$1}++;
			}
			if ($lines[2]=~/CDS/&&$lines[8]=~/ID\=(.*?)\.1\.Araport11\.447/) {
				$exonnum{$1}++;
			}
		}
	}
}
