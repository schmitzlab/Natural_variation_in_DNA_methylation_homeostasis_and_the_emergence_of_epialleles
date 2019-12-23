#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
#use Math::CDF qw(:all);
#use Statistics::Multtest qw(:all);
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $Time_Start = &sub_format_datetime(localtime($BEGIN_TIME));
print "Program Starts Time:$Time_Start\n";
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fa,$gff,$fOut,$tsv);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"fa:s"=>\$fa,
				"gff:s"=>\$gff,
				"tsv:s"=>\$tsv,
				) or &USAGE;
&USAGE unless ($fa and $gff and $fOut and $tsv);

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	zhangyw <yz46606\@uga.edu> 
Program Date:   2016/09/11
Description:	this program is used to find the methylation level of each gene
  Options:
  -fa <file>  input file,fa,forced  
  -gff <file>  input file,gff,forced 
  -tsv <file>  input dir,tsv,forced 
  -o <dir>  output file,forced  
  -h         Help

USAGE
	print $usage;
	exit;
}
mkdir $fOut if (! -d $fOut);
my %context=(
"CGA"=>"CG","CGC"=>"CG","CGG"=>"CG","CGT"=>"CG","CG"=>"CG","CGN"=>"CG","CHG"=>"CHG","CHH"=>"CHH",
"CAG"=>"CHG","CCG"=>"CHG","CTG"=>"CHG","CAA"=>"CHH","CAC"=>"CHH","CAT"=>"CHH",
"CCA"=>"CHH","CCC"=>"CHH","CCT"=>"CHH","CTA"=>"CHH","CTC"=>"CHH","CTT"=>"CHH",
);
my @code=("CG","CHG","CHH");

######################################1. reading gff and protein_primaryTranscriptOnly.fa files to get the cds location of primary transcripts
$/="\>";
open (IN, $fa) or die $!;
my %primaryTrans;
while (<IN>) {
	chomp;
	next if (/^$/);
	my @lines=split/\n/,$_,2;
	my $name=(split/\s+/,$lines[0])[0];
	$primaryTrans{$name}++;
}
close IN;

my %cdsloc;###rearrange the cds order based on strand state
my %bcdsloc;###recording cds in original order of gff files
my %ordergene;
my %gene;
$/="\n";
my %strand;
open (IN, $gff) or die $!;
while (<IN>) {
	chomp;
	next if (/^$/||/^\#/);
	my @lines=split/\s+/,$_;
	if ($lines[2]=~/CDS/&&$lines[0]=~/\d+/) {
		my $gene;
		if ($lines[8]=~/ID\=(.*?)\.Araport11.447.CDS/) {
			$gene=$1;
			if ($primaryTrans{$gene}) {
				if (!$gene{$gene}) {
					push @{$ordergene{$lines[0]}},$gene;
					$gene{$gene}++;
					$strand{$gene}=$lines[6];
				}
				$bcdsloc{$lines[0]}{$gene}.="|$lines[3] $lines[4]";
				#print Dumper %cdsloc;die;
			}
		}
	}
}
close IN;
foreach my $chr (sort {$a<=> $b}keys %bcdsloc) {
	foreach my $gene (@{$ordergene{$chr}}) {
		if (!$bcdsloc{$chr}{$gene}) {
			print "$chr\t$gene\n";
		}
		my @a=split/\|/,$bcdsloc{$chr}{$gene};
		if ($strand{$gene}eq"+") {
			for (my $i=1;$i<@a;$i++) {
				push @{$cdsloc{$chr}{$gene}},$a[$i];
			}
			#print Dumper @{$cdsloc{$chr}{$gene}};die;
		}
		elsif ($strand{$gene}eq"-") {
			#print $bcdsloc{$chr}{$gene};
			for (my $i=@a-1;$i>0;$i--){
				#print "$a[$i]\n";
				push @{$cdsloc{$chr}{$gene}},$a[$i];
			}
			#print Dumper @{$cdsloc{$chr}{$gene}};die;
		}
	}
}
%bcdsloc=();
my $name=basename($tsv);
my %record=();
my %pos=();
open (IN, $tsv) or die $!;
while (<IN>) {
	chomp;
	next if (/^$/||/^\#/||$.==1);
	my @lines=split/\s+/,$_;
	$lines[0]=~s/Chr//;
	if ($lines[5]>=3) {
		$record{$lines[0]}{$lines[1]}=$_;
		push @{$pos{$lines[0]}},$lines[1];
	}
}
close IN;

my %Mstat=();my %Tstat=();
my %MAstat=();my %TAstat=();
my %Mreads=();my %Treads=();
my %MAreads=();my %TAreads=();
foreach my $chr (sort {$a<=> $b}keys %cdsloc) {
	#print $chr;die;
	my $shift=0;
	foreach my $gene (@{$ordergene{$chr}}) {
		MM:foreach my $cds (@{$cdsloc{$chr}{$gene}}) {
			my ($start,$end)=split/\s+/,$cds;
			if ($pos{$chr}&&$shift<@{$pos{$chr}}) {
				for (my $i=$shift;$i<@{$pos{$chr}};$i++) {
					if ($pos{$chr}[$i]>=$start&&$pos{$chr}[$i]<=$end) {
						my @info=split/\s+/,$record{$chr}{$pos{$chr}[$i]};
						if ($context{$info[3]}) {
							$Tstat{$gene}{$context{$info[3]}}++;
							$TAstat{$context{$info[3]}}++;
							$Mreads{$gene}{$context{$info[3]}}+=$info[4];
							$Treads{$gene}{$context{$info[3]}}+=$info[5];
							$MAreads{$context{$info[3]}}+=$info[4];
							$TAreads{$context{$info[3]}}+=$info[5];
							if ($info[-1]==1) {
								$Mstat{$gene}{$context{$info[3]}}++;
								$MAstat{$context{$info[3]}}++;
							}
						}
						
						#print "$shift\t$i\t$pos{$chr}[$i]\t$cds\n";
						$shift=$i+1;
					}
					elsif ($pos{$chr}[$i]>$end) {
						#print "$shift\t$i\t$pos{$chr}[$i]\t$cds\n";
						$shift=$i+1;
						next MM;
					}
					elsif ($pos{$chr}[$i]<$start) {
						#print "$shift\t$i\t$pos{$chr}[$i]\t$cds\n";
						$shift=$i+1;
					}
				}
			}
		}
	}
}
open (OUT, ">$fOut/$name.1.raw.out") or die $!;
print OUT "TOTAL count";
foreach my $cont (@code) {
	if (!$MAstat{$cont}) {
		$MAstat{$cont}="0";
	}
	if (!$TAstat{$cont}) {
		$TAstat{$cont}="0";
	}
	print OUT "\t$cont\t$MAstat{$cont}\t$TAstat{$cont}\t$MAreads{$cont}\t$TAreads{$cont}";
}
print OUT "\n";
foreach my $chr (sort {$a<=> $b}keys %cdsloc) {
	foreach my $gene (@{$ordergene{$chr}}) {
		if ($Tstat{$gene}) {
			print OUT "$chr\t$gene\t$strand{$gene}";
#			foreach my $cds (@{$cdsloc{$chr}{$gene}}) {
#				print OUT "$cds|";
#			}
			foreach my $cont (@code) {
				#print $cont;die;
				if (!$Mstat{$gene}{$cont}) {
					$Mstat{$gene}{$cont}="0";
				}
				if (!$Tstat{$gene}{$cont}) {
					$Tstat{$gene}{$cont}="0";
				}
				if (!$Mreads{$gene}{$cont}) {
					$Mreads{$gene}{$cont}="0";
				}
				if (!$Treads{$gene}{$cont}) {
					$Treads{$gene}{$cont}="0";
				}
				print OUT "\t$cont\t$Mstat{$gene}{$cont}\t$Tstat{$gene}{$cont}\t$Mreads{$gene}{$cont}\t$Treads{$gene}{$cont}";
			}
			print OUT "\n";
		}
	}
}

#######################################################################################
my $Time_End   = sub_format_datetime(localtime(time()));
print STDOUT "Program Ends Time:$Time_End\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
#######################################################################################

sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

#######################################################################################

sub max{#&max(lists or arry);
	#求列表中的最大值
	my $max=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$max=$max>$temp?$max:$temp;
	}
	return $max;
}

#######################################################################################

sub min{#&min(lists or arry);
	#求列表中的最小值
	my $min=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$min=$min<$temp?$min:$temp;
	}
	return $min;
}

#######################################################################################

sub revcom(){#&revcom($ref_seq);
	#获取字符串序列的反向互补序列，以字符串形式返回。ATTCCC->GGGAAT
	my $seq=shift;
	$seq=~tr/ATCGatcg/TAGCtagc/;
	$seq=reverse $seq;
	return uc $seq;			  
}

#######################################################################################

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

#######################################################################################

sub sub_format_datetime {#Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


