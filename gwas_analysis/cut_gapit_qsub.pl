#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
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
my ($fIn,$fOUT,$split,$snp);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOUT,
				"i:s"=>\$fIn,
				"snp:s"=>\$snp,
				"split:i"=>\$split,
				) or &USAGE;
&USAGE unless ($fIn and $fOUT and $snp);

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
  -snp <file>  snp file,forced
  -split <num>  the number of files,default 50
  -o <dir>  OUTput file,forced  
  -h         Help

USAGE
	print $usage;
	exit;
}

mkdir $fOUT if (! -d $fOUT);
$fIn=&ABSOLUTE_DIR($fIn);
$fOUT=&ABSOLUTE_DIR($fOUT);
$snp=&ABSOLUTE_DIR($snp);


my $tfile="$fOUT/trait"; mkdir $tfile if (! -d $tfile);
my $result="$fOUT/result";mkdir $result if (! -d $result);

$split||=50;
my $num;
my @tfile;
open (IN, $fIn) or die $!;
my $count=0;
while (<IN>) {
	chomp;
	my @lines=split/\s+/,$_;
	#print $#lines;
	next if (/^$/);
	if ($.==1) {
		$num=int((@lines-1)/$split);
		#print $num;die;
	}	
	my $file=0;
	my $start=1;
	my $flag=1;
	for (my $i=1;$i<@lines ;$i++) {
		#print $lines[$i];die;
		if ($flag<=$num) {
			if ($start==1) {
				push @{$tfile[$file][$count]},$lines[0];
				push @{$tfile[$file][$count]},$lines[$i];
			}
			else {push @{$tfile[$file][$count]},$lines[$i]};
			$start++;
		}
		else {
			$file++;
			$flag=1;
			$start=1;
			push @{$tfile[$file][$count]},$lines[0];
			push @{$tfile[$file][$count]},$lines[$i];
			$start++;
		}
		$flag++;
	}
	$count++;
}
close IN;
#print Dumper @{$tfile[0][0]};die;
for (my $i=0;$i<@tfile ;$i++) {
	open (OUT, ">$tfile/trait.$i.txt") or die $!;
	for (my $j=0;$j<=$#{$tfile[$i]} ;$j++) {
		for (my $k=0;$k<=$#{$tfile[$i][$j]} ;$k++) {
			#print $tfile[$i][$j][$k];die;
			if ($k+1==1) {
				print OUT "$tfile[$i][$j][$k]";
			}
			else {
				print OUT "\t$tfile[$i][$j][$k]";
			}
		}
		print OUT "\n";
	}
	close OUT;
}

my @trait=glob"$tfile/*.txt";
open SH,">$result/qsub.sh";
foreach my $t (@trait) {
	my $name=basename($t);
	my $dir="$result/$name";
	mkdir $dir if (! -d $dir);
	open OUT,">$dir/gabit.R";
	print OUT "library('MASS')\n";
	print OUT "library(multtest)\n";
	print OUT "library(gplots)\n";
	print OUT "library(compiler)\n";#required for cmpfun
	print OUT "library(\"scatterplot3d\")\n";
	print OUT "library(ape)\n";
	print OUT "source(\"http://www.zzlab.net/GAPIT/emma.txt\")\n";
	print OUT "source(\"http://www.zzlab.net/GAPIT/gapit_functions.txt\")\n";
	print OUT "myY <- read.table(\"$t\", head = TRUE)\n";
	print OUT "myG <- read.table(\"$snp\", head = FALSE)\n";
	print OUT "setwd(\"$dir\")\n";
	print OUT "myGAPIT <- GAPIT(Y=myY,G=myG,PCA.total=3,kinship.cluster=c(\"average\", \"complete\", \"ward\"),kinship.group=c(\"Mean\", \"Max\"),group.from=200,group.to=1000000,group.by=10)";
	close OUT;
	print SH "R --slave < $dir/gabit.R\n";
}
close SH;
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

sub sub_format_datetime {#Time calculation subrOUTine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


