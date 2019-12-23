#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Math::CDF qw(:all);
use Statistics::Multtest qw(:all);
my $BEGIN_TIME=time();
my $Time_Start = &sub_format_datetime(localtime($BEGIN_TIME));
print "Program Starts Time:$Time_Start\n";
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"i:s"=>\$fIn,
				) or &USAGE;
&USAGE unless ($fIn and $fOut);

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	zhangyw <yz46606\@uga.edu> 
Program Date:   2016/09/11
Description:	this program is used to convert dmr_original list to standard gabit trait form
Usage:
  Options:
  -i <file>  input file,forced  tsv files
  -o <dir>  output file,forced  
  -h         Help

USAGE
	print $usage;
	exit;
}
##test##########################################
#print Dumper @code;die;
#my @p =(0.01, 0.02, 0.05,0.41,0.16,0.51);
#my $p=\@p;
##print @p;die;
#my $ref=bonferroni($p);
#print $ref->[0];
#die;
#my $cg=1-pbinom(16, 39, 0.24);
#$cg=sprintf("%.3f", $cg);
#print $cg;
#print Dumper @ratio;
#############################################

mkdir $fOut if (! -d $fOut);
my @ratio;
open (IN, $fIn) or die $!;
my @record;
my @rawcg;my @rawchg;my @rawchh;
while (<IN>) {
	chomp;
	next if (/^$/);
	if ($.==1) {
		my @lines=split/\t/,$_;
		my $cg=$lines[2]/$lines[3];
		#print "$lines[2]\t$lines[3]\t$cg";die;
		my $chg=$lines[7]/$lines[8];
		my $chh=$lines[12]/$lines[13];
		push @ratio,$cg;
		push @ratio,$chg;
		push @ratio,$chh;
		next;
	}
	my @lines=split/\s+/,$_;
	push @record,$_;
	my $pcg;my$pchg;my$pchh;

	if ($lines[5]<1) {
		$pcg=1;
	}
	else{
		if ($lines[4]-1>-1) {
			$pcg=1-pbinom($lines[4]-1,$lines[5],$ratio[0]);
		}
		else{
			$pcg=1-pbinom($lines[4],$lines[5],$ratio[0]);
		}
	}
	if ($lines[10]<1) {
		$pchg=1;
	}
	else{
		if ($lines[9]-1>-1) {
			$pchg=1-pbinom($lines[9]-1,$lines[10],$ratio[1]);
		}
		else{
			$pchg=1-pbinom($lines[9],$lines[10],$ratio[1]);
		}
		
	}
	if ($lines[14]<1) {
		$pchh=1
	}
	else{
		if ($lines[15]-1>-1) {
			$pchh=1-pbinom($lines[14]-1,$lines[15],$ratio[2]);
		}
		else {
			$pchh=1-pbinom($lines[14],$lines[15],$ratio[2]);
		}
	}
	
	push @rawcg,$pcg;
	push @rawchg,$pchg;
	push @rawchh,$pchh;
}
close IN;

my $rawcg=\@rawcg;
my $rawchg=\@rawchg;
my $rawchh=\@rawchh;
my $bhcg=BH($rawcg);
my $bhchg=BH($rawchg);
my $bhchh=BH($rawchh);
#print $bhchg->[0];die;

my $base=basename($fIn);
my $name=$base;
my $name2=$base;
$name=~s/1.raw.out/2.pvalue.out/;
$name2=~s/1.raw.out/3.gbm.genes/;

open (OUT, ">$fOut/$name") or die $!;
open (OUT2, ">$fOut/$name2") or die $!;
for (my $i=0;$i<@record ;$i++) {
	my @lines=split/\s+/,$record[$i];
	print OUT "$record[$i]";
	print OUT "\tCG\t$rawcg->[$i]\t$bhcg->[$i]";
	print OUT "\tCHG\t$rawchg->[$i]\t$bhchg->[$i]";
	print OUT "\tCHH\t$rawchh->[$i]\t$bhchh->[$i]\n";
	if ($lines[5]>=20&&$bhcg->[$i]<=0.05&&$bhchg->[$i]>0.05&&$bhchh->[$i]>0.05) {
		print OUT2 "$record[$i]"; 
		print OUT2 "\tCG\t$rawcg->[$i]\t$bhcg->[$i]";
		print OUT2 "\tCHG\t$rawchg->[$i]\t$bhchg->[$i]";
		print OUT2 "\tCHH\t$rawchh->[$i]\t$bhchh->[$i]\n";
	}
}

close OUT;
close OUT2;

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


