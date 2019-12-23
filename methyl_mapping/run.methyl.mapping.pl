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
my ($fIn,$fOut);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"i:s"=>\$fIn,
				) or &USAGE;
&USAGE unless ($fIn and $fOut);

sub USAGE {##
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	zhangyw <yz46606\@uga.edu> 
Program Date:   2018/03/12
Description:	this program is used to do the methylome mapping
Usage:
  Options:
  -i <file>  input file,forced  
  -o <dir>  output file,forced  
  -h         Help

USAGE
	print $usage;
	exit;
}
mkdir $fOut if (! -d $fOut);
$fIn=&ABSOLUTE_DIR($fIn);
$fOut=&ABSOLUTE_DIR($fOut);

my $forward="/lustre1/yz46606/02_reference/arab/TAIR10_chr_all.Chr_f";
my $reverse="/lustre1/yz46606/02_reference/arab/TAIR10_chr_all.Chr_r";
my $ref="/lustre1/yz46606/02_reference/arab/TAIR10_chr_all.Chr.fa";
my @file1=glob"$fIn/*.fastq";
my $count=1;
for (my $i=0;$i<@file1 ;$i++) {
	my $read=$file1[$i];
	my $name=basename($file1[$i]);
	$name=~s/\.fastq//;

	open (OUT, ">$fOut/methymapping.$count.sh") or die $!;
	print OUT "#!/bin/bash\n";
	print OUT "#PBS -N methyl_$count\n";
	print OUT "#PBS -l walltime=480:00:00\n";
	print OUT "#PBS -l nodes=1:ppn=6:AMD\n";
	print OUT "#PBS -l mem=20g\n";
	print OUT "#PBS -q batch\n";
	print OUT "\n";
	print OUT "cd \$PBS_O_WORKDIR\n";
	print OUT "module load methylpy\n";
	print OUT "module load picard/2.16.0-Java-1.8.0_144\n";
	print OUT "methylpy single-end-pipeline  --read-files $read --sample $name --forward-ref $forward --reverse-ref $reverse --ref-fasta  $ref --binom-test True --min-cov 3   --sort-mem 20G  --unmethylated-control ChrC --num-procs 3  --remove-clonal true --compress-output false --path-to-picard=\"/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144\"\n";
	$count++;

}
close OUT;
my @sh=glob"$fOut/*.sh";
foreach my $sh (@sh) {
	`qsub $sh`;
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


