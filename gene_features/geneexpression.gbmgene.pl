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
my ($fIn,$fOut,$expression,$gff);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"i:s"=>\$fIn,
				"gff:s"=>\$gff,
				"e:s"=>\$expression,
				) or &USAGE;
&USAGE unless ($fIn and $fOut and $expression and $gff);

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	zhangyw <yz46606\@uga.edu> 
Program Date:   2017/11/24
Description:	this program is used to gbM gene expression distribution
Usage:
  Options:
  -i <file>  input file,gbM gene list,forced
  -e <file>  input file,expression file,forced
  -gff <file>  input file,gff file,forced
  -o <dir>  output file,forced  
  -h         Help

USAGE
	print $usage;
	exit;
}
my $winsize=50;
my $binnum=20;
mkdir $fOut if (! -d $fOut);

my @ratio=("0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5","0.5-0.6","0.6-0.7","0.7-0.8","0.8-0.9","0.9-1");
#my @ratio=("0.9-0.91","0.91-0.92","0.92-0.93","0.93-0.94","0.94-0.95","0.95-0.96","0.96-0.97","0.97-0.98","0.98-0.99","0.99-1");
#my @ratio=("0.980-0.981","0.981-0.982","0.982-0.983","0.983-0.984","0.984-0.985","0.985-0.986","0.986-0.987","0.987-0.988","0.988-0.989","0.989-0.990","0.990-0.991","0.991-0.992","0.992-0.993","0.993-0.994","0.994-0.995","0.995-0.996","0.996-0.997","0.997-0.998","0.998-0.999","0.999-1");
########################################step 1. divide gene based on their stability
my %genegroup;
&group_gbM_gene(\%genegroup);
#print Dumper %genegroup;die;
########################################step 2. get the gene length for each gbM group
&gene_expression_distribution();

#############################sub function
sub group_gbM_gene{
	my ($group)=@_;
	my $total=927;
	open (IN, $fIn) or die $!;
	my $count;
	while (<IN>) {
		chomp;
		my @lines=split/\s+/,$_;
		$lines[0]=~s/\.\d+//;
		next if (/^$/);
		for (my $i=0;$i<@ratio ;$i++) {
			my ($s,$e)=split/\-/,$ratio[$i];
			my $fre=$lines[1]/$total;
			if ($fre>=$s&&$fre<=$e) {
				$$group{$lines[0]}=$ratio[$i];
			}
		}
		$count++;
	}
	close IN;
}
sub gene_expression_distribution{

	open (IN, $gff) or die $!;
	my %length;
	while (<IN>) {
		chomp;
		next if (/\#/);
		my @lines=split/\s+/,$_;
		if ($lines[2]=~/gene/&&$lines[8]=~/Name=(.*)/) {
			my $name=$1;
			$length{$name}=($lines[4]-$lines[3])/1000;
		}
	}
	close IN;

	open (IN, $expression) or die $!;
	my %expression;
	while (<IN>) {
		chomp;
		next if (/\#/||$.==1);
		my @lines=split/\s+/,$_;
		my $aver;
		for (my $i=1;$i<@lines ;$i++) {
			$aver+=$lines[$i];
		}
		if ($length{$lines[0]}) {
			$aver=int(($aver/(@lines-1))/$length{$lines[0]});
			$expression{$lines[0]}=$aver;
		}
	}
	close IN;

	open (OUT, ">$fOut/summary.txt") or die $!;
	print OUT "gene\tstatus\texpression\n";
	foreach my $gene (keys %expression) {
		if ($genegroup{$gene}) {
			print OUT "$gene\t$genegroup{$gene}\t$expression{$gene}\n";
		}
		else{
			print OUT "$gene\tUM\t$expression{$gene}\n";
		}
	}
	close OUT;
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


