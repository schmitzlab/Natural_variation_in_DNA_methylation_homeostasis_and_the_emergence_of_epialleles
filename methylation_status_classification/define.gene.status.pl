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
my ($fIn,$fOut,$fIn2);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"i:s"=>\$fIn,
				"n:s"=>\$fIn2,
				) or &USAGE;
&USAGE unless ($fIn and $fIn2 and $fOut);

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	zhangyw <yz46606\@uga.edu> 
Program Date:   2018/07/3
Description:	this program is used to generate a matrix for all accession to show their gbM conservation
Usage:
  Options:
  -i <dir>  input dir,each accession's methylation info forced  
  -n <file>  input file,number of gbM for accessions forced   
  -o <dir>  output file,forced  
  -h         Help

USAGE
	print $usage;
	exit;
}
mkdir $fOut if (! -d $fOut);
my @code=("CG","CHG","CHH");
#############STEP 1 define methylation status ##cg-gene=1,chg gene=2 chh gene=3,limitedcoverage=4(chh reads count<20 
#############since chh sites number is the most loose condition among 3 context),unmethylated=5(left genes)
my %genestatus;
my %genemlevel;
my @sam;
&define_status($fIn,\%genestatus,\%genemlevel,\@sam);


###############step 2 output matrix
&raw_output(\%genestatus,\%genemlevel);

##############################sub scripts
sub raw_output {
	my @order;
	open (IN, $fIn2) or die $!;
	my %gbMgene;
	while (<IN>) {
		chomp;
		my @lines=split/\s+/,$_;
		$lines[0]=~s/\.\d+//;
		$gbMgene{$lines[0]}=$lines[1];
	}
	close IN;
	foreach my $gbMgene (sort{$gbMgene{$b}<=>$gbMgene{$a}}keys %gbMgene) {
		push @order,$gbMgene;
	}
	foreach my $gene (keys %genestatus){
		if (!$gbMgene{$gene}) {
			push @order,$gene;
		}
	}
	open (OUT, ">$fOut/1.genestatus.txt") or die $!;
	print OUT "Gene";
	foreach my $sam (@sam) {
		print OUT "\t$sam";
	}
	print OUT "\n";
	foreach my $gene (@order) {
		print OUT "$gene";
		foreach my $sam (@sam) {
			if (!$genestatus{$gene}{$sam}) {
				$genestatus{$gene}{$sam}="4";
			}
			print OUT "\t$genestatus{$gene}{$sam}";
		}
		print OUT "\n";
	}
	my $count=2;
	foreach my $con (@code) {
		open (OUT, ">$fOut/$count.$con.level.txt") or die $!;
		print OUT "Gene";
		foreach my $sam (@sam) {
			print OUT "\t$sam";
		}
		print OUT "\n";
		foreach my $gene (@order) {
			print OUT "$gene";
			foreach my $sam (@sam) {
				if (!$genemlevel{$con}{$gene}{$sam}) {
					$genemlevel{$con}{$gene}{$sam}="NA";
				}
				else {
					$genemlevel{$con}{$gene}{$sam}=$genemlevel{$con}{$gene}{$sam}-1
				}
				print OUT "\t$genemlevel{$con}{$gene}{$sam}";
			}
			print OUT "\n";
		}
		close OUT;
		$count++;
	}
}
sub define_status {
	my @acces=glob"$fIn/*.out";
	foreach my $ac (@acces) {
		my $sam;
		if ($ac=~/(GSM\d+)\_/) {
			$sam=$1;
		}
		push @sam,$sam;
		my $cg;
		my $chg;
		my $chh;
		open (IN, $ac) or die $!;
		while (<IN>) {
			chomp;
			my @lines=split/\s+/,$_;
			next if (/^$/);
			$lines[1]=~s/\.\d+//;
			if ($lines[5]<1) {
				$genemlevel{$lines[3]}{$lines[1]}{$sam}=1;
			}
			else {
				$genemlevel{$lines[3]}{$lines[1]}{$sam}=$lines[4]/$lines[5]+1;
			}
			###record cg level
			$genemlevel{$lines[6]}{$lines[1]}{$sam}=$lines[17]+1;##record chg level
			$genemlevel{$lines[9]}{$lines[1]}{$sam}=$lines[20]+1;##record chh level
			if ($lines[11]>=20) {
				if ($lines[5]>=20&&$lines[14]<=0.05&&$lines[17]>0.05&&$lines[20]>0.05) {
					$genestatus{$lines[1]}{$sam}=1;
				}
				elsif ($lines[8>=20]&&$lines[17]<0.05&&$lines[20]>0.05) {
					$genestatus{$lines[1]}{$sam}=2;
				}
				elsif ($lines[11>=20]&&$lines[20]<0.05) {
					$genestatus{$lines[1]}{$sam}=3;
				}
				elsif ($lines[14]>0.05&&$lines[17]>0.05&&$lines[20]>0.05) {
					$genestatus{$lines[1]}{$sam}=5;
				}
				else {
					$genestatus{$lines[1]}{$sam}=4;###all undetermined situation due to limited reads counts are attribute to it status
				}
			}
			else {
				$genestatus{$lines[1]}{$sam}=4;
			}
		}
		close IN;
	}
}
#######################################################################################
my $Time_End   = sub_format_datetime(localtime(time()));
print "Program Ends Time:$Time_End\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
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


