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
my ($fOut,$list,$blast,$eu,$lyr);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"l:s"=>\$list,
				"b:s"=>\$blast,
				"fa1:s"=>\$eu,
				"fa2:s"=>\$lyr,
				) or &USAGE;
&USAGE unless ($fOut and $list and $blast and $eu and $lyr);

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	zhangyw <yz46606\@uga.edu> 
Program Date:   2016/09/11
Description:	this program is used to convert dmr_original list to standard gabit trait form
Usage:
  Options:
  -l <file>  input file,genelist,forced  
  -b <file>  input file,blast.out,forced  
  -fa1 <files>  input file,Eu fa files,forced 
  -fa2 <files>  input file,Alyrata fa files,forced
  -o <dir>  output file,forced  
  -h         Help

USAGE
	print $usage;
	exit;
}
mkdir $fOut if (! -d $fOut);
$fOut=&ABSOLUTE_DIR($fOut);
################################################1. find the orthologous gbm genes between lyrata and thalina
my %gene;
#&read_genes();
#&find_orthologous();

#################################################2. alignment prank
#&prank_alignment();

##############################################3.gblock
#&gblock();
#&convert_phylip();

###############################################4. yz00-palm estimate dn/ds ratio
&estimate_dnds();

################################################5.extract w 
&extract_w();

sub extract_w{
	$/="\n";
	my @files=glob"$fOut/yn00/*.out";
	my %omega;my %dn;my %ds;
	foreach my $file (@files) {
		open (IN, $file) or die $!;
		my $name=basename($file);
		(my $gene)=($name)=~/\_(.*?).\d+\_/;
		my $count;
		while (<IN>) {
			chomp;
			$count++;
			next if (/^$/);
			if ($count==91) {
				my @a=split/\s+/,$_;
				$omega{$gene}=$a[7];
				$dn{$gene}=$a[8];
				$ds{$gene}=$a[11];
			}
		}
		close IN;
	}
	open (OUT, ">$fOut/2.add.dnds.gene.features.txt") or die $!;
	open (IN, $list) or die $!;
	while (<IN>) {
		chomp;
		next if (/^$/);
		my @lines=split/\s+/,$_;
		if ($.==1){
			print OUT "$_\tdnds\n";
			next;
		}
		
		if ($omega{$lines[0]}) {
			print OUT "$_\t$omega{$lines[0]}\n";
		}
		else {
			print OUT "$_\tNA\n";
		}
	}
	close IN;
}

###################################################sub
sub estimate_dnds{
	my @files=glob"$fOut/phylip/*.nuc";
	my $yn="$fOut/yn00";
	mkdir $yn if (! -d $yn);
	foreach my $file (@files) {
		my $name=basename($file);
		$name=~s/\.nuc//;
		open (OUT, ">$yn/$name.ctl") or die $!;
		print OUT "seqfile = $file\n";
		print OUT "outfile = $yn/$name.out\n";
		print OUT "verbose = 0\n";
		print OUT "icode = 0\n";
		print OUT "weighting = 0\n";
		print OUT "commonf3x4 = 0\n";
		print OUT "ndata = 1\n";
		close OUT;
		`/lustre1/yz46606/01_bin/paml/paml4.9h/bin/yn00 $yn/$name.ctl`;
	}
}
sub convert_phylip{
	my $gblock="$fOut/gblock";
	mkdir $gblock if (! -d $gblock);
	`mv $fOut/alignment/*-gb $gblock/`;
	my @files=glob"$gblock/*-gb";
	mkdir "$fOut/phylip" if (! -d "$fOut/phylip");
	$/="\>";
	foreach my $file (@files) {
		my $name=basename($file);
		$name=~s/.fa.best.fas-gb//;
		open (IN, $file) or die $!;
		open (OUT, ">$fOut/phylip/$name.nuc") or die $!;
		my $block;
		my $sam;my $seq;
		while (<IN>) {
			chomp;
			next if (/^$/);
			my @lines=split/\n/,$_,2;
			$sam++;
			$block.="\n$lines[0]\n$lines[1]";
			if ($sam==1) {
				my @len=split//,$lines[1];
				foreach my $l (@len) {
					if ($l=~/\w/) {
						$seq++;
					}
				}
			}
		}
		close IN;
		print OUT " $sam $seq\n";
		print OUT "$block";
		close OUT;
	}
}
sub gblock {
	my @files=glob "$fOut/alignment/*.fas";
	foreach my $file (@files) {
		`/usr/local/apps/gb/gblocks/0.91b/Gblocks $file -t=C -b3=8 -b5=a -e gb`;
	}
}
sub prank_alignment {
	my $aligndir="$fOut/alignment";
	mkdir $aligndir if (! -d $aligndir);
	open (OUT, ">$aligndir/align.sh") or die $!;
	mkdir $aligndir if (! -d $aligndir);
	my @files=glob "$fOut/pairfasta/*.fa";
	foreach my $file (@files) {
		my $name=basename($file);
		my $outfile="$aligndir/$name";
		my $cmd="/lustre1/yz46606/01_bin/evoluation/prank/bin/prank -d=$file -o=$outfile -codon -F";
		print OUT "$cmd\n";
	}
	close OUT;
	
}
sub find_orthologous{
	open (IN, $blast) or die $!;
	my @pair;my %highest;my %eu;my %lyr;
	while (<IN>) {
		chomp;
		my @lines=split/\s+/,$_;
		next if (/^$/);
		my $gene=$lines[0];
		$gene=~s/\.\d+//;
		if (!$highest{$gene}&&$gene{$gene}&&$lines[-1]>500) {
			push @pair,"$lines[0]\t$lines[1]";
			$highest{$gene}++;
			$eu{$lines[0]}++;
			$lyr{$lines[1]}++;
		}
	}
	close IN;
	$/="\>";my %faeu;my %falyr;
	open (IN, $eu) or die $!;
	while (<IN>) {
		chomp;
		next if (/^$/);
		my @lines=split/\n/,$_,2;
		(my $id)=($lines[0])=~/(\S+)\s/;
		if ($eu{$id}) {
			$faeu{$id}=$lines[1];
		}
	}
	close IN;

	open (IN, $lyr) or die $!;
	while (<IN>) {
		chomp;
		next if (/^$/);
		my @lines=split/\n/,$_,2;
		(my $id)=($lines[0])=~/(\S+)\s/;
		if ($lyr{$id}) {
			$falyr{$id}=$lines[1];
		}
	}
	close IN;

	my $pairfasta="$fOut/pairfasta";
	mkdir $pairfasta if (! -d $pairfasta);
	my $count;
	foreach my $pair (@pair) {
		my ($euid,$lyrid)=split/\s+/,$pair;
		$count++;
		open (OUT, ">$pairfasta/$count\_$euid\_$lyrid.fa") or die $!;
		print OUT ">$euid\n";
		print OUT "$faeu{$euid}";
		print OUT ">$lyrid\n";
		print OUT "$falyr{$lyrid}";
		close OUT;
	}
}
sub read_genes{
	open (IN, $list) or die $!;
	while (<IN>) {
		chomp;
		next if ($.==1);
		next if (/^$/);
		my @lines=split/\s+/,$_;
		$gene{$lines[0]}=$_;
	}
	close IN;
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


