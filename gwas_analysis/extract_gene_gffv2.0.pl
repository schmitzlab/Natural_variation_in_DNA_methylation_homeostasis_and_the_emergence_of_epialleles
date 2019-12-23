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
my ($fIn1,$fIn2,$fIn3,$fIn4,$fOut);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"gff:s"=>\$fIn1,
				"anno:s"=>\$fIn2,
				"snpeff:s"=>\$fIn3,
				"pvalue:s"=>\$fIn4,

				) or &USAGE;
&USAGE unless ($fIn1 and $fIn2 and $fIn3 and $fIn4 and $fOut);

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	zhangyw <yz46606\@uga.edu>
Program Date:   2016/09/02
Description:	this program is used to annotated the snp around gwas peaks, to see if
				snp is within genes, snp annotation and gene annotation
Usage:
  Options:
  -gff <file>  input file,gff format,forced
  -anno <file> input file, gene annotation fileforced
  -snpeff <file>  input file,snp annotation information,forced
  -pvalue <dir>  gapit gwas result
  -o <dir>  output file,genelist forced
  -h         Help

USAGE
	print $usage;
	exit;
}

mkdir $fOut if (! -d $fOut) ;
$fOut=&ABSOLUTE_DIR($fOut);

##first extract the significant regions
my $limit=1e-5;
my $merge=10000;
my @gresult=glob"$fIn4/*.GWAS.Results.csv";
foreach my $gwas (sort @gresult) {
	my $name=basename($gwas);
	my %region=();
	my %sig_snp=();
	my %pvalue=();
	open (IN, $gwas) or die $!;
	while (<IN>) {
		chomp;
		next if ($.==1||/^$/);
		my @lines=split/\,/,$_;
		if ($lines[3]<$limit) {
			push @{$sig_snp{$lines[1]}},$lines[2];
			$pvalue{$lines[1]}{$lines[2]}=$lines[3];
		}
		else{last;}
	}
	close IN;
	my %sort_sigsnp;
	foreach my $key (sort {$a<=>$b} keys %sig_snp) {
		my $count=1;my $start;my $end;my $region;
		foreach my $loci (sort{$a<=>$b} @{$sig_snp{$key}}) {
			push @{$sort_sigsnp{$key}},$loci;
			#print $loci;die;
			if ($count==1) {
				$start=$loci;
				$end=$loci;
				$count++;
			}
			else {
				if ($loci-$end>$merge) {
					$region++;
					$region{$key}{$region}="$start $end";
					$start=$loci;
					$end=$loci;
					$count++;
				}
				else {
					$end=$loci;
					$count++;
				}
			}
		}
		$region++;
		$region{$key}{$region}="$start\t$end";
	}
	#print Dumper %region;die;
	my %gene;
	my %geneid;
	open (IN, $fIn1) or die $!;
	##second extract geneID from gff
	while (<IN>) {
		chomp;
		next if ($.==1||/^$/||/^\#/);
		my @lines=split/\s+/,$_;
		if ($lines[2]=~/gene/&&$lines[8]=~/Name=([A-Za-z0-9\.]{1,})\;{0,}/) {
			my $name=$1;
			push @{$gene{$lines[0]}},"$lines[3]\t$lines[4]";
			$geneid{$lines[0]}{"$lines[3]\t$lines[4]"}=$name;
		}
	}
	close IN;
	##from annotation
	my %anno;
	open (IN, $fIn2) or die $!;
	while (<IN>) {
		chomp;
		next if ($.==1||/^$/||/^\#/);
		my @lines=split/\s+/,$_,5;
		$anno{$lines[1]}=$lines[4];
	}
	close IN;
	##################################################define reliable mQTL region
	my %fRegion;###filterRegion;
	foreach my $key (sort {$a<=>$b} keys %region ) {
		foreach my $region (sort {$a<=>$b} keys %{$region{$key}}) {
			my ($start,$end)=split/\s+/,$region{$key}{$region};
			next if ($end-$start<100);
			push @{$fRegion{$key}},$region{$key}{$region};
		}
	}
	#print Dumper %fRegion;
	my %region_snp;
	my %snp_inregion;
	################################record snp point in each defined region;
	foreach my $chr (sort {$a<=>$b}keys %fRegion) {
		#print $chr,"\n";
		my $shift=0;
		MM:foreach my $region (@{$fRegion{$chr}}) {
			my ($start,$end)=split/\s+/,$region;
			for (my $i=$shift;$i<@{$sort_sigsnp{$chr}};$i++) {
				#print $sort_sigsnp{$chr}[$i];die;
				if ($sort_sigsnp{$chr}[$i]>=$start&&$sort_sigsnp{$chr}[$i]<=$end) {
					push @{$region_snp{$chr}{$region}},$sort_sigsnp{$chr}[$i];
					$snp_inregion{$chr}{$sort_sigsnp{$chr}[$i]}++;
					$shift=$i+1;
				}
				elsif ($sort_sigsnp{$chr}[$i]>$end) {
					$shift=$i;
					next MM;
				}
			}
		}
	}
	#################record SNPEFF files to get the snp annotation information
	my %snpeffinfo;
	open (IN, $fIn3) or die $!;
	while (<IN>) {
		chomp;
		next if (/\#/);
		my @lines=split/\s+/,$_;
		if ($snp_inregion{$lines[0]}{$lines[1]}) {
			if ($lines[7]=~/EFF\=(.*)/) {
				$snpeffinfo{$lines[0]}{$lines[1]}=$1;
			}
		}
	}
	close IN;

	###################find genes that having snp in it
	my %snp_gene;
	foreach my $chr (sort {$a<=>$b} keys %snp_inregion) {
		my $shift=0;
		MM:foreach my $snp (sort {$a<=>$b} keys %{$snp_inregion{$chr}}) {
			for (my $i=$shift;$i<@{$gene{$chr}};$i++) {
				my ($start,$end)=split/\s+/,$gene{$chr}[$i];
				if ($snp>=$start&&$snp<=$end) {
					$snp_gene{$chr}{$snp}=$gene{$chr}[$i];
					$shift=$i;
					next MM;
				}
			}
		}
	}
	close OUT;
	#######################output
	open (OUT, ">$fOut/$name.genelist.txt") or die $!;
	print OUT "Chromosome\tpeakRegion\tsnp\tsnpAnnotation\tgene\tgeneLocation\tgeneAnnotation\n";
	foreach my $chr (sort {$a<=>$b} keys %fRegion) {
		foreach my $region (@{$fRegion{$chr}}){
			foreach my $snp (@{$region_snp{$chr}{$region}}) {
				print OUT "$chr\t$region\t$snp\t$snpeffinfo{$chr}{$snp}\t$snp_gene{$chr}{$snp}\t$geneid{$chr}{$snp_gene{$chr}{$snp}}\t$anno{$geneid{$chr}{$snp_gene{$chr}{$snp}}}\n";
			}
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
	#���б��е�����ֵ
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
	#���б��е���Сֵ
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
	#��ȡ�ַ������еķ��򻥲����У����ַ�����ʽ���ء�ATTCCC->GGGAAT
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
