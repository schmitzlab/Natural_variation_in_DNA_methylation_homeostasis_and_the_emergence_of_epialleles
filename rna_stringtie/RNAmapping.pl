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
Program Date:   2016/09/11
Description:	this program is used to convert dmr_original list to standard gabit trait form
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
my $gff_file="/lustre1/yz46606/02_reference/arab/Athaliana_447_Araport11.gene.Chr.gff3";
my $index="/lustre1/yz46606/02_reference/arab/TAIR10_chr_all.Chr";
my $bin="/lustre1/yz46606/01_bin/rna_stringtie";
my $ref_genome="/lustre1/yz46606/02_reference/arab/TAIR10_chr_all.Chr.fa";
my $library="F";
my $intron_min_size=20;
my $intron_max_size=10000;
my $rRNA="/lustre1/yz46606/02_reference/arab/Araport11_tRNA_rRNA_etal_uniq.fa";
my $CPU=3;
$fIn=&ABSOLUTE_DIR($fIn);
$fOut=&ABSOLUTE_DIR($fOut);

my @data=glob "$fIn/*.fastq";
my $count;
foreach my $data (@data) {
	my $name=basename($data);
	$name=~s/\.fastq//;
	$count++;
	open (OUT, ">$fOut/$name.$count.sh") or die $!;
	print OUT "#!/bin/bash\n";
	print OUT "#PBS -N $count\n";
	print OUT "#PBS -l walltime=48:00:00\n";
	print OUT "#PBS -l nodes=1:ppn=$CPU:AMD\n";
	print OUT "#PBS -l mem=50g\n";
	print OUT "#PBS -q batch\n";
	print OUT "cd \$PBS_O_WORKDIR\n";
	print OUT "module load Trimmomatic/0.36-Java-1.8.0_144\n";
	print OUT "module load HISAT2/2.1.0-foss-2016b\n";
	print OUT "module load StringTie/1.3.4d-foss-2016b\n";
	print OUT "module load BEDTools/2.26.0-foss-2016b\n";
	print OUT "module load seqtk/1.2-foss-2016b\n";
	print OUT "module load BLAT/3.5-foss-2016b\n";
	print OUT "module load SAMtools/1.9-foss-2016b\n";

	## trim adaptor
	print OUT "java -jar \$EBROOTTRIMMOMATIC/trimmomatic-0.36.jar SE -phred33 $data $fOut/Step_1_$name\_trimmed.fastq ILLUMINACLIP:/usr/local/apps/trimmomatic/0.33/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36\n";
	
	## transform fastq to fasta
	print OUT "seqtk seq -a $fOut/Step_1_$name\_trimmed.fastq >$fOut/Step_2_$name\_trimmed.fasta\n";
	
	## blat with rRNA database
	print OUT "blat $rRNA $fOut/Step_2_$name\_trimmed.fasta -t=dna -q=dna -noHead $fOut/Step_3_$name.blat\n";

	## remove rRNA reads
	print OUT "perl $bin/Script_03_remove_rRNA.pl $fOut/Step_3_$name.blat $fOut/Step_1_$name\_trimmed.fastq $fOut/Step_4_$name\_trimmed_remove_rRNA.fastq\n";
	print OUT "rm -rf $fOut/Step_3_$name.blat\n";
	
	## hisat2 mapping
	print OUT "hisat2 -p 1 --dta -x $index --min-intronlen $intron_min_size --max-intronlen $intron_max_size --rna-strandness $library -k 2 -U $fOut/Step_4_$name\_trimmed_remove_rRNA.fastq -S $fOut/Step_5_$name.sam\n";
	
	## samtools
	print OUT "samtools view -Su $fOut/Step_5_$name.sam | samtools sort  -@ 4 -m 5G -o $fOut/Step_06_$name.sorted.bam\n";
	
	## stringtie
	print OUT "stringtie $fOut/Step_06_$name.sorted.bam --rf -G $gff_file -A $fOut/Step_07/$name/$name.gene_abund.tab -B -e -o $fOut/Step_07/$name/$name.gtf\n";
	close OUT;
}
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


