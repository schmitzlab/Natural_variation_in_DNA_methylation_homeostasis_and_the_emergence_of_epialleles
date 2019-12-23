#!/usr/bin/perl -w 

## This script is used for RNS_seq data analysis
## This script is used to remove rRNA seqs 
## This script was written by Lexiang Ji, 13-Mar-2015
## This script was updated by Lexiang Ji, 26-Nov-2018

$blast  = $ARGV[0];
$fastq  = $ARGV[1];
$output = $ARGV[2];

open AA, "$blast";
open BB, "$fastq";
open CC, ">$output";
select CC;

## capture the name of mapped reads
while ($line = <AA>) {
    $line =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
    $rRNA_seq{$10} = 'yes';
}

### remove r/tRNA reads 
while ($line1 = <BB>) {

    $line2 = <BB>;
    $line3 = <BB>;
    $line4 = <BB>;

    $line1 =~ /^@(\S+)/;
    $seq   = $1;
    if (!($rRNA_seq{$seq})){
    print "$line1";
    print "$line2";
    print "$line3";
    print "$line4";

    }
}

