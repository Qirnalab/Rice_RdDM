#!/usr/bin/perl

my ($anno,$set)=@ARGV;

open IM,$anno;

while(<IM>){
	chomp;
	@a=split /\t/,$_;
	$hash{@a[0]}=@a[1]."\t".@a[2];
}

open IN,$set;

while(<IN>){
	chomp;
	@b=split /\t/,$_;
	$tag= defined($hash{@b[0]}) ? $hash{@b[0]} : "NA\tNA";
	print "$_\t$tag\n";
}
