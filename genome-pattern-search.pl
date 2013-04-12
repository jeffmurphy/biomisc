#!/usr/bin/perl -w
use strict;
use Net::FTP;

my $cache = "/tmp/nihprots";
if (! -d $cache) {
	mkdir $cache;
}

my $host = "ftp.ncbi.nih.gov";
my $rdir = "/genomes";


my $ftp = Net::FTP->new($host, Debug => 0) 
	or die "Can not connect to $host: $@";

$ftp->login("anonymous", "anonymous@")
	or die "Can not login as anonymous to $host";

$ftp->cwd($rdir)
	or die "Can not change to remote dir $rdir ", $ftp->message;

my @ls = $ftp->ls();

chdir($cache);
for my $rd (@ls) {
	print "$rd\n";
	$ftp->get("$rd/protein/protein.fa.gz", "$rd-protein.fa.gz");
}
$ftp->quit;

