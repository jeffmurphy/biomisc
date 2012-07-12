#!/usr/bin/perl -w
# jeff murphy <jcmurphy@nickelsoft.com>
use strict;
use Data::Dumper;
use Math::Random::OO::Normal;
use Math::Random::OO::Uniform;

my $cfgfn = shift;
$cfgfn ||= "randomaa.txt";

my $cfg = read_config($cfgfn);
$cfg->{dtable} = make_dist_table($cfg->{aadist});

print Dumper($cfg) . "\n";

for(my $nf = 1; $nf <= $cfg->{files} ; $nf++) {
    my $ofn = $cfg->{fileprefix} . "${nf}.fasta";
    print "Writing file: $ofn\n";
    my $seqs = make_seqs($cfg);

    write_file($ofn, $seqs);
}

exit 0;

# write_file(filename, arrayref of sequences)

sub write_file {
    my $fn = shift;
    my $ss = shift;
    open (OFD, "> $fn") || die "cant open $fn for writing: $!";
    for (my $i = 0 ; $i <= $#$ss; $i++) {
	print OFD "> randomseq" . ($i+1) . "\n";
	print OFD $ss->[$i] . "\n";
    }
    close(OFD);
}

# make_seqs(confighash)

sub make_seqs {
    my $c = shift;
    my @s;

    my $rng;
    if ($c->{lendist} eq "random") {
	$rng = new Math::Random::OO::Uniform(
	    $c->{minlen}, $c->{maxlen}
	    );
    }
    elsif ($c->{lendist} eq "normal") {
	 $rng = new Math::Random::OO::Normal(
	     $c->{mean},
	     $c->{stddev}
	     );
    }

    $rng->seed(seedval());

    # for each sequence
    for (my $i = 0 ; $i < $c->{seqsperfile} ; $i++) {
	# sequence length.
	my $sl = 0;
	my $isgood = 0;
	while (!$isgood) {
	    $sl = int($rng->next());
	    # just incase 'normal' gives us an outlier due to a large
	    # stddev value
	    $isgood = 1 if ($sl >= $c->{minlen} && $sl <= $c->{maxlen});
	}
	push @s, make_seq($sl, $c->{dtable});
    }
    return \@s;
}

# make_seq(length, distribution table)

sub make_seq {
    my $l = shift;
    my $tbl = shift;

    my $s = '';
    for (my $i = 0 ; $i < $l ; $i++) {
	my $rn = int(rand(100));
	$s .= $tbl->[$rn];
    }
    return $s;
}


# given a distribution like:
#
# a 0.25
# b 0.50
# c 0.20
# d 0.10
#
# produce an array like:
#
# [ a, a, a (25 times),
#   b, b, b, (50 times)
#   c, c, c, (20 times),
#   d, d, d, (10 times) ]
#
# allowing us to pick a random number (1-100) and map it into
# our distribution
#

sub make_dist_table {
    my $d = shift;
    my $a;
    foreach my $aa (sort keys %$d) {
	my $v = $d->{$aa} * 100.0;
	$a .=  $aa x $v;
    }
    return [split(//, $a)];
}


sub read_config {
    my $fn = shift;
    my $c = { # default values
	'files' => 33,
	'fileprefix' => 'randomaa-',
	'seqsperfile' => 22,
	'minlen' => 100,
	'maxlen' => 10000,
	'mean' => 4900,
	'stddev' => 1000,
	'lendist' => 'random',
	'aadist' => {
	    'a' => '0.05',
	    'r' => '0.05',
	    'n' => '0.05',
	    'd' => '0.05',
	    'c' => '0.05',
	    'e' => '0.05',
	    'q' => '0.05',
	    'g' => '0.05',
	    'h' => '0.05',
	    'i' => '0.05',
	    'l' => '0.05',
	    'k' => '0.05',
	    'm' => '0.05',
	    'f' => '0.05',
	    'p' => '0.05',
	    's' => '0.05',
	    't' => '0.05',
	    'w' => '0.05',
	    'y' => '0.05',
	    'v' => '0.05',
	},
    };

    open(FD, $fn) || die "cant open config file: $fn error: $!";
    while(my $line = <FD>) {
	chomp $line;
	$line =~ s/#.*$//g;

	if ($line =~ /files:\s+(\d+)/) {
	    $c->{files} = $1;
	}
	elsif ($line =~ /fileprefix:\s+(\S+)/) {
	    $c->{fileprefix} = $1;
	}
	elsif ($line =~ /seqs.file:\s+(\d+)/) {
	    $c->{seqsperfile} = $1;
	}
	elsif ($line =~ /min.max length:\s+(\d+)\s+(\d+)/) {
	    $c->{minlen} = $1;
	    $c->{maxlen} = $2;
	} 
	elsif ($line =~ /normal mean.stddev:\s+(\d+)\s+(\d+)/) {
	    $c->{mean} = $1;
	    $c->{stddev} = $2;
	}
	elsif ($line =~ /length distribution:\s+(random|normal)/) {
	    $c->{lendist} = $1;
	}
	elsif ($line =~ /aa distribution:/) {
	    while($line = <FD>) {
		if ($line =~ /(\S)\s+([\d\.]+)/) {
		    $c->{aadist}->{$1} = $2;
		}
	    }
	}
	elsif ($line =~ /^\s*$/) {
	    next;
	}
	else {
	    die "not sure what this line from $fn means: $line";
	}
    }

    die "max len < min len in $fn"
	if ($c->{maxlen} < $c->{minlen});

    die "aa distribution doesnt sum to 1.00"
	unless (aadist_sanity($c->{aadist}));

    close(FD);
    return $c;
}

sub aadist_sanity {
    my $d = shift;
    my $s = 0;
    foreach my $aa (keys %$d) {
	$s += $d->{$aa};
    }
    print "s $s\n";
    return ($s != 1.0);
}

sub seedval {
        open(FD, "/dev/random") || die "$!";
        my $x;
        my $l = sysread(FD, $x, 4);
        close(FD);

        my $a = 0;
        for(my $i = 0 ; $i < $l ; $i++) {
                $a += unpack('C', substr($x, $i, 1));
        }
        return $a;
}
