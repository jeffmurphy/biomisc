#!/usr/bin/perl -w

use strict;
use Net::FTP;
use Getopt::Long qw(:config no_ignore_case bundling);
use URI;
use Bio::SeqReader::Fasta;
use Config::Simple;
use Data::Dumper;
use IO::File;
use IO::Uncompress::AnyUncompress;

my $VERSION = 2014011801;

my $_cfg = {
	'cache'               => "/tmp/nihprots",
	'debug'               => 0,
	'cache_lifespan'      => "7d",
	'url'                 => "ftp://ftp.ncbi.nih.gov/genomes",
	'user'                => 'anonymous',
	'pass'                => 'anonymous@',
	'patterns'            => [],
	'pattern_mode'        => 'any',
	'bacteria_species'    => [''],
	'nonbacteria_species' => ['.*'],
	'output_byspecies'    => 0,
	'output_bypattern'    => 1
};

my $config_file      = '';
my $help             = 0;
my $force_recache    = 1;

my $output_file      = undef;
my $output_byspecies = 0;
my $output_bypattern = 1;

my $version          = 0;

my $gor = GetOptions(
	"debug|d=i"    => \$_cfg->{debug},
	"config|C=s"   => \$config_file,
	"cache|c=s"    => \$_cfg->{cache},
	"url|U=s"      => \$_cfg->{url},
	"user|u=s"     => \$_cfg->{user},
	"pass|p=s"     => \$_cfg->{pass},
	"lifespan|l=s" => \$_cfg->{cache_lifespan},
	"mode|m=s"     => \$_cfg->{pattern_mode},
	"force|f"      => \$force_recache,
	"output|o=s"   => \$output_file,
	"byspecies"    => \$output_byspecies,
	"bypattern"    => \$output_bypattern,
	"version|v"    => \$version,
	"help|h|?"     => \$help
);

die "$VERSION\n" if $version;

$output_bypattern = 0 if $output_byspecies;

usage() if ( $help || !defined($output_file) );

if ( $config_file && -r $config_file ) {
	my $cfg = new Config::Simple($config_file);

	foreach my $_scalar (
		'cache', 'debug', 'cache_lifespan', 'url',
		'user',  'pass',  'pattern_mode', 'byspecies', 'bypattern'
	  )
	{
		if ( $cfg->param($_scalar) ) {
			D(1, "assigning $_scalar to cfg based on $config_file contents\n");
			$_cfg->{$_scalar} = $cfg->param($_scalar);
		}
	}

	foreach
	  my $_maybealist ( 'patterns', 'bacteria_species', 'nonbacteria_species' )
	{
		my @x = $cfg->param($_maybealist);
		$_cfg->{$_maybealist} = \@x;
	}
}

my @pats = @ARGV;
if ( $#pats > -1 ) {
	$_cfg->{'patterns'} = \@pats;
}

if ( !-d $_cfg->{cache} ) {
	mkdir $_cfg->{cache};
	$force_recache = 1;
}

my ( $host, $rdir ) = process_url( $_cfg->{url} );

clean_cache( $_cfg->{cache}, $_cfg->{cache_lifespan} );

if ( $force_recache
	|| cache_too_old( $_cfg->{cache}, $_cfg->{cache_lifespan} ) )
{
	D( 1, "refreshing cache in " . $_cfg->{cache} . "\n" );
	my $connection = connect_to_remote_host($host);
	cache_nonbacteria( $connection, $host, $rdir, $_cfg->{cache},
		$_cfg->{nonbacteria_species} );
	cache_bacteria( $connection, $host, $rdir . "/Bacteria",
		$_cfg->{cache}, $_cfg->{bacteria_species} );
	disconnect_from_remote_host($connection);
}
else {
	D( 1,
		    "not checking "
		  . $_cfg->{url}
		  . " because cache is less than "
		  . $_cfg->{cache_lifespan}
		  . " old\n" );
}

my $database = process_bacteria(process_nonbacteria());

if ($output_bypattern == 1) {
	D(1, "output_bypattern\n");
	write_bypattern($database);	
}
else {
	D(1, "output_byspecies\n");
	write_byspecies($database);
}

exit 0;

sub uniq {
    return keys %{{ map { $_ => 1 } @_ }};
}

sub write_bypattern {
	my $db = shift;
	
	my $fh = new IO::File $output_file, "w";
	my $output_buffer = {};
	
	if ( defined $fh ) {
		
		print $fh "species,pattern,#genes,#matches,percentage,species,species,pattern,#genes,#matches,percentage\n";

		foreach my $BorNB ('nonbac', 'bac') {
			my $match_count = $db->{$BorNB};
			my $species_count = $db->{"${BorNB}:species_count"};
			
			foreach my $species ( sort keys %$match_count ) {
				foreach my $pattern_matched ( sort keys %{ $match_count->{$species} } )
				{
					my $p = sprintf( "%8.8f",
						$match_count->{$species}->{$pattern_matched} /
						  $species_count->{$species} );
					
					$output_buffer->{$pattern_matched} = "" unless exists $output_buffer->{$pattern_matched};
					
					$output_buffer->{$pattern_matched} .= join( ',',
						$species,
						$species_count->{$species},
						$match_count->{$species}->{$pattern_matched}, $p ) . ",";
				}
			}
		}
		
		foreach my $pat (sort keys %$output_buffer) {
			print $fh $pat . "," . $output_buffer->{$pat} . "\n";
		}
		
		undef $fh;
	}
	else {
		die "failed to open $output_file for writing $!"
	}	
}

sub write_byspecies {
	my $db = shift;
	
	D(1, "writing to $output_file\n");
	
	my $fh = new IO::File $output_file, "w";
	if ( defined $fh ) {
		print $fh "species,pattern,#genes,#matches,percentage\n";

		foreach my $BorNB ('bac', 'nonbac') {
			my $match_count = $db->{$BorNB};
			my $species_count = $db->{"${BorNB}:species_count"};
			
			foreach my $species ( sort keys %$match_count ) {
				foreach my $pattern_matched ( sort keys %{ $match_count->{$species} } )
				{
					my $p = sprintf( "%8.8f",
						$match_count->{$species}->{$pattern_matched} /
						  $species_count->{$species} );
					print $fh join( ',',
						$species, $pattern_matched,
						$species_count->{$species},
						$match_count->{$species}->{$pattern_matched}, $p )
					  . "\n";
				}
			}
		}
		undef $fh;
	}
	else {
		die "failed to open $output_file for writing $!"
	}
}

sub extract_species {
	my $desc = shift;
	my $fn = shift;
	
	# match simple case where it's encoded as [...]
	if ( $desc =~ /\[([^\[\]]+)\]$/ ) {
		return $1;
	}
	
	# match [...[...]...] case
	if ( $desc =~ /\[(.*)\]$/g ) {
		return $1;
	}
	
	die "cant parse out [species] from \n\"" . $desc . "\"\n in file $fn\n";
}

sub process_nonbacteria {
	my $cache_dir = $_cfg->{cache};
	my $_seqpats  = $_cfg->{patterns};

	my $db = { };
	
	return $db if file_list_empty( $_cfg->{nonbacteria_species} );

	my $fnpat =
	  '^nonbac:(' . join( '|', @{ $_cfg->{nonbacteria_species} } ) . ')$';

	opendir( my $dh, $cache_dir ) || die "can't opendir $cache_dir: $!";
	my @fns = grep { /protein\.fa\.gz$/ && -f "$cache_dir/$_" } readdir($dh);
	closedir $dh;

	my $file_count    = 0;
	my $species_count = {};
	my $seq_count     = 0;
	my $match_count   = {};

	my $seqpats = join( '|', @$_seqpats );
	my @seqpats = ($seqpats);
	@seqpats = @$_seqpats if ( $_cfg->{pattern_mode} eq "each" );

	foreach my $fn ( sort @fns ) {
		if ( $fn =~ /$fnpat/ ) {
			$file_count++;
			my $fh3 =
			  new IO::Uncompress::AnyUncompress( $cache_dir . "/" . $fn )
			  or die "uncompress failed $!";
			my $in2 = new Bio::SeqReader::Fasta( fh => $fh3 );
			print "NB: $fn\n";
			while ( my $so = $in2->next() ) {
				my $speciesid = "";
				$seq_count++;
				my $desc = $so->desc();

				$speciesid = extract_species($desc, $fn);

				$species_count->{$speciesid}++;
				my $seq = $so->seq();

				foreach my $_pat (@seqpats) {
					$match_count->{$speciesid}->{$_pat} += 0;
					if ( $seq =~ /$_pat/ ) {
						$match_count->{$speciesid}->{$_pat}++;
					}
				}
			}
		}
	}

	D( 1, "nonbac filecount=$file_count seq_count=$seq_count\n" );

	$db->{'nonbac'} = $match_count;
	$db->{'nonbac:species_count'} = $species_count;
	
	return $db;
}

sub process_bacteria {
	my $cache_dir = $_cfg->{cache};
	my $_seqpats  = $_cfg->{patterns};

	my $db = shift;
	
	return $db if file_list_empty( $_cfg->{bacteria_species} );

	my $fnpat =
	  '^bacteria:(' . join( '|', @{ $_cfg->{bacteria_species} } ) . ')$';

	opendir( my $dh, $cache_dir ) || die "can't opendir $cache_dir: $!";
	my @fns = grep { /^bacteria:/ && -f "$cache_dir/$_" } readdir($dh);
	closedir $dh;

	my $file_count    = 0;
	my $species_count = {};
	my $seq_count     = 0;
	my $match_count   = {};

	my $seqpats = join( '|', @$_seqpats );
	my @seqpats = ($seqpats);
	@seqpats = @$_seqpats if ( $_cfg->{pattern_mode} eq "each" );

	foreach my $fn ( sort @fns ) {
		if ( $fn =~ /$fnpat/ ) {
			$file_count++;
			my $fh3 =
			  new IO::Uncompress::AnyUncompress( $cache_dir . "/" . $fn )
			  or die "uncompress failed $!";

			my $in2 = new Bio::SeqReader::Fasta( fh => $fh3 );

			print "BAC: $fn\n";

			while ( my $so = $in2->next() ) {
				my $speciesid = "";
				$seq_count++;
				my $desc = $so->desc();

				$speciesid = extract_species($desc, $fn);

				$species_count->{$speciesid}++;

				my $seq = $so->seq();

				foreach my $_pat (@seqpats) {
					$match_count->{$speciesid}->{$_pat} += 0;
					if ( $seq =~ /$_pat/ ) {
						$match_count->{$speciesid}->{$_pat}++;
					}
				}
			}
		}
	}

	D( 1, "bac filecount=$file_count seq_count=$seq_count\n" );
	
	$db->{'bac'} = $match_count;
	$db->{'bac:species_count'} = $species_count;
	
	return $db;
}

sub file_list_empty {
	my $a = shift;
	return 1 if ( $#$a == -1 || ( $#$a == 0 && $a->[0] eq '' ) );
	return 0;
}

sub connect_to_remote_host {
	my $host = shift;

	D( 1, "connecting to $host\n" );

	my $ftp = Net::FTP->new( $host, Debug => 0 )
	  or die "Can not connect to $host: $@";

	$ftp->login( "anonymous", "anonymous@" )
	  or die "Can not login as anonymous to $host";

	$ftp->binary;

	return $ftp;
}

sub disconnect_from_remote_host {
	D( 1, "disconnecting from remote host\n" );
	shift->quit;
}

sub cache_lifespan_to_seconds {
	my $cache_lifespan = shift;

	my $scale_factor = {
		'm' => 60,
		'h' => 60 * 60,
		'd' => 60 * 60 * 24,
		'w' => 60 * 60 * 24 * 7
	};

	D( 2, "cache_lifespan = $cache_lifespan\n" );

	if ( $cache_lifespan =~ /(\d+)([mhdw])$/ ) {
		$cache_lifespan = $1 * $scale_factor->{$2};
	}
	else {
		$cache_lifespan = int($cache_lifespan);
	}

	D( 2, "cache_lifespan = $cache_lifespan (seconds)\n" );
	return $cache_lifespan;
}

sub cache_too_old {
	my $cache_dir      = shift;
	my $cache_lifespan = cache_lifespan_to_seconds(shift);

	my $cache_dir_time = ( stat($cache_dir) )[9];

	D( 2, "cache_dir_time = $cache_dir_time\n" );

	return 1 if ( time() - $cache_dir_time > $cache_lifespan );
	return 0;
}

=pod

cache_normal(connection, remote_host, remote_dir, local_cache_dir)

1. connect to remote_host
2. fetch listing of dirs in remote_dir
3. foreach dir
   a. fetch size of remote_dir/dir/protein.fa.gz
   b. if size differs from local copy update local copy

=cut

sub cache_nonbacteria {
	my ( $ftp, $host, $rdir, $cache, $species ) = @_;
	$| = 1;

	return if file_list_empty($species);

	$ftp->cwd($rdir)
	  or die "Can not change to remote dir $rdir ", $ftp->message;

	my @ls = $ftp->ls();

	chdir($cache);

	my $fnpat = '^(' . join( '|', @{$species} ) . ')$';
	D( 2, "only fetching non-bacterial species that match: $fnpat\n" );

	for my $rd ( sort @ls ) {

		next if $rd =~ /Bacteria/;    # skip bacteria folders
		next if $rd !~ /$fnpat/;      # dont fetch things we arent interested in

		my $rprotfile = $rdir . "/" . $rd . "/protein/protein.fa.gz";
		my $lprotfile = $cache . "/nonbac:" . $rd . "-protein.fa.gz";

		print "Checking for $rprotfile\t";
		my $rsize = $ftp->size($rprotfile);

		if ( $rsize && $rsize > 0 ) {
			my $lsize = 0;
			$lsize = ( stat($lprotfile) )[7] if ( -r $lprotfile );
			print "[local $lsize remote $rsize]\t";

			if ( $lsize != $rsize ) {
				print "[updating cached copy]\t";
				$ftp->get( "$rd/protein/protein.fa.gz", $lprotfile );
				print "[done]\n";
			}
			else {
				print "[cache still valid]\n";
			}
		}
		else {
			print "[remote doesnt exist]\n";
		}
	}

}

=pod

cache_bacteria(connection, remote_host, remote_dir, local_cache_dir)

1. connect to remote_host
2. fetch listing of dirs in remote_dir
3. foreach file in above listing
   a. fetch listing of /that/ dir
   b. foreach file in /that/ listing
      i. if file ends in .faa and size differs from local copy, update local copy

=cut

sub cache_bacteria {
	my ( $ftp, $host, $rdir, $cache, $species ) = @_;
	$| = 1;

	return if file_list_empty($species);

	$ftp->cwd($rdir)
	  or die "Can not change to remote dir $rdir ", $ftp->message;

	my @ls = $ftp->ls();

	chdir($cache);

	my $fnpat = '^(' . join( '|', @{$species} ) . ')$';

	D( 2, "only fetching bacterial species that match: $fnpat\n" );

	for my $rd ( sort @ls ) {
		next if $rd !~ /$fnpat/;    # dont fetch things we arent interested in

		my @innerdir = $ftp->ls( $rdir . "/" . $rd );

		foreach my $innerfile ( sort @innerdir ) {
			if ( $innerfile =~ /\.faa$/ ) {
				my $rbacfile = $innerfile;
				my @rbacparts = split( '/', $rbacfile );
				my $lbacfile =
				  $cache . "/"
				  . join( ":",
					"bacteria", $rbacparts[ $#rbacparts - 1 ],
					$rbacparts[$#rbacparts] );

				print "Checking for $rbacfile\t";

				my $rsize = $ftp->size($rbacfile);

				if ( $rsize && $rsize > 0 ) {
					my $lsize = 0;
					$lsize = ( stat($lbacfile) )[7] if ( -r $lbacfile );
					print "[local $lsize remote $rsize]\t";

					if ( $lsize != $rsize ) {
						print "[updating cached copy]\t";
						$ftp->get( $rbacfile, $lbacfile );
						print "[done]\n";
					}
					else {
						print "[cache still valid]\n";
					}
				}
				else {
					print "[remote doesnt exist]\n";
				}
			}
		}
	}

}

sub clean_cache {
	my ( $cache_dir, $age ) = ( shift, shift );
	return unless ( -d $cache_dir );

	D( 2, "purging cache of files older than $age\n" );

	$age = cache_lifespan_to_seconds($age);

	opendir( my $dh, $cache_dir ) || die "can't opendir $cache_dir: $!";
	foreach my $fn ( grep { /\.(faa|gz)$/ && -f "$cache_dir/$_" } readdir($dh) )
	{
		my $lmtime = 0;
		$lmtime = ( stat("$cache_dir/$fn") )[9] if ( -r "$cache_dir/$fn" );
		my $fileage = time() - $lmtime;
		if ( $fileage > $age ) {
			D( 3, "\tdeleting $fn\n" );
			unlink("$cache_dir/$fn") || warn "\t\tfailed to delete $fn -- $!";
		}
	}
	closedir $dh;
}

sub D {
	my ( $level, $msg ) = ( shift, shift );
	print $msg if $_cfg->{debug} >= $level;
}

sub process_url {
	my $url = shift;
	return undef unless $url;
	my $u = URI->new($url);
	die "we can only handle FTP right now\n" if ( $u->scheme ne "ftp" );
	D( 1, "host: " . $u->host . " path: " . $u->path . "\n" );
	return ( $u->host, $u->path );
}

sub usage {
	print "$0 [options] <-o outputfile> [pattern] [pattern] [...]\n
	--output <file>   file to write results to (required; shortcut -o) **REQUIRED

	--bypattern		  group output by pattern (default)
	--byspecies		  group output by species
	
	--debug #         verbose output (shortcut: -d)
	--config <file>   use this file instead of command line options (shortcut: -C)
	--cache <dir>     cache folder (shortcut -c)
	--url <url>       url location of genome data (ftp:// only; shortcut -U)
	--user <user>     login with this username (default: anonymous shortcut: -u)
	--pass <passwd>   password to login with (default: anonymous@ shortcut: -p)
	--lifespan T      maximum cache age before we will refresh the cache (default: 7d shortcut: -l)
	    T = #[mhdw]     specify cache age in minutes/hours/days/weeks
	--force           force a recache (same as -l 0; shortcut -f)
	--version         print version (shortcut -v)
	--help            this message (shortcut: -h -?)
	
	version: $VERSION
	\n";
	exit 0;
}

# I tried using Bio::SeqReader::Fasta but it didnt parse things cleaning, threw errors
# and didnt work at all on the bac files for some reason so it was just quicker to roll
# an analog

package FastaReader;
use Fcntl qw(:flock SEEK_CUR);

sub new {
	my $class = shift;
	my $self  = {};
	bless $self, $class;

	$self->{fh} = shift;

	return $self;
}

sub next {
	my $self = shift;

	my $line = "";

	my $done = 0;
	my $desc = '';

	my $fh = $self->{fh};

	while ( !$done && !eof( $self->{fh} ) ) {
		$line = <$fh>;
		if ( $line =~ /^>(.*)/ ) {
			$done = 1;
			$desc = $1;
		}
	}

	$done = 0;
	my @s = ();

	while ( !$done && !eof( $self->{fh} ) ) {
		$line = <$fh>;
		if ( $line =~ /^>/ ) {
			$done = 1;
			seek( $self->{fh}, SEEK_CUR, -length($line) );
		}
		else {
			push @s, $line;
		}
	}

	my $fr = new FastaRecord( $desc, join( '', @s ) );

}

package FastaRecord;

sub new {
	my $class = shift;
	my $self  = {};
	bless $self, $class;
	$self->{desc} = shift;
	$self->{seq}  = shift;
	return $self;
}

sub desc {
	my $self = shift;
	return $self->{desc};
}

sub seq {
	my $self = shift;
	return $self->{seq};
}

1;
