#!/usr/bin/perl

use File::Find::Rule;
use Bio::DB::EUtilities;
use DB_File;

my $dir = shift @ARGV || die "Provide PUTs repository root dir";

my %gi_acc_map;

tie %gi_acc_map, 'DB_File', "$dir/.puts.acc.db", O_CREAT|O_RDWR;

if (scalar(keys(%gi_acc_map)) == 0) {
    print STDERR "Fetching GIs...\n";
    get_gis(\%gi_acc_map);
}

my $gi_count = scalar(keys(%gi_acc_map));
my @gis = ();
my $counter = 0;
foreach my $gi(keys(%gi_acc_map)) {
    $counter++;
    if ($gi_acc_map{$gi} eq '') {
        push(@gis, $gi);
    }
    if (scalar(@gis) == 30 || (scalar(@gis) > 0 && $counter == $gi_count)) {
        untie %gi_acc_map;
        tie %gi_acc_map, 'DB_File', "$dir/.puts.acc.db", O_CREAT|O_RDWR;
        my @accs = get_accs(@gis);
        my $max = scalar(@accs);
        print "Fetched accessions: ";
        for (my $i = 0; $i < $max; $i++) {
            print "$gis[$i]|$accs[$i];";
            $gi_acc_map{$gis[$i]} = $accs[$i];
        }
        print "\n";
        
        @gis = ();
    }
}
untie %gi_acc_map;

sub get_gis {
    my ($hash_ref) = @_;
        
    my $rule = new File::Find::Rule;
    $rule->file();
    $rule->name('*.PUT_member.txt');
    my @member_files = $rule->in($dir);

    foreach my $file(@member_files) {
        open (IN, $file) || die "Failed to open '$file' for reading: $!";
        while (<IN>) {
            chomp;

            my @t = split("\t");

            $hash_ref->{$t[2]} = '';
        }
    }
}
 
sub get_accs {
    my @ids = @_;

    my $retry = 1;
    my @accs;
    while ($retry) {
            my $factory = Bio::DB::EUtilities->new(-eutil => 'efetch',
                         -db => 'nucleotide',
                         -id => \@ids,
                         -rettype => 'acc');


            @accs = split(m{\n},$factory->get_response->content);
            
            if (scalar(@accs) != scalar(@ids)) {
                print STDERR "Couldn't get accessions for all of the GIs, retrying\n";
            } else {
                $retry = 0;
            }
      
    }

    return @accs;
}

END {
    untie %gi_acc_map;
}
