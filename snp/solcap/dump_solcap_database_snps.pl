#!/opt/rocks/bin/perl

use Carp;
use DBI;
use Getopt::Long;
use Bio::DB::Fasta;

my ($db_user, $db_pass);

GetOptions(
    'user|u=s'      =>  \$db_user,
    'password|p=s'  =>  \$db_pass,
);

if (! $db_user || ! $db_pass) {
    croak "Must provide --user and --password";
}

my $db_name     = 'solcap_snp_paper';
my $db_server   = 'pg-dev';

my $dbh = DBI->connect("dbi:Pg:dbname=$db_name;host=$db_server", $db_user, $db_pass);

my $ref_orgs = {
    'potato'    =>  {
        'dm'    =>'/projects/solcap/PGSC0003DMS.fa',
    },
    'tomato'    =>  {
        'sbm'   =>'/projects/solcap/sanger_snp_analysis/microtom_sbm_snp/tomato_sbm_contigs_031610.fasta',
        'tgi'   =>'/projects/solcap/sanger_snp_analysis/ta496_tgi_snp/tomato_tgi_genome.fasta',
    },
};

my $snp_orgs = {
    'potato'    =>  [
        'atlantic',
        'bintje',
        'kennebec',
        'premier',
        'shepody',
        'snowden',
    ],
    'tomato'    =>  [
        'fl7600',
        'microtom',
        'nc84173',
        'oh086405',
        'oh9242',
        'pi114490',
        'pi128216',
        'ta496',
    ],
};


foreach my $ref(keys(%{$ref_orgs})) {
    foreach my $ref_type(keys %{$ref_orgs->{$ref}}) {

        my $db = new Bio::DB::Fasta($ref_orgs->{$ref}->{$ref_type});
        
        foreach my $snp_type(@{$snp_orgs->{$ref}}) {

            my $table_name = $ref."_".$snp_type."_".$ref_type."_snp";

#            ## check to see if the table has a genome_ori column
#            my $table_type_query = "select * from $table_name limit 1";
#            my $tqh = $dbh->prepare($table_type_query);
#            $tqh->execute();
            
#            my $thref = $tqh->fetchrow_hashref();
#            my $has_genome_ori = 0;
#            if (defined($thref->{'genome_ori'})) {
#                $has_genome_ori = 1;
#            }
            
            my $query = '';

#            if ($has_genome_ori) {
#                $query = "select ref_seq, ref_pos, ref_base, snp_base, genome_ori from $table_name order by ref_seq, ref_pos";
#            } else {
                $query = "select ref_seq, ref_pos, ref_base, snp_base, '1' as genome_ori from $table_name order by ref_seq, ref_pos";
#            }

            #print $query."\n";    
            my $sth = $dbh->prepare($query);
            $sth->execute();
            while (my $row = $sth->fetchrow_hashref()) {

                my $snp = uc($db->seq($row->{'ref_seq'}, $row->{'ref_pos'}, $row->{'ref_pos'}));

                if ($snp ne uc($row->{'ref_base'}) && $snp eq reverse_complement_dna($row->{'ref_base'})) {
                    $row->{'genome_ori'} = -1;
                }

                if ($row->{'genome_ori'} == -1) {
                    $row->{'ref_base'} = reverse_complement_dna($row->{'ref_base'});
                    $row->{'snp_base'} = reverse_complement_dna($row->{'snp_base'});
                }
                print join("\t", (
                    $ref,
                    $ref_type,
                    $snp_type,
                    $row->{'ref_seq'},
                    $row->{'ref_pos'},
                    $row->{'ref_base'},
                    $row->{'snp_base'},
                ))."\n";
            }
        }
    }
}


sub reverse_complement_dna {
    my ($r_seq) = @_;

    $r_seq =~ tr/AaCcGgTtMmRrWwSsYyKkVvHhDdBb/TtGgCcAaKkYyWwSsRrMmBbDdHhVv/;
    $r_seq = reverse($r_seq);

    return $r_seq;
}
