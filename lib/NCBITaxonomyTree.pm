#!/usr/bin/perl

package NCBITaxonomyTree;

### 
### Test using Bio::Taxon to look up taxonomy IDs by species name using the NCBI taxonomy dumps
###
### Brett Whitty
### whitty@msu.edu
###

use strict;
use warnings;

use Carp;
use Bio::Taxon;
use File::Basename qw{ dirname };
use Cwd qw{ abs_path }; 


my $no_cache = 0;

sub new {
            my ($class, %args) = @_;

            ## argument 'no_cache' will prevent the object
            ## from caching the results of lookups for
            ## name_has_ancestor_name
            ## id_has_ancestor_name
            ## in case less memory usage is preferable to faster run time
            if (defined($args{'no_cache'})) {
                $no_cache = 1;                
            }

            if (! defined($args{'path'})) {
               
                ## by default look for an NCBI taxonomy dump directory where the module is located
                my $self_path = abs_path(dirname($INC{'NCBITaxonomyTree.pm'}));
                my $dump_path = (-e "$self_path/taxonomy_dump") ? "$self_path/taxonomy_dump" : '';
                
                if ($dump_path) {
                    $args{'path'} = $dump_path;
                } else {
                    confess "Must provide path to NCBI taxonomy dump files with constructor 'path' argument";
                }
            }
            
            _init_taxon_db($args{'path'});

            return bless [], $class;
}    

{
    my $taxon_db;
    my $cache = {};
                               
    sub get_children_species {
        #my $pot = $taxon_db->get_taxon(-name => 'Solanum tuberosum');
        #my @pot_c = $taxon_db->each_Descendent($pot);
        ## iterate through children and return array of species level tax ids

        croak "Not implemented yet";
    }

    sub id_has_ancestor_name {
        my ($self, $node_id, $name) = @_;

        my $node = $taxon_db->get_taxon(-taxonid => $node_id) || confess "Taxon id '$node_id' not found in DB";
        
        my $cache_key = "ihan|$node_id|$name";
        
        return (defined($cache->{$cache_key})) ?
              $cache->{$cache_key} 
            : _cache_result($cache_key, has_ancestor_name($self, $node, $name));
    }
    
    sub name_has_ancestor_name {
        my ($self, $node_name, $name) = @_;

        my $node = $taxon_db->get_taxon(-name => $node_name) || confess "Taxon name '$node_name' not found in DB";
        
        my $cache_key = "nhan|$node_name|$name";

        #use Data::Dumper;
        #print Dumper $cache;

        return (defined($cache->{$cache_key})) ?
               $cache->{$cache_key} 
            :  _cache_result($cache_key, has_ancestor_name($self, $node, $name));
    }
    
    
    sub has_ancestor_name {
        ## option for $node->rank cutoff?
        my ($self, $node, $name) = @_;
    
        my $ancestor = $taxon_db->get_taxon(-name => $name);

        confess "Taxon name '$name' not found in DB" unless defined $ancestor;

        _has_ancestor($node, $ancestor->id);
    } 

    sub has_ancestor_id {
        my ($self, $node, $id) = @_;
    
        my $ancestor = $taxon_db->get_taxon(-taxonid => $id);

        confess "Taxon id '$id' not found in DB" unless defined $ancestor;    

        _has_ancestor($node, $ancestor->id);
    }

    sub get_taxon_id {
        my ($self, $name) = @_;
       
        my $taxon = $taxon_db->get_taxon(-name => $name); 
         
        carp "Failed to find taxonomy node for '$name' in taxonomy db dump" unless defined $taxon;
        
        return $taxon->id;
    }
 
    sub valid_scientific_name {
        my ($self, $name) = @_;
        
        if ($taxon_db->get_taxon(-name => $name)) {
           return 1;
        } else {
           return 0;
       }
    }
   
    sub get_scientific_name {
        my ($self, $taxon_id) = @_;
        
         my $taxon = $taxon_db->get_taxon(-taxonid => $taxon_id);
         
         return $taxon->scientific_name;
    }
     
    sub get_common_name {
        my ($self, $taxon_id) = @_;
        
         my $taxon = $taxon_db->get_taxon(-taxonid => $taxon_id);
         
         return $taxon->common_name;
    }

    sub is_valid_taxon_id {
        my ($self, $taxon_id) = @_;

        my $taxon = $taxon_db->get_taxon(-taxonid => $taxon_id);

        return (defined($taxon)) ? 1 : 0;
    }

    sub get_taxon_id_for_rank_by_taxon_id {
        my ($self, $node_id, $rank) = @_;

        my $node = $taxon_db->get_taxon(-taxonid => $node_id) || confess "Taxon id '$node_id' not found in DB";

        my $taxon_id = undef;
        while ($node = $node->ancestor) {
            if ($node->rank =~ /^$rank$/i) {
                return $node->id;
            }
        }

        return undef;
    }

    sub get_name_for_rank_by_taxon_id {
        my ($self, $node_id, $rank) = @_;

        my $taxon_id = $self->get_taxon_id_for_rank_by_taxon_id($node_id, $rank);

        return $self->get_scientific_name($taxon_id);
    }

    sub _has_ancestor {
        my ($node, $target_id) = @_;

        my $ancestor = $node->ancestor || return 0;

        if ($node->id eq $target_id) {
            return 1;
        } else {
            _has_ancestor($node->ancestor, $target_id);
        }
    }

    sub _init_taxon_db {
        my ($path) = @_;

        $taxon_db = new Bio::DB::Taxonomy(                                
                                -source        => 'flatfile',
                                -directory     =>  $path.'/index',
                                -nodesfile     =>  $path.'/nodes.dmp',
                                -namesfile     =>  $path.'/names.dmp',
                                         );
    
        return $taxon_db;
    }

    sub _cache_result {
        my ($key, $result) = @_; 
       
        unless ($no_cache) { 
            $cache->{$key} = $result;
        }

        return $result;
    }
}

1;
