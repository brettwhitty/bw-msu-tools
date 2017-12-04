#!/usr/bin/perl

package PBS::Tool::AAT;

use warnings;
use strict;
use Carp;

## PBS::Client::Job is embedded in the PBS::Client module file
require PBS::Client;

my %global_opts = (
                    queue       =>  'brody_main',
                    name        =>  'AAT_na',
                    nodes       =>  1,
                    ppn         =>  1,
                    wallt       =>  '00:05:00',
                    maillist    =>  'whitty@msu.edu',
                    mailopt     =>  'e',
                  );

sub new {
    my ($class, %args) = @_;

    unless ($args{'input_file'} && $args{'dds_database'}) {
        confess "input_file and dds_database are required by the constructor";
    }
    
    my $root_job = init_commands($args{'input_file'}, $args{'dds_database'}, $args{'gap2_database'});

    return $root_job;
}

sub init_commands {

    my ($input_file, $dds_database, $gap2_database) = @_;
    
    ## default to use the same database for both programs
    unless ($gap2_database) {
        $gap2_database = $dds_database;
    }
    
    ## get a consistent file prefix for use in the output files whether we're dealing with masked or unmasked
    ## (To avoid using a flag I want to use the file extensions for detecting masked aat searches)
    ## if a masked input file was used for dds, then unmasked will be used for subsequent steps
    my $prefix = $input_file;
    if ($input_file =~ /^(.*)\.masked$/) {
        $prefix = $1;
    }
    
    my $dds = { 
                cmd     =>  "dds $input_file $dds_database -f 100 -i 30 -o 75 -p 70 -a 2000",
                stdout  =>  "$prefix.dds.raw",
                opts    =>  {
                                %global_opts,
                                name    =>  'dds',   
                            },
              };
    

    my $ext = { 
                cmd     =>  "ext $prefix.dds.raw",
                stdout  =>  "$prefix.ext.raw",
                opts    =>  {
                                %global_opts,
                                name    =>  'ext',   
                            },
              };
              
    my $extcollapse = {
                cmd     =>  "extCollapse.pl $ext->{stdout}",
                stdout  =>  "$prefix.ext.collapsed",
                opts    =>  {
                                %global_opts,
                                name    =>  'extCollapse',
                            },
                      };

    my $filterext = {
                cmd     =>  "filter_ext $extcollapse->{stdout} -c 10",
                stdout  =>  "$prefix.filter.raw",
                opts    =>  {
                                %global_opts,
                                name    =>  'filter_ext',
                            },
                    };

    my $gap2 = {
                cmd     => "gap2 $prefix $gap2_database $filterext->{stdout} $prefix.gap2.btab",
                stdout  =>  "$prefix.gap2.raw",
                opts    =>  {
                                %global_opts,
                                name    =>  'gap2',   
                            },
               };

    # create jobs           
    my $job_dds = new PBS::Client::Job(
                                        %{$dds->{'opts'}},
                                        cmd     =>  $dds->{'cmd'},
                                        ofile   =>  $dds->{'stdout'},
                                      );
    my $job_ext = new PBS::Client::Job( 
                                        %{$ext->{'opts'}},
                                        cmd     =>  $ext->{'cmd'},
                                        ofile   =>  $ext->{'stdout'},
                                      );
    my $job_extcollapse = new PBS::Client::Job( 
                                        %{$extcollapse->{'opts'}},
                                        cmd     =>  $extcollapse->{'cmd'},
                                        ofile   =>  $extcollapse->{'stdout'},
                                      ); 
    my $job_filterext = new PBS::Client::Job( 
                                        %{$filterext->{'opts'}},
                                        cmd     =>  $filterext->{'cmd'},
                                        ofile   =>  $filterext->{'stdout'},
                                      );
    my $job_gap2 = new PBS::Client::Job( 
                                        %{$gap2->{'opts'}},
                                        cmd     =>  $gap2->{'cmd'},
                                        ofile   =>  $gap2->{'stdout'},
                                      );
                                      
    # set dependencies                                  
    $job_dds->next({ ok => $job_ext});                                      
    $job_ext->next({ ok => $job_extcollapse});                                      
    $job_extcollapse->next({ ok => $job_filterext});                                      
    $job_filterext->next({ ok => $job_gap2});                                      
                                          
    return $job_dds;                                      
}

1;
