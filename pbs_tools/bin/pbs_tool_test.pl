#!/usr/bin/perl

use strict;
use warnings;
use Carp;

use lib '/mnt/lustre/whitty/perl5lib/lib/perl5/site_perl/5.8.8';
use lib '/mnt/home/whitty/lib';

use Data::GUID;
use PBS::Client; 
use PBS::Tool::AAT; 

my $guid = new Data::GUID();

my $client= new PBS::Client; 

my $aat = new PBS::Tool::AAT(
                                input_file      =>  '/mnt/lustre/whitty/projects/test/fasta_repository/AC206935.fsa.masked',
                                dds_database    =>  '/mnt/lustre/whitty/projects/test/fasta_repository/Solanum_lycopersicum.mRNA.PUT.fasta',
                            ); 

$client->qsub($aat);
