#!/bin/bash

#$ -N parse_orthomcl
#$ -cwd
#$ -e /dev/null
#$ -o /dev/null

export PERL5LIB=/opt/rocks/lib/perl5/site_perl/5.10.1/

/run/whitty/potato/v3_2/orthomcl/parse_orthomcl_results.pl -i /projects/potato/hlin/orthoMCL_04202010/potato_n_all_14_plants/potato_n_all_14_plants_all_orthomcl_parsed -F /run/whitty/potato/v3_2/orthomcl/all_14_plants_peps_04202010.fasta -d trees/ -P 2 -p 14_plants -I 0 -o 14_plants.mappings.list

