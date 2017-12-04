Exonerate Workflow
==================
For aligning transcript sequences (eg: PlantGDB PUTs) to genome assemblies and producing GFF3 for display in a genome browser.

### Split long assembly (target) sequence: [optional]

    fastaoverlap input_assemblies.fsa >split_assemblies.fsa

will default to 100k chunks with 10k overlap.


### Run 'exonerate' filtered / unfiltered

#### Run the filtered exonerate script:
First index the query database(s) you are going to search with:

    bio_db_fasta_indexer.pl /path/to/puts_database.fsa 
    
Then run the script:

    multi_filter_exonerate.pl -q puts_database.fsa -t split_assemblies.fsa
        --query_type e --output_dir ./output --work_dir ./work_dir 
        --percent_score 50 --filter_type [blast|nucmer]
  
Query is the PUTs, target is the assemblies. query_type flag affects the model used by exonerate, options are:

    e - est2genome, c - cdna2genome, t - coding2genome, p - protein2genome

#### Run exonerate directly:

Exonerate must be run with the following parameters if not using the filtering script:

    exonerate -q puts_database.fsa -t split_assemblies.fsa -Q dna -T dna
        --model [est2genome|cdna2genome|coding2genome] --verbose 0 --showalignment no
        --showsugar no --showquerygff no --showtargetgff yes --showcigar yes --showvulgar no
        --ryo "%ti\t%qi\t%ql\t%qal\t%r\t%s\t%pi\t%ps\n" --percent SOME_PERCENT_SCORE_CUTOFF
        --softmasktarget 1 --maxintron SOME_MAX_INTRON_LENGTH_CUTOFF
        --subopt 0 --fsmmemory 1000 --dpmemory 1000 > output_file.raw
  
It's important to pick a percent cutoff (which is percent of maximal alignment score possible for the query --- see exonerate docs) to filter garbage alignments, and also some reasonable maximum intron length parameter. (I've used --percent 50 and --maxintron 2000) These also affect run time.


### Remap the coordinates of the raw output if the sequence was split in step (1): [optional]

    adjust_fastaoverlap_coords_in_exonerate_raw.pl <exonerate_output.raw >remapped_output.raw
  

### Convert raw exonerate output to GFF3:

    convert_exonerate_gff_to_gff3.pl -i [output_file.raw|remapped_output.raw]
        -o output.gff3 -s source_name
  
The script has a hash of feature types to skip outputting to GFF (eg: splice5, splice3, intron, utr5, utr3), but you can force these to be added to the output using the --no_skip flag.
 

### Remove overlapping alignments of the same query sequence from the GFF3:
(in the overlapping regions of the split assemblies, or produced by coding2genome)

    filter_overlapping_exonerate_gff3_features.pl <output.gff3 >filtered.gff3
