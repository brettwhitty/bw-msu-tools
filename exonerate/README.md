# Exonerate utility scripts

"Exonerate is a generic tool for pairwise sequence comparison. It allows you to align sequences using a many alignment models, either exhaustive dynamic programming or a variety of heuristics."

See: <https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate>

## Uses

My primary use of Exonerate was in aligning transcript assemblies (eg: PlantGDB PUTs) or polypeptide sequences (eg: UniRef100) to genomic assemblies (eg: Potato), mostly for display in GBrowse as part of the Solanaceae Comparative Genomics Resource web portal.

I developed the Exonerate-based workflow describe below (and supporting scripts) to take advantage of Exonerate's high-quality alignments, optionally allowing for filtering / heuristics to speed up the searching where time / compute resources required.

In my previous annotation work at TIGR/JCVI the tool [AAT](http://aatpackage.sourceforge.net/) was normally used for this purpose, and by design this workflow will produce results of equal or greater quality.

## Exonerate Transcript/Polypeptide to Genome Workflow
For aligning transcript assembly DBs (eg: PlantGDB PUTs) or polypeptide DBs (eg: UniRef100) to genome assemblies and producing GFF3 for display in a genome browser.

### (1) Split long assembly (target) sequence: <br /> *[optional]*

Uses 'fastaoverlap' utility from the Exonerate distribution:

    /exonerate/bin/fastaoverlap
        input_assemblies.fsa 
        > split_assemblies.fsa

...which will default to **100k** *chunks* with **10k** *overlap*.

### (2) Run 'exonerate': <br />(2a) filtered, or (2b) unfiltered

#### (2a) Run the filtered exonerate script:
First index the query database(s) you are going to search with:

    ./bio_db_fasta_indexer.pl /path/to/query_database.fsa
        

Then run the filtering script:

    ./multi_filter_exonerate.pl
        --query          query_database.fsa
        --target         split_assemblies.fsa
        --query_type     e 
        --output_dir     ./output
        --work_dir       ./work_dir 
        --percent_score  50
        --filter_type    [blast|nucmer]
  
Query is the PUTs (or UniRef100), target is the assemblies.
'query_type' flag affects the model used by exonerate, options are:

    e - est2genome
    c - cdna2genome
    t - coding2genome
    p - protein2genome

#### (2b) Run exonerate directly: <br /> *[no heuristic pre-filter]*

Exonerate must be run with the following parameters if not using the filtering script:

    /exonerate/bin/exonerate
        --query         query_database.fsa
        --target        split_assemblies.fsa
        -Q              dna
        -T              dna
        --model         [est2genome|cdna2genome|coding2genome]
        --verbose       0
        --showalignment no
        --showsugar     no 
        --showquerygff  no
        --showtargetgff yes
        --showcigar     yes
        --showvulgar    no
        --ryo           "%ti\t%qi\t%ql\t%qal\t%r\t%s\t%pi\t%ps\n"
        --percent       {SOME_PERCENT_SCORE_CUTOFF}
        --softmasktarget 1 
        --maxintron     {SOME_MAX_INTRON_LENGTH_CUTOFF}
        --subopt        0
        --fsmmemory     1000
        --dpmemory      1000
        > output_file.raw
  
It's important to pick a *percent cutoff* --- which is percent of maximal alignment score possible for the query ([see exonerate docs](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate-user-guide)) --- to filter garbage alignments, and also some reasonable *maximum intron length* parameter.

I've used **--percent** *50*, **--maxintron** *2000*.

The cutoffs chosen will also affect run time.


### (3) Remap the coordinates of the raw output: <br /> *[only if assembly sequence was split in step (1)]*

    ./adjust_fastaoverlap_coords_in_exonerate_raw.pl
        < exonerate_output.raw 
        > remapped_output.raw
  
### (4) Convert raw exonerate output to GFF3:

    ./convert_exonerate_gff_to_gff3.pl
        -i [ output_file.raw | remapped_output.raw ]
        -o output.gff3 
        -s source_name
  
The script has a hash of feature types to skip outputting to GFF (eg: splice5, splice3, intron, utr5, utr3), but you can force these to be added to the output using the '**--no\_skip**' flag.

### (5) Remove overlapping alignments of the same query sequence from the GFF3:

These features may be present because of annotations in overlapping regions of the split assembly sequences, or in the normal output from Exonerate's coding2genome model (which aligns to both strands, and will output the suboptimal alignment on the opposite strand).

    ./filter_overlapping_exonerate_gff3_features.pl
        < output.gff3
        > filtered.gff3

### (6) Done

The '**filtered.gff3**' file created in the last step should now be ready to load to any GFF3-supporting genome browser, or for whatever other purpose it's needed.

### Notes

#### Parallel Execution on Split Assemblies

Splitting the assemblies is an easy way to parallelize the execution on a compute grid. 
