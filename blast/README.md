## BLAST Utility Scripts

### blast\_output\_to\_simple_table.pl

Uses Bio::SearchIO to output a very simple table:
	
	query_name	hit_name	rank	score	sig
	
from raw BLAST output.

	./blast_output_to_simple_table.pl (infile.raw|STDIN) [e_value_cutoff|default=1e-5]

### blast\_output\_to\_table.pl

Like simple table but adds 'identity', 'similarity' and 'coverage'; adds filtering on identity and coverage.

### blast\_table\_interval\_collapse.pl

"Takes NCBI BLAST tabular output (-m 8) and generates collapsed intervals with the specified number of tiers in an output format similar to filter.pl from the AAT component"

There's a lot of repurposed code in here from scripts in TIGR's ["Analysis and Annotation Tool" (AAT)](http://aatpackage.sourceforge.net/).

### do\_paralog\_blastp.pl

Uses BioPerl's Bio::Tools::Run::StandAloneBlast to do a BLASTP for the purpose of indentifying putative paralogs in a subject DB.

### do\_rbh\_blastp.pl

Uses BioPerl's Bio::Tools::Run::StandAloneBlast to run BLASTP and parse out reciprocal best hits (RBH) as a rough approximation of putative orthologs between query and subject DBs.
