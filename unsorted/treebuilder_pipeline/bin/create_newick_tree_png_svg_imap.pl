#!/usr/bin/perl

# 2 November 2009

# This script is adapted from the fetch_ortholog_diagram.cgi script that
# had been slated to be used directly within the annotation_report page.
# The script will fetch from the bfgr_orthologs database the newick tree 
# data for every ortholog cluster.  The newick file will be converted 
# into a phylogenetic tree and written to a png file.  An image map will 
# be created that will link each leaf node to the annotation report page 
# for the gene at that leaf node.

# For each ortholog cluster, both the png file and image map will be stored
# in the bfgr_orthologs database.

use strict;
use warnings;

use Getopt::Std;
#use HTML::Template;
use GD::SVG;
use File::Temp;
#use DBI;
#use DBI qw(:sql_types);
use Bio::TreeIO;
use XML::Twig;
use File::Basename;
use Cwd;
use Getopt::Long;

my $usage = "\n$0    -t /path/to/tree.newick\n\n";

our ($opt_t, $opt_h, $opt_n);

GetOptions(
    'n:s'   =>  \$opt_n,
    't:s'   =>  \$opt_t,
    'h!'    =>  \$opt_h,
);

if ($opt_h) {
    die $usage;
}

unless (-e $opt_t) {
    die $usage;
}

my $annot = {};
if (defined($opt_n) && -e $opt_n) {
    open (my $infh, '<', $opt_n) || die "$!";
    while (<$infh>) {
        chomp;

        my @t = split("\t", $_);

        grep { $_ =~ s/\%([A-Fa-f0-9]{2})/pack('C', hex($1))/seg; } @t;

        my @annot = ();
        foreach my $col(@t[1 .. $#t]) {
            if ($col ne '') {
                push(@annot, $col);
            }
        }

        $annot->{uc($t[0])} = join("; ", @annot);
    }
}

my $newick_tree = $opt_t;
my $cluster_name = basename($newick_tree, '.final.newick');
#my $database = $opt_d;

#if (!defined($database)) {
#    die $usage;
#}

# This is the width for the entire image.  It allows extra room for the leaf-node labels.
my $IMAGE_WIDTH = 900;

# A hash that is used to find the bfgr database name for each species.
my %db_hash = (
           'Arabidopsis_thaliana' => 'bfgr_3702',
           'Sorghum_bicolor' => 'bfgr_4558',
           'Oryza_sativa' => 'bfgr_4530',
           'Populus_trichocarpa' => 'bfgr_3694',
           'Vitis_vinifera' => 'bfgr_29760',

           'Helianthus_annuus' => 'bfgr_4232',
           'Medicago_sativa' => 'bfgr_3879',
           'Medicago_truncatula' => 'bfgr_3880',
           'Populus_tremula_x_Populus_tremuloides' => 'bfgr_47664',

           'Brachypodium_distachyon' => 'bfgr_15368',
           'Cenchrus_ciliaris' => 'bfgr_35872',
           'Leymus_cinereus_x_Leymus_triticoides' => 'bfgr_407946',
           'Hordeum_vulgare' => 'bfgr_4513',
           'Panicum_virgatum' => 'bfgr_38727',
           'Secale_cereale' => 'bfgr_4550',      
           'Sorghum_propinquum' => 'bfgr_132711',
           'Triticum_aestivum' => 'bfgr_4565',
           'Triticum_turgidum_subsp__durum' => 'bfgr_4567',
           'Triticum_monococcum' => 'bfgr_4568',
           'Zea_mays' => 'bfgr_4577',

           'Picea_glauca' => 'bfgr_3330',
           'Picea_sitchensis' => 'bfgr_3332',
           'Pinus_taeda' => 'bfgr_3352',
          );

# A hash the defines colors for the sequences for each species.
# Acceptable color names are defined in the DrawTree subroutine.
my %db_to_color_hash = (
            # Model genomes
            'bfgr_3702' => 'red',
            'bfgr_4558' => 'blue',
            'bfgr_4530' => 'green',
            'bfgr_3694' => 'brown',
            'bfgr_29760' => 'black',

            # Dicot species.
            'bfgr_4232' => 'dgray',
            'bfgr_3879' => 'dblue',
            'bfgr_3880' => 'orange',
            'bfgr_47664' => 'chocolate',

            # Monocot species.
            'bfgr_15368' => 'dgray',
            'bfgr_35872' => 'dblue',
            'bfgr_407946' => 'orange',
            'bfgr_4513' => 'chocolate',
            'bfgr_38727' => 'orangered',
            'bfgr_4550' => 'dred',
            'bfgr_132711' => 'gold',
            'bfgr_4565' => 'purple',
            'bfgr_4567' => 'cyan',
            'bfgr_4568' => 'lightslategray',
            'bfgr_4577' => 'gray',

            # Gymnosperm species.
            'bfgr_3330' => 'dgray',
            'bfgr_3332' => 'dblue',
            'bfgr_3352' => 'orange',
               );

# Open a connection to the sqlite database.
#my $dbh = DBI->connect("dbi:SQLite:dbname=$database", "", "") || die "\nUnable to open database, $database\n\n";

# Only three different queries will be needed.
#my $get_cluster_and_newick_data_sth = $dbh->prepare(qq/select c.cluster_id, cluster_name, c.newick_tree, c.phyla from cluster c/) || die "Unable to prepare cluster selection statement.";
#my $get_gene_ids_and_db_names_sth = $dbh->prepare(qq/select name, functional_annotation, db_id from cluster_member where cluster_id = ?/) || die "There was a problem preparing the gene and database id select statement.";
#my $load_image_data_sth = $dbh->prepare(qq/update cluster set tree_png = ?, image_map = ? where cluster_id = ?/) || die "\nUnable to prepare update sql statement.\n\n";

# Since we may loop through several clusters (model genome proteins may belong to several clusters),
# the html output will be a bit convoluted.  Print the html after finishing with the processing of
# cluster data.

#my $num_tests = 0;
my $links_to_images = qq|<a name="orthologs_top">\n|;
#$get_cluster_and_newick_data_sth->execute();
#while (my ($cluster_id, $cluster_name, $newick, $phyla) = $get_cluster_and_newick_data_sth->fetchrow_array()) {
#    if (!defined($cluster_id) || !defined($cluster_name) || !defined($newick) || !defined($phyla)) {
    # This should not be a fatal result.  Some genes will not be members of a cluster.
#    die "\nFailed to retrieve complete cluster information from bfgr_ortholog database.\n\n";
#    }

#    $get_gene_ids_and_db_names_sth->execute($cluster_id);
#    my (%cluster_member_to_db_id, %gene_to_func_annotation);
#    while (my ($gene_id, $func_annot, $db_id) = $get_gene_ids_and_db_names_sth->fetchrow_array) { 
#    $cluster_member_to_db_id{$gene_id} = $db_id;
#    $gene_to_func_annotation{$gene_id} = $func_annot;
#    }

    # Write the newick data to a temporary file.
    # Set unlink to 1 to erase 
    # the file when the $tmp_xml_fh variable goes out of scope.
#    my $tmp_newick_fh = File::Temp->new( TEMPLATE => "bfgr_cluster_newick-XXXXXXXX",
#                      SUFFIX => '.newick',
#                      DIR => '/tmp',
#                      UNLINK => 1);
#    my $temp_newick_name = $tmp_newick_fh->filename();
#    print $tmp_newick_fh $newick;

    # It is unclear why the TreeIO method won't read this if it is only 600, but 
    # that's what happens.
#    system("chmod 666 $tmp_newick_fh");

open(my $tree_fh, '<', $newick_tree) || die "$!";

    # Read the newick data into TreeIO object.
    # Extract edge length and parent/child information into hashes.

my ($child_node_to_parent_node, $branch_lengths_ref) = get_hashes_from_tree($newick_tree);

use Data::Dumper;

#print Dumper $child_node_to_parent_node;
#print Dumper $branch_lengths_ref;
#die();

# Prepare some values for the tree drawing function.
my $tree_width = 400;  # Width of the tree portion of the image, not including the labels.
my $image_height = 300;
    
####    if (scalar(keys(%cluster_member_to_db_id)) > 8 ) {
####        $image_height = scalar(keys(%cluster_member_to_db_id)) * 40;
####    }
    my $line_size = 4;
    my $text_border = 3;
    my $border_size = 40;

    # Draw the tree.
    # Write an svg version of the tree.
    # Modifiy the svg xml.
    # Create a png version of the tree.
    # All performed in DrawTree.

    #my $png_file_dir = "/tmp/";
    my $png_file_dir = dirname($newick_tree);
   
    my %gene_to_func_annotation = ();
    #my %db_to_color_hash = ();
    my %cluster_member_to_db_id = ();
    my $phyla = 'Potato';


    my ($image_coords_html, $new_image_height, $png_file) = DrawTree(
        $line_size,                     ## ok 
        $text_border,                   ## ok
        $image_height,                  ## ok
        $tree_width,                    ## ok
        $border_size,                   ## ok
        $child_node_to_parent_node,     ## ok
        $branch_lengths_ref,            ## ok
        $png_file_dir,                  ## ok
        \%cluster_member_to_db_id,      ## ???
        \%gene_to_func_annotation,      ## ???
        \%db_to_color_hash,             ## ???
        $phyla,                         ## ???test 
        $cluster_name                   ## ???test
    );

    # When the annotation_report script grabs these data from the bfgr_orthologs database,
    # it will write them to this file.
    my $png_url = "images/trees/" . $cluster_name . ".png";

    # Write html to describe the imagemap of the tree.
    # Make the width of the image map just a bit wider than the standard width of the image.
    my $div_width = $IMAGE_WIDTH + 20;
    my $image_map;
##    $image_map .= qq|<a name="map_|;
##    $image_map .= $phyla;
##    $image_map .= qq|"></a>\n|;
#    $image_map .= '<div style="height:400px;border:1px solid black;width:'.$div_width.'px;overflow-x:auto;overflow-y:auto;">'."\n";
    $image_map .= qq|    <img src="$png_url" style="border:0" width="$IMAGE_WIDTH" height="$new_image_height" usemap="#map_| . $phyla . qq|"/>\n|;
    $image_map .= qq|    <map name="map_| . $phyla . qq|">\n|;
    $image_map .= $image_coords_html;
    $image_map .= qq|    </map>\n|;
#    $image_map .= qq|</div>\n|;

    # Additional links will be added below each image during the creation of the annotation report page.

    # Load the database with the png file and the image map.
    #my $png_blob = `cat $png_file`;
    #$load_image_data_sth->bind_param(1, $png_blob, SQL_BLOB);
    #$load_image_data_sth->bind_param(2, $image_map);
    #$load_image_data_sth->bind_param(3, $cluster_id);
    #$load_image_data_sth->execute();

    open(my $outfh, '>', $png_file_dir."/".$cluster_name.".inc.html") || die "$!";
    print $outfh $image_map;

    # Be polite to /tmp/ and remove the png file.
#    system("rm $png_file");

#    ++$num_tests;
#    if ($num_tests < 20) {
    # This code tested the ability to read a png from the db and to write the png to file.
#    my $sth = $dbh->prepare("SELECT tree_png FROM cluster WHERE cluster_id = $cluster_id");
#    $sth->execute();
#    my $row = $sth->fetch();
#    my $png_blob_2 = $row->[0];
    # now $png_blob_2 == $png_blob

    # For testing purposes.
#    my $output_file = "/tmp/kchilds/" . $cluster_name . ".png";
#    open OUTF, ">$output_file" || die "\nCan't open $output_file for writing\n";
#    binmode OUTF;
#    print OUTF $png_blob_2;
#    close OUTF || die "Can't close $output_file\n";
#    }
    #
#    if (!($num_tests % 1000)) {
#    print "Processed $num_tests clusters.\n";
#    }
#}

#print "\nFinished $0\n\n";
exit;

sub get_hashes_from_tree {

    my ($tree_file) = @_;

    # The Bio::TreeIO::phyloxml.pm module is not working in the latest Bioperl release.
    #my $treeio = Bio::TreeIO->new(-format => 'phyloxml', -file => $tree_file);
    my $treeio = Bio::TreeIO->new(-format => 'newick', -file => $tree_file);

    my $tree = $treeio->next_tree;

    my $rootnode = $tree->get_root_node;

    # Hashes that will be used by the graphics function.
    my (%branch_lengths, %child_node_to_parent_node);

    # Get the name of the root node.  It's length is 0.
    my @leaf_names;
    foreach my $sub_node ($rootnode->get_all_Descendents()) {
    # Presumably, this is a breadth-first search.
    # This breadth-first ordering is expected when the tree is drawn.
    # So, do not order the @leaf_names array.
    if ($sub_node->is_Leaf()) {
        push @leaf_names, $sub_node->id();
    }
    }
    my $node_name = join ":", @leaf_names;
    $branch_lengths{$node_name} = 0;

    # Examine all descendents of the root node.
    # Determine their names and the length of the branches from
    # each nodes' parent to the node.
    # Also, for each child node, determine it's parent's name,
    # and save the name in the child_node_to_parent_node hash.
    foreach my $node ( $rootnode->get_all_Descendents() ) {
    my $node_name;
    if ( $node->is_Leaf() ) { 
        $node_name = $node->id();
        $branch_lengths{$node_name} = $node->branch_length;
    }
    else {
        # Node is an internal node.
        # Get the names of all of the leaves below this node.
        my @leaf_names;
        foreach my $sub_node ($node->get_all_Descendents()) {
        if ($sub_node->is_Leaf()) {
            push @leaf_names, $sub_node->id();
        }
        }
        $node_name = join ":", @leaf_names;
        $branch_lengths{$node_name} = $node->branch_length;
    }
    my $parent_node = $node->ancestor();
    my @leaf_names;
    foreach my $sub_node ($parent_node->get_all_Descendents()) {
        if ($sub_node->is_Leaf()) {
        push @leaf_names, $sub_node->id();
        }
    }
    my $parent_node_name = join ":", @leaf_names;
    
    $child_node_to_parent_node{$node_name} = $parent_node_name;
    }

    return \%child_node_to_parent_node, \%branch_lengths;
}

sub DrawTree {
    # This subroutine is taken from  Julian Gough's JgoughPhylo perl module.
    # The original code was a hack, and it is only been partially cleaned up.

    # ARGS: DrawTree($outfile,$linesize,$textborder,$image_height,$imagewidth,$bordersize,
    #                 \%nodeup,\%nodex)
    # The imagewidth is the width for purposes of creating the tree.

    my ($linesize, $textborder, $image_height, $imagewidth, $bordersize, $nodeup_ref, $nodex_ref, $png_file_dir, $cluster_member_to_db_id_ref, $gene_to_func_annotation_ref, $db_to_color_hash_ref, $phyla, $cluster_id) = @_;

    # The %nodeup hash provides the names of the parent node for each child node.
    my %nodeup = %{$nodeup_ref};

    # The %nodex hash provides the lengths of each tree branch from parent to child.
    # These values will be used for determining the x coordinate for each node.
    my %nodex = %{$nodex_ref};

    my @temp;
    my @members;
    my $node;
    # The %nodey hash will provide the y coordinate for each node and line.
    # These coordinates are determined based on the number of leaf nodes in the tree
    # and on the total desired height of the tree.
    my %nodey;
    my (%labely, %ymin, %ymax);  
    my ($i, $j, $x, $y);
    my $max = 0; 
    my $min = 999999;
    my ($black, $image, $white);
    my $labelheight = 0;
    my $labelwidth = 0; 
    my (@len, @order, @nodes);
    my ($line, $lineup);
    my $databoxheight = 0;
    my $databoxwidth = 0;
    my $extraheight;
    my %yminno;

    # Can databoxheight and labelheight have static assignments???
    foreach $node (keys(%nodex)) {
    $i = 0;
    if ($databoxheight < $i * 7 + 2 * $textborder) {
        $databoxheight = $i * 7 + 2 * $textborder;
    }
    }
    if ($labelheight - $databoxheight < $linesize) {
    $labelheight = $databoxheight + $linesize;
    }

    # The @order_node_names array contains ordered %nodex keys.
    # This array orders child nodes before parent nodes which
    # have names that are composed of a concatenation of their child node names.
    # So, parent node names are necessarily longer than child node names.
    my @ordered_node_names = sort sort_by_length keys(%nodex);
    my @ordered_node_names_rev = reverse @ordered_node_names;

    # X-COORDINATES
    # Find-coordintes by adding up node lengths.
    # Start at the root node (which has the longest name).
    foreach $node (@ordered_node_names_rev) {

    # For our purposes, all branch lengths must be positive values.
    if ($nodex{$node} < 0) {
        $nodex{$node} = -$nodex{$node};
    }

    # The x coordinate of a node is the sum of the branch lengths
    # of the node and all of its ancestors.
    if (exists($nodeup{$node})) {
        $nodex{$node} = $nodex{$node} + $nodex{$nodeup{$node}};
    }
    }

    # Find the length of the root node and the max and min node distances within the entire tree.
    my $root = '';
    foreach $node (keys(%nodex)) {
    if (length($node) > length($root)) {
        $root = $node;
    }
    if ($max < $nodex{$node}) {
        $max = $nodex{$node};
    }
    if ($min > $nodex{$node} or $min == 999999) {
        $min = $nodex{$node};
    }
    }

    @members = split /:/, $root;
    my $number_leaves = scalar(@members);

    # The image size is based on the number of leaf nodes and various extra buffers.
    $extraheight = $image_height - (($databoxheight - 1) * $number_leaves + $labelheight * $number_leaves + 2 * $bordersize) - 40;
    if ($extraheight < 0) {
    $image_height = $databoxheight * ($number_leaves - 1) + $labelheight * $number_leaves + 2 * $bordersize;
    $extraheight = 0;
    }
    if ($imagewidth < 2 * $bordersize + $linesize + $labelwidth + $databoxwidth) {
    $imagewidth = 2 * $bordersize + $linesize + $labelwidth + $databoxwidth;
    }

    my $max_min_diff = $max - $min;
    if ($max_min_diff == 0) {
    print "cluster\t$cluster_id\t$max\t$min\n";
    $max_min_diff = 1;
    }

    # Now, adjust the x coordinates so that they take into account borders, total image width
    # and line size.
    foreach $node (keys(%nodex)) {
    $nodex{$node} = $bordersize + $nodex{$node} * 
      (($imagewidth - ($bordersize * 2) - $databoxwidth - $labelwidth) /
       ($max_min_diff)) + $linesize / 2;
    }

    # Find the new max and min x coordinates.
    my $new_max = 0;
    my $new_min = 999999999;
    foreach $node (keys(%nodex)) {
    if ($nodex{$node} > $new_max) {
        $new_max = $nodex{$node};
    }
    if ($nodex{$node} < $new_min) {
        $new_min = $nodex{$node};
    }
    }

    # Determine the scale for the branches of the tree.
    my $distance = $max_min_diff;
    my $new_distance = $new_max - $new_min;
    my $scaling_factor = $new_distance / $distance;
    my $length_legend_bar = 10 * $scaling_factor;
    my $legend_bar_start_x = $bordersize + $linesize / 2;
    my $legend_bar_end_x = $legend_bar_start_x + $length_legend_bar;
    my $legend_bar_y = $image_height - 10;
    $y = $bordersize - $labelheight / 2 - $databoxheight - $extraheight / ($number_leaves - 1) - 1;

    # y coordinates for leaf children.
    foreach $node (@members) {
    $y = $y + $databoxheight + $labelheight + $extraheight / ($number_leaves - 1);
    $nodey{$node} = $y;
    $labely{$node} = $y;
    if (!exists($ymax{$nodeup{$node}})) {
        $ymax{$nodeup{$node}} = $y;
    }
    elsif ($ymax{$nodeup{$node}} < $y) {
        $ymax{$nodeup{$node}} = $y; 
    }
    if (!exists($ymin{$nodeup{$node}})) {
        $ymin{$nodeup{$node}} = $y;
        $yminno{$nodeup{$node}} = $node;
    }
    elsif ( $ymin{$nodeup{$node}} > $y) {
        $ymin{$nodeup{$node}} = $y;
        $yminno{$nodeup{$node}} = $node;
    }
    
    }

    # Scale all of the y coordinates to account for the desired image height.
    # And calculate y coordinates for the non-leaf nodes.
    foreach $node (@ordered_node_names) {
    if (!exists($nodey{$node})) {

        $nodey{$node} = ($ymax{$node} - $ymin{$node}) / 2 + $ymin{$node};
        my @child_nodes = split /:/, $yminno{$node};
        $min = 0;
        foreach my $child_node (@child_nodes) {
        if ($min < $nodey{$child_node}) {
            $min = $nodey{$child_node};
        }
        }
        $labely{$node} = $min + ($databoxheight + $labelheight) / 2 + $extraheight / (2 * ($number_leaves - 1));
        if (exists($nodeup{$node})) {
        unless (exists($ymax{$nodeup{$node}})) {
            $ymax{$nodeup{$node}} = $nodey{$node};
        }
        elsif ($ymax{$nodeup{$node}} < $nodey{$node}) {
            $ymax{$nodeup{$node}} = $nodey{$node};
        }
        unless (exists($ymin{$nodeup{$node}})) {
            $ymin{$nodeup{$node}} = $nodey{$node};
            $yminno{$nodeup{$node}} = $node;
        }
        elsif ( $ymin{$nodeup{$node}} > $nodey{$node}) {
            $ymin{$nodeup{$node}} = $nodey{$node};
            $yminno{$nodeup{$node}} = $node;
        }
        }
    }
    }

    # Prepare a name for the svg file.
    my $svg_name = $png_file_dir . "/" . $cluster_id . ".svg";
    my $final_png_file = $png_file_dir . "/" . $cluster_id . ".png";

    # Make the image.
    # $IMAGE_WIDTH is a width that takes into account the expected length of long node names.
    $image = new GD::SVG::Image($IMAGE_WIDTH, $image_height);

    # Allocate some colours
    $white = $image->colorAllocate(255,255,255);
    $black = $image->colorAllocate(0,0,0);

    my %color_codes = (
               # The GD library may accept stock color names, but GD::SVG does not.
               'black' => $image->colorAllocate(0,0,0),
               'red' => $image->colorAllocate(255,0,0),
               'blue' => $image->colorAllocate(0,0,255),
               'cyan' => $image->colorAllocate(0,255,255),
               'subcyan' => $image->colorAllocate(0,183,235),
               'darkcyan' => $image->colorAllocate(0,139,139),
               'seagreen' => $image->colorAllocate(32,178,170),
               'green' => $image->colorAllocate(0,128,0),
               'yellowgreen' => $image->colorAllocate(64,128,0),
               'brown' => $image->colorAllocate(165,42,42),
               'dgray' => $image->colorAllocate(49,49,49),
               'dblue' => $image->colorAllocate(0,0,139),
               'orange' => $image->colorAllocate(225,135,0),
               'chocolate' => $image->colorAllocate(210,105,30),
               'orangered' => $image->colorAllocate(255,69,0),
               'dred' => $image->colorAllocate(139,0,0),
               'gold' => $image->colorAllocate(235,195,0),
               'purple' => $image->colorAllocate(128,0,128),
               'lightslategray' => $image->colorAllocate(99,106,133),
               'gray' => $image->colorAllocate(88,88,88),
              );

    # Make the background transparent and interlaced
    $image->transparent($white);
    $image->interlaced('true');

    # Set the thickness of all lines.
    $image->setThickness(4);

    # Draw horizontal and vertical lines of the tree.
    foreach $node (@ordered_node_names) {
    $line = $linesize;
    if (exists($nodeup{$node})) {
        $lineup = $linesize;
        # Create horizontal lines.
        $image->line($nodex{$nodeup{$node}}, $nodey{$node}, $nodex{$node}, $nodey{$node}, $black);
    }
    if (exists($ymax{$node})) {
        # Create vertical lines.
        $image->line($nodex{$node}, $ymin{$node}, $nodex{$node}, $ymax{$node}, $black);
    }
    }

    # The node identifiers are created here.
    # The stringFT() is not implemented.  The simple string() function
    # must be used.
    # Also, create the html with coordinates for the image map.
    # <area shape='rect' coords='0,0,82,126' href='sun.htm' alt='Sun' />
    my $image_coords_html;
    my $offset_from_line = 4;
    foreach $node (@members) {
    @temp = split /\n/, $node;
    $i = 0;
    foreach $line (@temp) {
        $y = ($labely{$node} - 12 * scalar(@temp) / 2 + 12 * $i - 2);
        $x = $nodex{$node} + $textborder + $offset_from_line;
        # This is the largest font size possible, but font-size can be manipulated
        # in the .svg xml file.
        # In the .svg file, this is Helvetica, 15 point font.
        # Safari does not correctly interpret this as a Helvetica 15 point font.

        #######################################################################
        # The color of the leaf-node label should be determined by the species.
        my $string_color = $color_codes{'black'};
        if ($line =~ /^PGSC/i) {
            $string_color = $color_codes{'blue'};
        } elsif ($line =~ /^At/i) {
            $string_color = $color_codes{'green'};
        } elsif ($line =~ /^GSVIV/i) {
            $string_color = $color_codes{'purple'};
        } elsif ($line =~ /^Cucsa/i) {
            $string_color = $color_codes{'yellowgreen'};
        } elsif ($line =~ /^\d+/i) {
            $string_color = $color_codes{'subcyan'};
        } elsif ($line =~ /^evm.TU/i) {
            $string_color = $color_codes{'darkcyan'};
        } elsif ($line =~ /^POPTR/i) {
            $string_color = $color_codes{'seagreen'};
        } elsif ($line =~ /^LOC_Os/i) {
            $string_color = $color_codes{'brown'};
        } elsif ($line =~ /^Sb/i) {
            $string_color = $color_codes{'chocolate'};
        } elsif ($line =~ /^GR|^AC\d+|^AF\d+|^AY\d+|^EF\d+/i) {
            $string_color = $color_codes{'gold'};
        } elsif ($line =~ /^Bradi/i) {
            $string_color = $color_codes{'orangered'};
        } elsif ($line =~ /^jgi\|Selmo/i) {
            $string_color = $color_codes{'lightslategray'};
        } elsif ($line =~ /^jgi\|Phypa/i) {
            $string_color = $color_codes{'dgray'};
        } elsif ($line =~ /^Au9/i) {
            $string_color = $color_codes{'gray'};
        }

        $image->string(
            gdGiantFont, 
            $x, 
            $y, 
            $line, 
            $string_color, ##$color_codes{$db_to_color_hash_ref->{$cluster_member_to_db_id_ref->{$line}}}
        );
        $i++;

        # Round the coordinates down.
        $x = sprintf("%u", $x);
        $y = sprintf("%u", $y);

        # Trying to estimate the point-based length of the string.
        my $new_x = $x + (length($line) * 8);  
        # Trying to estimate the point-based height of the string.
        my $new_y = $y + 14;
        my $coords = join ",", ($x, $y, $new_x, $new_y);
        my $link = "#"; ## CHANGE THIS LINK
        $link = get_link($line);
        my $alt = $line;
        my $title = $line;
        if (defined($annot->{uc($line)})) {
            $title = $annot->{uc($line)};
        }
        $image_coords_html .= "        <area shape='rect' coords='$coords' target='phylo_link' href='$link' class='treeImageTip' title='$title' alt='$alt' />\n";
    }
    }

    # Put a title on this image.
    $image->string(
        gdGiantFont, 
        40, 
        8, 
        "Gene tree for $cluster_id", 
        $black
    );

    # Print the SVG file.
    open OUTSVG, ">$svg_name";
    print OUTSVG $image->svg;
    close OUTSVG;

    # Use XML::Twig to modify the svg file.
    my $twig = new XML::Twig( twig_handlers => { text => \&text_handler } );
    $twig->parsefile($svg_name);
    # Write changes to the file.
    open OUTSVG, ">$svg_name";
    print OUTSVG $twig->sprint();
    close OUTSVG;

    # Convert the file to .png format using batik-rasterizer.jar
    #my $rasta = "java -Djava.awt.headless=true -jar /share/apps/batik/batik-rasterizer.jar $temp_name >& /dev/null";
    my $rasta = "java -Djava.awt.headless=true -jar /share/apps/batik-1.7/batik-rasterizer.jar $svg_name >& /dev/null";
    system($rasta);

    # Clean up the svg file.
#    system("rm $temp_name");

    return ($image_coords_html, $image_height, $final_png_file);
}

sub sort_by_length {
    length($a) <=> length($b);
}


sub text_handler {
    my ($twig, $text) = @_;

    $text->set_att('font-family', 'Arial,Verdana');
    $text->del_att('font');
    $text->set_att('font-weight', 'normal');

    my $content = $text->text();
    if ($content =~ /^Phylogeny of cluster/) {
    # This is the title of the image.
    $text->set_att('font-size', '18');
    }
}

sub get_link {
    my ($acc) = @_;

    if ($acc =~ /^PGSC/) {
        return 'http://potatogenomics.plantbiology.msu.edu/cgi-bin/gbrowse/potato/?name=Sequence:'.$acc;
    } elsif ($acc =~ /^At/i) {
        return 'http://gbrowse.arabidopsis.org/cgi-bin/gbrowse/arabidopsis/?name=Sequence:'.$acc;
    } elsif ($acc =~ /^POPTR/i) {
        return 'http://solanaceae.plantbiology.msu.edu/cgi-bin/gbrowse/poplar/?name=Sequence:'.$acc;
    } elsif ($acc =~ /^LOC_Os/i) {
        return 'http://rice.plantbiology.msu.edu/cgi-bin/gbrowse/rice/?name=Sequence:'.$acc;
    } elsif ($acc =~ /^GSVIV/i) {
        return 'http://solanaceae.plantbiology.msu.edu/cgi-bin/gbrowse/grape/?name=Sequence:'.$acc;
    } else {
        return '#';
    }
}
