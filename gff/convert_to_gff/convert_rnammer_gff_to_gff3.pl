#!/opt/rocks/bin/perl

my $counter = 0;

while (<>) { 
    if (/^#/) { 
        print;
        next;
    }
    chomp;
    @t = split("\t");
    $t[8] = "ID=".$t[0]."-rnammer".++$counter.";Name=".$t[8];
    print join("\t", @t)."\n";
}
