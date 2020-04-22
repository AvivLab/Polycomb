use strict;

my $counter=1;
while (my $line=<STDIN>) {
    chomp $line;
    my $folder= "2020_"."$counter\n";
    print $line." $folder";
    $counter = $counter + 1;
    
}
