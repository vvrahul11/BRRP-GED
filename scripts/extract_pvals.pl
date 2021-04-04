use strict;

my @r;
while (my $l=<>) {
        $l=~s/<\/td>/<\/td>\n/g;
        push(@r,split(/\n/,$l));
}

for (my $i=0; $i<(@r-1); $i++) {
        if ($r[$i+1]=~/Survival\:/) {
                if ($r[$i]=~/top\"\>(.+)\<\/a\>\<\/td\>$/) {
                        #print $1,"\t";
                } else {
                        die($r[$i]);
                }
        }
        if ($r[$i]=~/P value\:/) {
                $r[$i+1]=~s/<td>//;
                $r[$i+1]=~s/\<\/td\>//;
                print $r[$i+1],"\n";
        }
}
