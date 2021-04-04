use strict;

my %counts;
while (my $l=<>) {
	chomp $l;

	@_=split(/\t/,$l);

	++$counts{$_[0]};
	++$counts{$_[1]};
}


foreach my $gene (sort {$counts{$b}<=>$counts{$a}} keys %counts) {
	print join("\t",$gene,$counts{$gene}),"\n";
}
