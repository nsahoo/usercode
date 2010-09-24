#! /usr/bin/perl

$arg = "";
for($i = 0; $i < 26; $i++) {
    if (-e "zee_$i.root") {
        $arg = $arg."zee_$i.root ";
    }
}
$system = `hadd new.root $arg`;
print "\n";
