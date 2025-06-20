#!/usr/bin/perl
use strict;
use warnings;

# Ρυθμίσεις αρχικές για τη προσομοίωση ακολουθιών DNA
my $num_seq = shift || 1000;         # Πλήθος τυχαίων ακολουθιών DNA (default: 1000)
my $seq_len = shift || 1_000_000;    # Μήκος κάθε ακολουθίας (default: 1.000.000 βάσεις)
# Έβαλα shift για να επιταχύνω τις δοκιμές όπως (perl erg4.pl 10 10000)

# Ορισμός κατανομών βάσεων
my @katanomes = (
  { name => 'iso',       freqs => { A=>0.25, T=>0.25, G=>0.25, C=>0.25 } }, #"ισοπίθανες" βάσεις
  { name => 'at30_gc70', freqs => { A=>0.15, T=>0.15, G=>0.35, C=>0.35 } }, #Α+Τ=30%  G+C=70%
  { name => 'gc30_at70', freqs => { A=>0.35, T=>0.35, G=>0.15, C=>0.15 } }, #Α+Τ=70%  G+C=30%
);
# Πίνακας κωδικονίων λήξης
my %stop = map { $_=>1 } qw(TAA TAG TGA);

# Συμπληρωματική Αλυσίδα DNA (complement DNA sequences)
my %comp = (A=>'T', T=>'A', G=>'C', C=>'G');

# Συνάρτηση για reverse complement (για εύρεση ORFs και στην άλλη αλυσίδα)
sub revcomp {
    my $seq = shift;
    $seq = reverse $seq;
    $seq =~ tr/ACGT/TGCA/;
    return $seq;
}

# Κύρια διαδικασία ανά κατανομή βάσεων
for my $k (@katanomes) {
    print "\n=== Κατανομή: $k->{name} ===\n";
    my @lengths;  # Λίστα για τα μήκη των ORFs

    # Προσομοίωση DNA ακολουθιών
    for (1..$num_seq) {
        my $dna = '';
        for (1..$seq_len) {
            my $r = rand();
            $dna .= $r < $k->{freqs}->{A} ? 'A'
                  : $r < $k->{freqs}->{A} + $k->{freqs}->{T} ? 'T'
                  : $r < $k->{freqs}->{A} + $k->{freqs}->{T} + $k->{freqs}->{G} ? 'G'
                  : 'C';
        }

        # Εύρεση ORFs και στις δύο κατευθύνσεις
        foreach my $strand ($dna, revcomp($dna)) {
            for my $frame (0,1,2) {  # 3 δυνατά reading frames ανά κατεύθυνση
                for (my $i = $frame; $i < length($strand)-2; $i += 3) {
                    if (substr($strand,$i,3) eq 'ATG') {
                        for (my $j = $i+3; $j < length($strand)-2; $j += 3) {
                            my $codon = substr($strand,$j,3);
                            if ($stop{$codon}) {
                                push @lengths, $j+3-$i;  # Καταγραφή μήκους ORF
                                last;
                            }
                        }
                    }
                }
            }
        }
    }

    # Στατιστικά ανάλυσης
    my $n = scalar @lengths;
    my $sum = 0; $sum += $_ for @lengths;
    my $mean = $n ? $sum/$n : 0;
    my $sq = 0; $sq += ($_-$mean)**2 for @lengths;
    my $var = $n>1 ? $sq/($n-1) : 0;

    printf "Σύνολο ORFs: %d\nΜέσος όρος: %.2f\nΔιασπορά: %.2f\n", $n, $mean, $var;

    # Υπολογισμός για Ιστόγραμμα
    my ($min,$max) = ($lengths[0], $lengths[0]);
    for (@lengths) {
        $min = $_ if $_ < $min;
        $max = $_ if $_ > $max;
    }

    my $bins = 50;

    # Αν όλα τα μήκη είναι ίδια, δείξε μια μπάρα και προχώρα
    if ($max == $min) {
        print "\n Ιστόγραμμα μηκών ORFs: \n";
        printf "%6s - %6s | %6d %s \n", $min, $max, $n, '*' x 50;
        next;
    }

    my $bin_size = ($max-$min)/$bins;
    my @hist = (0) x $bins;
    $hist[int(($_-$min)/$bin_size)]++ for @lengths;

    # Εμφάνιση Ιστογράμματος σε ASCII
    print "\n Ιστόγραμμα μηκών ORFs: \n";
    for my $i (0..$bins-1) {
        my $low  = sprintf("%.0f", $min + $i*$bin_size);
        my $high = sprintf("%.0f", $min + ($i+1)*$bin_size);
        my $bar  = '*' x int(($hist[$i]/$n)*50 + 0.5);
        printf "%6s - %6s | %6d %s\n", $low, $high, $hist[$i], $bar;
    }
}
