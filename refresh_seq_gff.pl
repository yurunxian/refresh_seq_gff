################
#  Runxian Yu  #
#   2025/2/8   #
################

use strict;
use warnings;
use Modern::Perl;
use Getopt::Long;

my $original_fasta = "";
my $original_gff = "";
my $modification = "";
my $modified_fasta = "";
my $modified_gff = "";
my $warnings = "";
my $help = 0;

GetOptions ("f=s" => \$original_fasta,
            "g=s" => \$original_gff,
            "m=s" => \$modification,
            "fm=s" => \$modified_fasta,
            "gm=s" => \$modified_gff,
            "o=s" => \$warnings,
            "h" => \$help);
if ($original_fasta eq "" or 
    $original_gff eq "" or
    $modification eq "" or
    $modified_fasta eq "" or
    $modified_gff eq "" or
    $warnings eq "" or $help eq 1
    ) {printHelp(); exit;}

open (FASTA,"<$original_fasta") or die("Cannot find <original fasta>.");
open (GFF,"<$original_gff") or die("Cannot find <original gff>.");
open (MODIFICATION,"<$modification") or die("Cannot find <modification file>.");
open (FASTAO,">$modified_fasta") or die("Cannot open <modified sequence>.");
open (GFFO,">$modified_gff") or die("Cannot open <modified gff>.");
open (ERROR,">$warnings") or die("Cannot open <summary>.");

my $head = "";
my $count = 0;
my $count_gene = 0;
my @changed_feature = ();
my @original_gff = <GFF>;
my %all_site = ();
my %original_gff = ();
my %modify = ();
my %fasta = ();

say "Loading genome sequence......";
while (<FASTA>) {
    $_ =~ s/[\r\n]+//gi;
    if (/>(\S+)/) {$head = $1}
    else {$fasta{$head} .= $_}
}

say "Loading modification info......";
while (<MODIFICATION>) {
    $_ =~ s/[\r\n]+//gi;
    my @tmp = split /\t/,$_;
    if (scalar @tmp eq 0) {next}
    my $change = length($tmp[3])-length($tmp[2]);
    $all_site{$tmp[0]}{$tmp[1]}{"modified"} = $change;
    $modify{$tmp[0]}{$tmp[1]-1} = $tmp[2]."-".$tmp[3];
}

say "Loading gff file......";
foreach my $line (@original_gff) {
    $line =~ s/[\r\n]+//gi;
    if ($line =~ /#/) {next}
    my @tmp = split /\t/,$line;
    if (scalar @tmp eq 0) {next}
    $all_site{$tmp[0]}{$tmp[3]}{"boundary"} = $tmp[3];
    $all_site{$tmp[0]}{$tmp[4]}{"boundary"} = $tmp[4];
}

say "Processing......";
foreach my $chr (sort keys %all_site) {
    my $trace = 0;
    foreach my $i (sort {$a <=> $b} keys %{$all_site{$chr}}) {
        if (exists $all_site{$chr}{$i}{"boundary"}) {
            $all_site{$chr}{$i}{"boundary"} += $trace;
        }
        if (exists $all_site{$chr}{$i}{"modified"}) {
            $trace += $all_site{$chr}{$i}{"modified"};
        }
    }
    foreach my $i (sort {$b <=> $a} keys %{$modify{$chr}}) {
        my @change = split /-/,$modify{$chr}{$i};
        my $seq = substr($fasta{$chr},$i,length($change[0]));
        if ($seq eq $change[0]) {
            substr($fasta{$chr},$i,length($change[0])) = $change[1];
            $count++;
        }
        else {my $pos = $i+1; say "Warning: Sequence incongruence at $chr:$pos. \"$change[0]\" supposed but \"$seq\" found."}
    }
}

say "Preparing output......";
foreach my $gff (@original_gff) {
    if ($gff =~ /#/) {print GFFO "$gff\n";next}
    my @tmp = split /\t/,$gff;
    if (scalar @tmp eq 0) {next}
    my $name = "";
    if ($tmp[8] =~ /ID=(\S+?);/) {$name = $1}
    elsif ($tmp[8] =~ /ID=(\S+)/) {$name = $1}
    my ($original_pos,$modified_pos) = ($all_site{$tmp[0]}{$tmp[3]}{"boundary"}."..".$all_site{$tmp[0]}{$tmp[4]}{"boundary"},$tmp[3]."..".$tmp[4]);
    my $change = ($all_site{$tmp[0]}{$tmp[4]}{"boundary"}-$all_site{$tmp[0]}{$tmp[3]}{"boundary"})-($tmp[4]-$tmp[3]);
    if ($change ne 0 and $tmp[2] =~ /CDS/i) {
        $count_gene++;
        if (abs($change) % 3 ne 0) {push @changed_feature,$tmp[0]."\t$name\t$original_pos\t$modified_pos\t$name\t$change\tFrameshift"}
        else {push @changed_feature,$tmp[0]."\t$name\t$original_pos\t$modified_pos\t$name\t$change"}
    }
    $tmp[3] = $all_site{$tmp[0]}{$tmp[3]}{"boundary"};
    $tmp[4] = $all_site{$tmp[0]}{$tmp[4]}{"boundary"};
    print GFFO join("\t",@tmp),"\n";
}

foreach my $i (sort keys %fasta) {print FASTAO ">$i\n",$fasta{$i},"\n"}

print ERROR "Chr\tID\tOriginal_pos\tModified_pos\tLength change (bp)\tNote\n";
print ERROR join("\n",@changed_feature),"\n";

$count =~ s/(?<=^\d)(?=(\d\d\d)+$)|(?<=^\d\d)(?=(\d\d\d)+$)/,/gx;
$count_gene =~ s/(?<=^\d)(?=(\d\d\d)+$)|(?<=^\d\d)(?=(\d\d\d)+$)/,/gx;
say "A total of $count modifications applied, influencing $count_gene CDS regions.";
say "Done!";

sub printHelp {
print << "END";
    Usage: perl refresh_seq_gff.pl -f <sequence> -g <gff file> -m <modification file> -fm <modified sequence> -gm <modified gff> -o <summary>

    <sequence>          Original sequence file.
    <gff file>          Original gff3 file.
    <modification file> SNP or INDEL required to be modified.
                        Format:  sequence_name Position Original_seq Modified_seq # tab-separated; Position is 1-based
                        Example: Scaffold_27    9106783    C    CA  # INDEL
                                 Scaffold_1     12690342   A    G   # SNP
    <modified sequence> Output modified sequence file.
    <modified gff>      Output modified gff file.
    <summary>           Output influenced CDS regions. # require CDS features in <gff file>

    -h                  Print this page.
END
}
