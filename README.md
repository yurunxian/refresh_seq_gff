# refresh_seq_gff
This script is designed for integrating SNP and INDEL modifications into existing reference fasta and gff3 files and generating a new version. Any modifications related to CDS regions will be separately output into a summary file.

# Usage
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

# Example
    perl refresh_seq_gff.pl -f example/sequence.fasta -g example/gene.gff -m example/modification.txt -o example/summary.txt -fm example/modifed_sequence.fasta -gm example/modified_gene.gff > example/log.txt
