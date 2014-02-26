#!/bin/perl -w

use Bio::Seq;
use Bio::SeqIO;
use lib qw(/mnt/PDanalysis/scripts/vcftools_0.1.10/perl);
use Vcf;

if ($#ARGV != 2) {
    print "usage: perl annotate_snps.pl myReference.fasta mySNPs.vcf myOrfs.txt\n";
    exit;
}

my $seqio_obj = Bio::SeqIO->new(-file=> "ARGV[0]", -format=> "fasta" );
my $vcf= Vcf->new(file=> "ARGV[1]");
#import the reference FASTA file and the VCF file


my %refhash = ();
#initiate the hash

while (my $seq_obj = $seqio_obj->next_seq) {
#load the data fields from each line into variables
    my $id = $seq_obj->display_id;
    my $desc = $seq_obj->desc;
    my $refseq = $seq_obj->seq;
    $refhash{ $id } = $refseq;
#add key (contig id) and value (reference sequences)
}

my $seqio_orf =Bio::SeqIO->new(-file=>"ARGV[2]", -format=> "fasta", -alphabet=> "protein" );
#import the output file from OrfPredictor
my %orfseqhash = ();
my %orfdeschash = ();
while (my $prot_obj = $seqio_orf->next_seq) {
    my $id_orf = $prot_obj->display_id;
    my $desc_orf = $prot_obj->desc;
    my $orfseq = $prot_obj->seq;
    $orfseqhash{ $id_orf } = $orfseq;
    $orfdeschash{ $id_orf } = $desc_orf;

}
#load a hash of the OrfPredictor output file, and store its contents (including predicted translated ORF peptide sequence) as variables to be compared later



$vcf -> parse_header();
#remove and store the VCF header lines


print $vcf->format_header();


while (my $x=$vcf->next_data_array()) {
#read through the VCF SNP-by-SNP
    my $contig=$$x[0];
    my $pos=$$x[1];
    my $ref=$$x[3];
    my $alt=$$x[4];
    my @splitref = split( //, $refhash{ $contig } );
    $splitref[$pos-1]=$alt;
#substitute the alternate allele into the reference sequence, -1 since perl is 0-based
    my $altcontig=join('', @splitref);
#transform the array into a new string
    if (exists $orfdeschash{ $contig }) {
	my @splitorfdesc = split( /\s/, $orfdeschash{ $contig } );
	my $frame = $splitorfdesc[0];
	my $start = $splitorfdesc[1];
	my $stop = $splitorfdesc[2];
	my $fs = $splitorfdesc[3];
    #attempt to slice and translate reference and make it look like the OrfPredictor entry
        my $altpeptide = '';
	#reset the altpeptide variable each loop, just in case
	if ($frame > 0 ) {
#if the ORF is on the coding strand, translate it with the SNP substituted in
	    my @translate = split( //, $altcontig );
	    my @altorf = @translate[$start-1..$stop-1];
	    my $altorfstring = join('', @altorf );
	    my $altorfseq = Bio::Seq->new(-seq => $altorfstring, -alphabet => 'dna');
	    $altpeptide = $altorfseq->translate;
	    $altpeptide = $altpeptide->seq;
	}
	else {
#if the ORF is on the noncoding strand, reverse complement it and translate it with the SNP substituted in
	    my $altcontig_obj = Bio::Seq->new(-seq => $altcontig, -alphabet => 'dna');
	    my $altcontig_revcom = $altcontig_obj->revcom;
	    $altcontig_revcom = $altcontig_revcom->seq;
	    my @translate = split( //, $altcontig_revcom );
	    my @altorf = @translate[$start-1..$stop-1];
	    my $altorfstring = join('', @altorf );
            my $altorfseq = Bio::Seq->new(-seq => $altorfstring, -alphabet => 'dna');
            $altpeptide = $altorfseq->translate;
	    $altpeptide= $altpeptide->seq;
	}
	
	if (defined $fs) {
	    if ( $fs eq "FS" ) {
		$$x[7]=$vcf->add_info_field($$x[7],'EF'=>'fs'); print join ("\t",@$x)."\n";
	    }
	}
	
	elsif ( ($pos < $start) || ($pos > $stop) ) {
	    $$x[7]=$vcf->add_info_field($$x[7],'EF'=>'utr'); print join ("\t",@$x)."\n";
	}

	elsif ( $altpeptide eq $orfseqhash{ $contig } ) {
#then compare the peptide sequence that is the product of the alternate allele to the reference
	    $$x[7]=$vcf->add_info_field($$x[7],'EF'=>'syn'); print join ("\t",@$x)."\n";
	}
	else {
	    $$x[7]=$vcf->add_info_field($$x[7],'EF'=>'ns'); print join ("\t",@$x)."\n";
	}
    }
}


