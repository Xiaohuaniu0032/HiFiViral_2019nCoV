use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use FindBin qw/$Bin/;

my ($fastq,$name,$samtools,$bcftools,$python3,$outdir);

GetOptions(
	"fq:s"        => \$fastq,      # Need
	"n:s"         => \$name,       # Need
	"samtools:s"  => \$samtools,   # Need
	"bcftools:s"  => \$bcftools,   # Need
	"py3:s"       => \$python3,    # Need
	"od:s"        => \$outdir,     # Need
	) or die "unknown args\n";


if (not defined $fastq || not defined $name || not defined $outdir){
	die "please specify: -fq <fq> -n <name> -od <outdir>\n";
}

if (not defined $samtools || not defined $bcftools || not defined $python3){
	die "please specify: -mm2 <minimap2_bin> -samtools <samtools_bin> -bcftools <bcftools> -py3 <python3>\n";
}

# main steps
# 1. Variant Calling using bcftools
# 2. Variant Calling using pbaa
# 3. Generating consensus sequence using VCFCons


if (!-d "$outdir/$name"){
	`mkdir $outdir/$name`;
}

if (!-d "$outdir/$name/bcftools"){
	`mkdir $outdir/$name/bcftools`;
}

if (!-d "$outdir/$name/pbaa"){
	`mkdir $outdir/$name/pbaa`; # PB AMPLICON ANALYSIS
}


my $runsh = "$outdir/$name/run\_$name\.sh";

open SH, ">$runsh" or die;

# 1. variant calling using bcftools
# bamtools convert -format fastq -in per_patient.primer_trimmed.bam > per_patient.primer_trimmed.fastq
print SH "\# minimap2 align\n";
my $ref = "$Bin/CoSA/data/NC_045512.2.fasta";
my $minimap2 = "$Bin/bin/minimap2-2.24_x64-linux/minimap2";

my $sam = "$outdir/$name\.sam";
# fastq -> SAM
my $cmd = "$minimap2 -a $ref $fastq >$sam";
print SH "$cmd\n\n";

# SAM -> BAM
my $sort_bam = "$outdir/$name\.sort.bam";
$cmd = "$samtools view -b -S $sam \| $samtools sort >$sort_bam";
# -b: output BAM
# -S: input is SAM
print SH "$cmd\n\n";

# samtools index
$cmd = "$samtools index $sort_bam";
print SH "$cmd\n\n";

# cal pileup
my $pileup = "$outdir/$name\.mpileup";
$cmd = "$samtools mpileup --min-BQ 1 -f $ref $sort_bam >$pileup";
print SH "$cmd\n\n";

# cal depth
my $depth = "$outdir/$name\.bam.depth"; # .bam.depth will be used by VCFCons.py
$cmd = "$samtools depth -q 0 -Q 0 $sort_bam >$depth";
print SH "$cmd\n\n";

# call variants using bcftools
print SH "\# call variants using bcftools\n";
my $bcftools_vcf = "$outdir/bcftools/$name\.bcftools.vcf";
$cmd = "$bcftools mpileup -f $ref $sort_bam | $bcftools call -mv -Ov -o $bcftools_vcf";
print SH "$cmd\n\n";


# filter variants & generate consensus using VCFCons
# https://github.com/Magdoll/CoSA/blob/master/vcf/run_VCFCons_bcftools.sh
print SH "\# generate consensus using VCFCons.py\n";
my $VCFCons_py = "$Bin/CoSA/vcf/VCFCons.py";
my $prefix = "$outdir/$name";
my $MIN_COVERAGE = 4;
my $MIN_ALT_FREQ = 0.5;
$cmd = "$python3 $VCFCons_py $ref $prefix -c $MIN_COVERAGE -f $MIN_ALT_FREQ --input_depth $depth --input_vcf $bcftools_vcf --vcf_type bcftools";
print SH "$cmd";

close SH;

`chmod 755 $runsh`;







# 2. variant calling using pbaa
#my $guide_fa = "$Bin/ref/sarscov2_guide_for_pbaa.fasta"; # need to have a samtools faidx file

# Note that unlike bcftools and DeepVariant, pbaa does not take an input aligned BAM file.
# Instead we convert the primer-trimmed (post-lima) BMA file into FASTQ format and run pbaa on it.


# bamtools convert -format fastq -in per_patient.primer_trimed.bam > sample.fastq
# samtools faidx sample.fastq

#my $cmd = "$pbaa cluster --min-cluster-read-count 2 --trim-ends 0 $guide_fa "
