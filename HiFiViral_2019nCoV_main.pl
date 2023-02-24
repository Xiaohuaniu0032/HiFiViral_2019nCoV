use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use FindBin qw/$Bin/;

my ($fastq,$name,$samtools,$bcftools,$python3,$pbaa,$outdir);

GetOptions(
	"fq:s"        => \$fastq,      # Need
	"n:s"         => \$name,       # Need
	"samtools:s"  => \$samtools,   # Default: /root/miniconda3/bin/samtools
	"bcftools:s"  => \$bcftools,   # Default: /root/miniconda3/bin/bcftools
	"py3:s"       => \$python3,    # Default: /usr/bin/python3
	"pbaa:s"      => \$pbaa,       # Default: /root/miniconda3/bin/pbaa
	"od:s"        => \$outdir,     # Need
	) or die "unknown args\n";


if (not defined $fastq || not defined $name || not defined $outdir){
	die "please specify: -fq <fq> -n <name> -od <outdir>\n";
}


# default value

if (not defined $samtools){
	$samtools = "/root/tools/biosoft/samtools-1.16.1/samtools";
}

if (not defined $bcftools){
	$bcftools = "/root/miniconda3/bin/bcftools";
}

if (not defined $python3){
	$python3 = "/root/miniconda3/bin/python3";
}

if (not defined $pbaa){
	$pbaa = "/root/miniconda3/bin/pbaa";
}




########################## Main Steps ##########################
# 1.Generate CCS data
# 2.(if multiplexed) de-multiplex different patient samples
# 3.Trim amplicon primers
# 4.Variant Calling using bcftools
# 5.Generating consensus sequence using VCFCons.py



if (!-d "$outdir/$name"){
	`mkdir $outdir/$name`;
}

if (!-d "$outdir/$name/bcftools"){
	`mkdir $outdir/$name/bcftools`;
}

if (!-d "$outdir/$name/pbaa"){
	`mkdir $outdir/$name/pbaa`; # PB AMPLICON ANALYSIS
}


# abs path
my $ref = "$Bin/ref/NC_045512.2.fasta";
my $minimap2 = "$Bin/bin/minimap2-2.24_x64-linux/minimap2";




###################### run shell ######################
my $runsh = "$outdir/$name/run\_$name\.sh";

open SH, ">$runsh" or die;

###################################### Step.1 variant calling using bcftools
print SH "\# minimap2 align\n";

my $sam = "$outdir/$name/$name\.sam";
# fastq -> SAM
my $cmd = "$minimap2 -a $ref $fastq >$sam";
print SH "$cmd\n\n";

# SAM -> BAM
my $sort_bam = "$outdir/$name/$name\.sort.bam";
$cmd = "$samtools view -b -S $sam \| $samtools sort >$sort_bam";
# -b: output BAM
# -S: input is SAM
print SH "$cmd\n\n";

# samtools index
$cmd = "$samtools index $sort_bam";
print SH "$cmd\n\n";

# cal pileup
my $pileup = "$outdir/$name/$name\.mpileup";
$cmd = "$samtools mpileup --min-BQ 1 -f $ref $sort_bam >$pileup";
print SH "$cmd\n\n";

# cal depth
my $depth = "$outdir/$name/$name\.bam.depth"; # .bam.depth will be used by VCFCons.py
$cmd = "$samtools depth -q 0 -Q 0 $sort_bam >$depth";
print SH "$cmd\n\n";

# call variants using bcftools
print SH "\# call variants using bcftools\n";
my $bcftools_vcf = "$outdir/$name/bcftools/$name\.bcftools.vcf";
$cmd = "$bcftools mpileup -f $ref $sort_bam | $bcftools call -mv -Ov -o $bcftools_vcf";
print SH "$cmd\n\n";


# filter variants & generate consensus using VCFCons
# https://github.com/Magdoll/CoSA/blob/master/vcf/run_VCFCons_bcftools.sh
print SH "\# generate consensus using VCFCons.py\n";

my $VCFCons_py = "$Bin/bin/VCFCons.py";
my $prefix_bcftools = "$outdir/$name/bcftools";

my $MIN_COVERAGE = 4;
my $MIN_ALT_FREQ = 0.5;
$cmd = "$python3 $VCFCons_py $ref $prefix_bcftools -c $MIN_COVERAGE -f $MIN_ALT_FREQ --input_depth $depth --input_vcf $bcftools_vcf --vcf_type bcftools";
print SH "$cmd\n\n";






#################################### Step.2 variant calling using pbaa
my $guide_fa = "$Bin/ref/sarscov2_guide_for_pbaa.fasta"; # need to have a samtools faidx file

# Note that unlike bcftools and DeepVariant, pbaa does not take an input aligned BAM file.
# Instead we convert the primer-trimmed (post-lima) BMA file into FASTQ format and run pbaa on it.


# bamtools convert -format fastq -in per_patient.primer_trimed.bam > sample.fastq
$cmd = "$samtools fqidx $fastq";
print SH "$cmd\n\n";


$cmd = "$pbaa cluster --min-cluster-read-count 2 --trim-ends 0 $guide_fa $fastq $outdir/$name/pbaa";
print SH "$cmd\n\n";

#my $pbaa_dir = "$outdir/$name/pbaa";
#my $read_info = "$outdir/$name/pbaa/sample_pbaa_read_info.txt";
#my $pass_cluster = "$outdir/$name/pbaa/passed_cluster_sequences.fasta";

#$cmd = "$python3 $Bin/bin/consensusVariants.py -r $pbaa_dir -p $pbaa_dir --read_info $read_info --hifiSupport $fastq $pass_cluster";
#print SH "$cmd\n\n";


#my $pbaa_vcf = "$pbaa_dir/$name\.pbaa.vcf";
#my $pbaa_alleles = "$pbaa_dir/$name\.pbaa_alleles.csv";
#my $pbaa_variants = "$pbaa_dir/$name\.pbaa_variants.csv";

#$cmd = "$python3 $Bin/bin/pbaa2vcf.py --passOnly -s barcode -o $pbaa_vcf $pbaa_alleles $pbaa_variants $ref";
#print SH "$cmd\n\n";

# run VCFCons.py
# https://github.com/Magdoll/CoSA/blob/master/vcf/run_VCFCons_pbaa.sh
#my $prefix_pbaa = "$outdir/$name/pbaa";

#$cmd = "$python3 $VCFCons_py $ref $prefix_pbaa -c $MIN_COVERAGE -f $MIN_ALT_FREQ --input_depth $depth --vcf_type pbaa --input_vcf $pbaa_vcf";
#print SH "$cmd\n\n";





#minimap2 -a $REF ${SAMPLE}.vcfcons.frag.fasta > ${SAMPLE}.vcfcons.frag.fasta.sam
#$SAMTOOLS view -bS ${SAMPLE}.vcfcons.frag.fasta.sam > ${SAMPLE}.vcfcons.frag.fasta.aligned.bam
#$SAMTOOLS sort ${SAMPLE}.vcfcons.frag.fasta.aligned.bam > ${SAMPLE}.vcfcons.frag.fasta.sorted.aligned.bam
#$SAMTOOLS index ${SAMPLE}.vcfcons.frag.fasta.sorted.aligned.bam


close SH;

`chmod 755 $runsh`;