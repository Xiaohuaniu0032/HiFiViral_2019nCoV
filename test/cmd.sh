fq="C_D3R7.hifi_reads.fastq"
n="C_D3R7"

samtools_bin="/path/to/samtools"
bcftools_bin="/path/to/bcftools"
py3_bin="/path/to/python3"

perl ../HiFiViral_2019nCoV_main.pl -fq $fq -n $n -samtools $samtools_bin -bcftools $bcftools_bin -py3 $py3_bin -od $PWD
