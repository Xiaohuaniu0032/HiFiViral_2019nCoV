# minimap2 align
/root/tools/git_repo/HiFiViral_2019nCoV/bin/minimap2-2.24_x64-linux/minimap2 -a /root/tools/git_repo/HiFiViral_2019nCoV/ref/NC_045512.2.fasta /root/tools/git_repo/HiFiViral_2019nCoV/data/C_D3R7.hifi_reads.fastq >/root/tools/git_repo/HiFiViral_2019nCoV/test/C_D3R7/C_D3R7.sam

/root/tools/biosoft/samtools-1.16.1/samtools view -b -S /root/tools/git_repo/HiFiViral_2019nCoV/test/C_D3R7/C_D3R7.sam | /root/tools/biosoft/samtools-1.16.1/samtools sort >/root/tools/git_repo/HiFiViral_2019nCoV/test/C_D3R7/C_D3R7.sort.bam

/root/tools/biosoft/samtools-1.16.1/samtools index /root/tools/git_repo/HiFiViral_2019nCoV/test/C_D3R7/C_D3R7.sort.bam

/root/tools/biosoft/samtools-1.16.1/samtools mpileup --min-BQ 1 -f /root/tools/git_repo/HiFiViral_2019nCoV/ref/NC_045512.2.fasta /root/tools/git_repo/HiFiViral_2019nCoV/test/C_D3R7/C_D3R7.sort.bam >/root/tools/git_repo/HiFiViral_2019nCoV/test/C_D3R7/C_D3R7.mpileup

/root/tools/biosoft/samtools-1.16.1/samtools depth -q 0 -Q 0 /root/tools/git_repo/HiFiViral_2019nCoV/test/C_D3R7/C_D3R7.sort.bam >/root/tools/git_repo/HiFiViral_2019nCoV/test/C_D3R7/C_D3R7.bam.depth

# call variants using bcftools
/root/miniconda3/bin/bcftools mpileup -f /root/tools/git_repo/HiFiViral_2019nCoV/ref/NC_045512.2.fasta /root/tools/git_repo/HiFiViral_2019nCoV/test/C_D3R7/C_D3R7.sort.bam | /root/miniconda3/bin/bcftools call -mv -Ov -o /root/tools/git_repo/HiFiViral_2019nCoV/test/C_D3R7/bcftools/C_D3R7.bcftools.vcf

# generate consensus using VCFCons.py
/root/miniconda3/bin/python3 /root/tools/git_repo/HiFiViral_2019nCoV/bin/VCFCons.py /root/tools/git_repo/HiFiViral_2019nCoV/ref/NC_045512.2.fasta /root/tools/git_repo/HiFiViral_2019nCoV/test/C_D3R7/bcftools -c 4 -f 0.5 --input_depth /root/tools/git_repo/HiFiViral_2019nCoV/test/C_D3R7/C_D3R7.bam.depth --input_vcf /root/tools/git_repo/HiFiViral_2019nCoV/test/C_D3R7/bcftools/C_D3R7.bcftools.vcf --vcf_type bcftools

/root/tools/biosoft/samtools-1.16.1/samtools fqidx /root/tools/git_repo/HiFiViral_2019nCoV/data/C_D3R7.hifi_reads.fastq

/root/miniconda3/bin/pbaa cluster --min-cluster-read-count 2 --trim-ends 0 /root/tools/git_repo/HiFiViral_2019nCoV/ref/sarscov2_guide_for_pbaa.fasta /root/tools/git_repo/HiFiViral_2019nCoV/data/C_D3R7.hifi_reads.fastq /root/tools/git_repo/HiFiViral_2019nCoV/test/C_D3R7/pbaa

