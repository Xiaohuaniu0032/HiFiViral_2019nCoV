# minimap2 align
/root/tools/git_repo/HiFiViral_2019nCoV/bin/minimap2-2.24_x64-linux/minimap2 -a /root/tools/git_repo/HiFiViral_2019nCoV/CoSA/data/NC_045512.2.fasta C_D3R7.hifi_reads.fastq >/root/tools/git_repo/HiFiViral_2019nCoV/test/C_D3R7.sam

/path/to/samtools view -b -S /root/tools/git_repo/HiFiViral_2019nCoV/test/C_D3R7.sam | /path/to/samtools sort >/root/tools/git_repo/HiFiViral_2019nCoV/test/C_D3R7.sort.bam

/path/to/samtools index /root/tools/git_repo/HiFiViral_2019nCoV/test/C_D3R7.sort.bam

/path/to/samtools mpileup --min-BQ 1 -f /root/tools/git_repo/HiFiViral_2019nCoV/CoSA/data/NC_045512.2.fasta /root/tools/git_repo/HiFiViral_2019nCoV/test/C_D3R7.sort.bam >/root/tools/git_repo/HiFiViral_2019nCoV/test/C_D3R7.mpileup

/path/to/samtools depth -q 0 -Q 0 /root/tools/git_repo/HiFiViral_2019nCoV/test/C_D3R7.sort.bam >/root/tools/git_repo/HiFiViral_2019nCoV/test/C_D3R7.bam.depth

# call variants using bcftools
/path/to/bcftools mpileup -f /root/tools/git_repo/HiFiViral_2019nCoV/CoSA/data/NC_045512.2.fasta /root/tools/git_repo/HiFiViral_2019nCoV/test/C_D3R7.sort.bam | /path/to/bcftools call -mv -Ov -o /root/tools/git_repo/HiFiViral_2019nCoV/test/bcftools/C_D3R7.bcftools.vcf

# generate consensus using VCFCons.py
/path/to/python3 /root/tools/git_repo/HiFiViral_2019nCoV/CoSA/vcf/VCFCons.py /root/tools/git_repo/HiFiViral_2019nCoV/CoSA/data/NC_045512.2.fasta /root/tools/git_repo/HiFiViral_2019nCoV/test/C_D3R7 -c 4 -f 0.5 --input_depth /root/tools/git_repo/HiFiViral_2019nCoV/test/C_D3R7.bam.depth --input_vcf /root/tools/git_repo/HiFiViral_2019nCoV/test/bcftools/C_D3R7.bcftools.vcf --vcf_type bcftools