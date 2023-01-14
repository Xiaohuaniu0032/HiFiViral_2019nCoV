# HiFiViral_2019nCoV



# Test data
`https://www.pacb.com/connect/datasets/`
Click Dataset `Omicron samples` and you will see `https://downloads.pacbcloud.com/public/dataset/HiFiViral/Jan_2022/`

## Fastq
Download `samples.hifi_reads.fastq.zip`: Processed HiFi reads: demultiplexed, trimmed, and QCed (96 samples)

## Aligned BAM
Download `samples.mapped.bam.zip`: Processed HiFi reads mapped to Wuhan reference NC.045512.2 (96 samples)

> We use 


# PacBio Documents

https://github.com/PacificBiosciences/CoSA/wiki/Variant-calling-using-PacBio-HiFi-CCS-data

This tutorial includes:
1. generate ccs data [Optional]
2. (if multiplexed) de-multiplex different patient samples [Optional]
3. amplicon primer trimming
4. read mapping using minimap2 or pbmm2
5. variant calling using bcftools or DeepVariant or pbaa
6. generating consensus sequence using VCFCons

