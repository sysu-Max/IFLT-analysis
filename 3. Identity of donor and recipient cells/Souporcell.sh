##-----note: Those three files including reference genome, bam, and known vcf files, must in same folder. 
~/Miniconda3/envs/Singularity/bin/singularity exec ~/RAW_data/3_Bei_RAW_Fastq/LT/souporcell_latest.sif souporcell_pipeline.py \
-i possorted_genome_bam.bam \
-b sample_cDNA.RData_cleanBarcode.txt \
-f refdata-cellranger-GRCh38_18Y_NC_007605.1.fa \
--skip_remap SKIP_REMAP \
--common_variants filtered_2p_1kgenomes_unchr_for_Soupcells.vcf  \
-t 15 \
-o sample_output \
-k 2



####
