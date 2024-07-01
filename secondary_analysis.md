#### Secondary analysis

Quality Control
- detonate to assess quality of assembly - compared
- filter out contigs with <10 RSEM estimate
- ARSyNseq for batch correction
- TMM normalization

Differential Expression
- NOISeq

DE plots
- volcano plots

Pathway Analysis
- network
- bar plots


##### Detonate to assess quality of assemblies

module load gcc/5.2.0
module load detonate/1.11

rsem-eval-calculate-score --no-qualities --bowtie2 --num-threads 20 ~/metatranscriptome/assembly/trinity/trinity_non_rrna_p_up_75_out/both.fa ~/metatranscriptome/assembly/trinity/trinity_non_rrna_p_up_75_out/Trinity.fasta trinity_non_rrna_p_up_75_out 150
----
rsem-eval-calculate-score --no-qualities --bowtie2 --num-threads 20 ~/metatranscriptome/assembly/trinity/trinity_non_rrna_p_up_bigmem_out/both.fa ~/metatranscriptome/assembly/trinity/trinity_non_rrna_p_up_bigmem_out/Trinity.fasta trinity_non_rrna_p_up_bigmem_fasta_out 150
----
rsem-eval-calculate-score --no-qualities --bowtie2 --num-threads 10 ~/metatranscriptome/assembly/trinity/trinity_out/both.fa ~/metatranscriptome/assembly/trinity/trinity_out/Trinity.fasta trinity_non_rrna_p_up_fasta_out 150
----
rsem-eval-calculate-score --no-qualities --paired-end ~/metatranscriptome/seqs/trimmed/merged/sortmerna/non_rRNA/mt_non_rrna_f_up.fasta ~/metatranscriptome/seqs/trimmed/merged/sortmerna/non_rRNA/mt_non_rrna_r.fasta --bowtie2 --num-threads 10 ~/metatranscriptome/assembly/idba-mt/idbaud_out_copy/contig.fa idbaud_out_copy 150
----
rsem-eval-calculate-score --no-qualities ~/metatranscriptome/seqs/trimmed/merged/sortmerna/non_rRNA/mt_non_rrna_merged.fasta --bowtie2 --num-threads 20 ~/metatranscriptome/assembly/idba-mt/idbaud_out_memmod/contig.fa idbaud_out_memmod 150
----
rsem-eval-calculate-score --no-qualities --bowtie2 --num-threads 20 ~/metatranscriptome/assembly/trinity/trinity_non_rrna_p_up_75_out/both.fa ~/metatranscriptome/assembly/trinity/trinity_non_rrna_p_up_75_out/Trinity.fasta trinity_non_rrna_p_up_75_out 150
----

#####



##### END #####
