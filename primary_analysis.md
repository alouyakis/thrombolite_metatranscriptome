#### [A year in the life of a thrombolite: comparative metatranscriptomics reveals dynamic metabolic changes over diel and seasonal cycles](https://enviromicro-journals-onlinelibrary-wiley-com.ezproxy.lib.uconn.edu/doi/full/10.1111/1462-2920.14029)  
Artemis S. Louyakis, Hadrien Gourlé, Giorgio Casaburi, Rachelle M. E. Bonjawo, Alexandrea A. Duscher, Jamie S. Foster  

Samples were collected at 6am, 12pm, 12am, and 6pm over three days in three seasons/months March, August, and October and sequenced (2x150, 550 bp insert) in duplicate with each duplicate being a separate day. Additional samples are stored in RNAlater at -80&deg in the Foster Lab. This study examined 24 samples.  

##### Primary analysis

Quality check
- FastQC v0.11.4

```bash
parallel -j 10 'fastqc rawseqs/{1}_r1.fastq rawseqs/{1}_r2.fastq -o quality -t 5 &&
  echo "FastQC Complete for {1}" && date' ::: $( basename -a rawseqs/*_r1.fastq | cut -f 1 -d '_')
```

Trim reads
- Trimmomatic v0.36
- ran through Trinity pipeline - ran first on raw reads, then again on rRNA-free reads; no need to do that

rRNA removal
- Sortmerna v2.1

```bash
list=$(for i in samps.txt; do awk '{print $1}' $i; done)
for SAMPLE in ${list}; do
  ./merge-paired-reads.sh ${SAMPLE}.qtrim.r1.fq ${SAMPLE}.qtrim.r2.fq ${SAMPLE}.qtrim.merge.fq;
done

parallel -j 5 'rrna_db="/databases/sortmerna/data/rRNA_databases" &&
  sortmerna --ref ${rrna_db}/rfam-5.8s-database-id98.fasta \
    --ref ${rrna_db}/rfam-5s-database-id98.fasta \
    --ref ${rrna_db}/silva-arc-16s-id95.fasta \
    --ref ${rrna_db}/silva-arc-23s-id98.fasta \
    --ref ${rrna_db}/silva-bac-16s-id90.fasta \
    --ref ${rrna_db}/silva-bac-23s-id98.fasta \
    --ref ${rrna_db}/silva-euk-18s-id95.fasta \
    --ref ${rrna_db}/silva-euk-28s-id98.fasta \
    --reads {1}.qtrim.r1.fq \
    --reads {1}.qtrim.r2.fq \
    --aligned sortmerna/rrna_{1} \
    --other sortmerna/mrna_{1} --out2 \
    --paired_in --fastx --num_alignments 1 \
    --threads 20 -m 20000 --workdir sortmerna/{1} &&
  echo "rRNA removal complete for {1}" && date' ::: $( basename -a *.qtrim.r1.fq | cut -f 1 -d '_')
wait

mkdir -p sortmerna/unmerge
list=$(for i in samps.txt; do awk '{print $1}' $i; done)
for i in $list; do \
  ./unmerge-paired-reads.sh sortmerna/${i}_non_rrna.fq sortmerna/unmerge/${i}_non_rrna_r1.fq sortmerna/unmerge/${i}_non_rrna_r2.fq &
done
```

Assembly
- Trinity v2.2.0

```bash
## run assembly
Trinity --seqType fq --max_memory 400G --trimmomatic --normalize_reads --CPU 20 --output trinity_out \
  --left \
    ad1t12_non_rrna_r1.fastq,ad1t18_non_rrna_r1.fastq,ad1t24_non_rrna_r1.fastq,ad1t6_non_rrna_r1.fastq,ad3t12_non_rrna_r1.fastq,ad3t18_non_rrna_r1.fastq,ad3t24_non_rrna_r1.fastq,ad3t6_non_rrna_r1.fastq,md2t12_non_rrna_r1.fastq,md2t18_non_rrna_r1.fastq,md2t24_non_rrna_r1.fastq,md2t6_non_rrna_r1.fastq,md3t12_non_rrna_r1.fastq,md3t18_non_rrna_r1.fastq,md3t24_non_rrna_r1.fastq,md3t6_non_rrna_r1.fastq,od2t12_non_rrna_r1.fastq,od2t18_non_rrna_r1.fastq,od2t24_non_rrna_r1.fastq,od2t6_non_rrna_r1.fastq,od3t12_non_rrna_r1.fastq,od3t18_non_rrna_r1.fastq,od3t24_non_rrna_r1.fastq,od3t6_non_rrna_r1.fastq \
  --right \
    ad1t12_non_rrna_r2.fastq,ad1t18_non_rrna_r2.fastq,ad1t24_non_rrna_r2.fastq,ad1t6_non_rrna_r2.fastq,ad3t12_non_rrna_r2.fastq,ad3t18_non_rrna_r2.fastq,ad3t24_non_rrna_r2.fastq,ad3t6_non_rrna_r2.fastq,md2t12_non_rrna_r2.fastq,md2t18_non_rrna_r2.fastq,md2t24_non_rrna_r2.fastq,md2t6_non_rrna_r2.fastq,md3t12_non_rrna_r2.fastq,md3t18_non_rrna_r2.fastq,md3t24_non_rrna_r2.fastq,md3t6_non_rrna_r2.fastq,od2t12_non_rrna_r2.fastq,od2t18_non_rrna_r2.fastq,od2t24_non_rrna_r2.fastq,od2t6_non_rrna_r2.fastq,od3t12_non_rrna_r2.fastq,od3t18_non_rrna_r2.fastq,od3t24_non_rrna_r2.fastq,od3t6_non_rrna_r2.fastq

/apps/trinity/r20160329-2.2.0/util/support_scripts/get_Trinity_gene_to_trans_map.pl Trinity.fasta > Trinity.fasta.gene_trans_map

TrinityStats.pl Trinity.fasta

## gff for downstream plot
gmhmmp -d -f g -m /apps/genemark/metagenemark/3.26/bin/MetaGeneMark_v1.mod Trinity.filtered.fasta -o Trinity.filtered.gff
```

Alignment
- Bowtie2 (via Trinity pipeline)
Abundance estimation
- RSEM (via Trinity pipeline)

```bash
## run as array with slurm
RUN=${SLURM_ARRAY_TASK_ID}
DATA_DIR="~/metatranscriptome/rawseqs"
ASSEMBLY="~/metatranscriptome/trinity_out"

INPUT_PATH=$(ls ${DATA_DIR}/*non_rrna_r1.fq | sed -n "${RUN} p" | sed 's/_non_rrna_r1.fq//' )
SAMPLE=$(basename ${INPUT_PATH})

/isg/shared/apps/trinity/2.6.6/util/align_and_estimate_abundance.pl \
  --transcripts ${ASSEMBLY}/Trinity.fasta \
  --seqType fq \
  --left ${DATA_DIR}/${SAMPLE}_non_rrna_r1.fq \
  --right ${DATA_DIR}/${SAMPLE}_non_rrna_r2.fq \
  --est_method RSEM \
  --output_dir ${SAMPLE} \
  --aln_method bowtie2 \
  --thread_count 10 \
  --trinity_mode \
  --output_prefix ${SAMPLE} \
  --include_rsem_bam

## create matrix and estimate counts
abundance_estimates_to_matrix.pl \
  --est_method RSEM \
  --out_prefix iso \
  --name_sample_by_basedir \
  ad1t6/ad1t6.isoforms.results \
  ad1t12/ad1t12.isoforms.results \
  ad1t18/ad1t18.isoforms.results \
  ad1t24/ad1t24.isoforms.results \
  ad3t6/ad3t6.isoforms.results \
  ad3t12/ad3t12.isoforms.results \
  ad3t18/ad3t18.isoforms.results \
  ad3t24/ad3t24.isoforms.results \
  md2t6/md2t6.isoforms.results \
  md2t12/md2t12.isoforms.results \
  md2t18/md2t18.isoforms.results \
  md2t24/md2t24.isoforms.results \
  md3t6/md3t6.isoforms.results \
  md3t12/md3t12.isoforms.results \
  md3t18/md3t18.isoforms.results \
  md3t24/md3t24.isoforms.results \
  od2t6/od2t6.isoforms.results \
  od2t12/od2t12.isoforms.results \
  od2t18/od2t18.isoforms.results \
  od2t24/od2t24.isoforms.results \
  od3t6/od3t6.isoforms.results \
  od3t12/od3t12.isoforms.results \
  od3t18/od3t18.isoforms.results \
  od3t24/od3t24.isoforms.results
```

Annotations
- Trinotate v3.0.1 pipeline
- gcc v5.2.0
- gnuplot v5.0.3

```bash
## transdecoder
TransDecoder.LongOrfs -t Trinity.fasta
TransDecoder.Predict -t Trinity.fasta

## annotate contigs against uniprot swissprot db
blastx -query Trinity.fasta -db uniprot_sprot.pep -num_threads 15 -max_target_seqs 1 -outfmt 6 > blastx.outfmt6
blastp -query Trinity.fasta.transdecoder.pep -db uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 > blastp.outfmt6

## hmmscan
hmmscan --cpu 20 --domtblout trinotate.pfam.out Pfam-A.hmm Trinity.fasta.transdecoder.pep

## tmhmm
tmhmm --short Trinity.fasta.transdecoder.pep tmhmm.out

## rnammer on contigs
RnammerTranscriptome.pl --org_type bac --transcriptome Trinity.fasta --path_to_rnammer /apps/trinotate/rnammer/1.2/rnammer
RnammerTranscriptome.pl --org_type arc --transcriptome Trinity.fasta --path_to_rnammer /apps/trinotate/rnammer/1.2/rnammer
RnammerTranscriptome.pl --org_type euk --transcriptome Trinity.fasta --path_to_rnammer /apps/trinotate/rnammer/1.2/rnammer

## signalp for various taxa groups on transdecoder and pfam protein files
signalp -f short -g png+eps -t euk -n signalp_euk.out Trinity.fasta.transdecoder.pep
signalp -f short -g png+eps -t gram- -n signalp_gramn.out Trinity.fasta.transdecoder.pep
signalp -f short -g png+eps -t gram+ -n signalp_gramp.out Trinity.fasta.transdecoder.pep

## get trinotate sql database template
wget "https://data.broadinstitute.org/Trinity/Trinotate_v3_RESOURCES/Trinotate_v3.sqlite.gz" -O Trinotate.sqlite.gz

## initialize db
Trinotate Trinotate.sqlite init --gene_trans_map Trinity.filtered.fasta.gene_trans_map \
  --transcript_fasta Trinity.fasta --transdecoder_pep Trinity.fasta.transdecoder.pep

## load db with annotations
Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.outfmt6
Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6
Trinotate Trinotate.sqlite LOAD_pfam trinotate.pfam.out
Trinotate Trinotate.sqlite LOAD_tmhmm tmhmm.out
Trinotate Trinotate.sqlite LOAD_signalp signalp.euk.out
Trinotate Trinotate.sqlite LOAD_signalp signalp.gramn.out
Trinotate Trinotate.sqlite LOAD_signalp signalp.gramp.out
Trinotate Trinotate.sqlite LOAD_rnammer rnammer.bac/Trinity.filtered.fasta.rnammer.gff
Trinotate Trinotate.sqlite LOAD_rnammer rnammer.arc/Trinity.filtered.fasta.rnammer.gff
Trinotate Trinotate.sqlite LOAD_rnammer rnammer.euk/Trinity.filtered.fasta.rnammer.gff

## output reports filtering for evalue
Trinotate Trinotate.sqlite report > trinotate_annotation_report.xls
Trinotate Trinotate.sqlite report -E 0.00001 > trinotate_annotation_report_1e5.xls
Trinotate Trinotate.sqlite report -E 0.0001 > trinotate_annotation_report_1e4.xls


## extract GO assignments (if needed)
extract_GO_assignments_from_Trinotate_xls.pl  \
--Trinotate_xls trinotate_annotation_report.xls \
-G --include_ancestral_terms \
> go_annotations.txt


## get name map
Trinotate_get_feature_name_encoding_attributes.pl trinotate_annotation_report.xls > trinotate_annotation_report.xls.name_mappings
```

#### Secondary analysis

Quality Control
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
