#!/bin/bash
###############################################
# genome assembly for Fuscarium avenaceum v1.0#
###############################################

#ran one DNA sample on MinIon with a R10.3 Version flowcell. Power outages occured twice during sequencing run requiring the use of a read recovery.

#read recover example
/opt/ont/minknow/bin/recover_reads /var/lib/minknow/data/queued_reads/complete_reads_79098d43-4d8e-4210-ba20-0ed38e317a79/ --output-directory ./genomes/fuscarium/run1_restart/

#basecalling of recovered reads
guppy_basecaller -i ../run1/ -c dna_r10.3_450bps_hac.cfg -r -s ../run1/hac_called/ --device cuda:0
guppy_basecaller -i ../run1_restart/ -c dna_r10.3_450bps_hac.cfg -r -s ../run1_restart/hac_called/ --device cuda:1

#merge fastq files
cat ./genomes/fuscarium/run1/hac_called/pass/*.fastq > ./genomes/fuscarium/fuscarium_run1.fastq
cat ./genomes/fuscarium/run1_restart/hac_called/pass/*.fastq > ./genomes/fuscarium/fuscarium_run1_restart.fastq

cat ./genomes/fuscarium/*.fastq > ./genomes/nanopore_hac.fastq
mv ./genomes/nanopore_hac.fastq ./genomes/fuscarium/

#checking quality of nanopore reads
nanoQC -o ./ ./genomes/fuscarium/nanopore_hac.fastq


NanoPlot -t 20 --fastq ./genomes/fuscarium/nanopore_hac.fastq


cat nanopore_hac.fastq | fastq-scan -g 42000000 > fuscarium_hacq-scan_out.txt &
"""{
    "qc_stats": {
        "total_bp": 7833220732,
        "coverage": 186.505,
        "read_total": 298551,
        "read_min": 10,
        "read_mean": 26237.5,
        "read_std": 21054.1,
        "read_median": 21250,
        "read_max": 261485,
        "read_25th": 10743,
        "read_75th": 36453,
        "qual_min": 11,
        "qual_mean": 19.293,
        "qual_std": 2.91784,
        "qual_max": 88,
        "qual_median": 19,
        "qual_25th": 17,
        "qual_75th": 21
    }"""


#compress read file
gzip ./genomes/fuscarium/nanopore_hac.fastq

#remove adapter sequences
porechop -i ./genomes/fuscarium/nanopore_hac.fastq.gz  -o ./genomes/fuscarium/nanopore_hac_trimmed.fastq.gz --discard_middle --threads 30

"""
Adapter seqs found:
SQK-NSK007_Y_Top
SQK-NSK007_Y_Bottom
1D2_part_2_start
1D2_part_2_end

195,186/298,551 reads had adapters trimmed from their start (7,177,988 bp removed)
103,749/298,551 reads had adapters trimmed from their end (904,118 bp removed)

171/298,551 reads were discarded based on middle adapters
"""

#remove short seqs of 10Kb or less
gunzip -c nanopore_hac_trimmed.fastq.gz | NanoFilt -l 10000 | gzip > nanopore_hac_filtered.fastq.gz

#assembly with canu
mkdir ONT_asm
cd ONT_asm
canu \
 -p fuscarium_v0.0.1 -d /./ont_asm_fuscarium \
 genomeSize=42m \
 -nanopore ./genomes/fuscarium/nanopore_hac_filtered.fastq.gz

gt seqstat fuscarium_v0.0.1.contigs.fasta
# number of contigs:     17
# total contigs length:  42517872
# mean contig size:      2501051.29
# contig size first quartile: 420365
# median contig size:         2626751
# contig size third quartile: 4726737
# longest contig:             7141040
# shortest contig:            89058
# contigs > 500 nt:           17 (100.00 %)
# contigs > 1K nt:            17 (100.00 %)
# contigs > 10K nt:           17 (100.00 %)
# contigs > 100K nt:          16 (94.12 %)
# contigs > 1M nt:            10 (58.82 %)
# N50                4726737
# L50                4
# N80                2626751
# L80                8

busco -m genome -i ./ont_asm_fuscarium/fuscarium_v0.0.1.contigs.fasta -o canu_busco -l hypocreales_odb10

chris@matrix-reloaded:~/programs/merqury-1.3$ ln -s ./ont_asm_fuscarium/fuscarium_v0.0.1.contigs.fasta canu_asm.fasta
./merqury.sh fuscarium.meryl canu_asm.fasta fuscarium_out_canu

#assembly with NECAT
necat.pl config fuscarium_config.txt

nano fuscarium_config.txt
"""
PROJECT=fuscarium v.0.1
ONT_READ_LIST=read_list.txt
GENOME_SIZE=46000000
THREADS=24
MIN_READ_LENGTH=10000
PREP_OUTPUT_COVERAGE=40
OVLP_FAST_OPTIONS=-n 500 -z 20 -b 2000 -e 0.5 -j 0 -u 1 -a 1000
OVLP_SENSITIVE_OPTIONS=-n 500 -z 10 -e 0.5 -j 0 -u 1 -a 1000
CNS_FAST_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0
CNS_SENSITIVE_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0
TRIM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 1 -a 400
ASM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 0 -a 400
NUM_ITER=2
CNS_OUTPUT_COVERAGE=30
CLEANUP=1
USE_GRID=false
GRID_NODE=0
GRID_OPTIONS=
SMALL_MEMORY=0
FSA_OL_FILTER_OPTIONS=
FSA_ASSEMBLE_OPTIONS=
FSA_CTG_BRIDGE_OPTIONS=
POLISH_CONTIGS=true
"""
necat.pl correct fuscarium_config.txt

necat.pl assemble fuscarium_config.txt

necat.pl bridge fuscarium_config.txt

gt seqstat
# number of contigs:     9
# total contigs length:  41834346
# mean contig size:      4648260.67
# contig size first quartile: 4694329
# median contig size:         4958227
# contig size third quartile: 6489478
# longest contig:             7127630
# shortest contig:            711171
# contigs > 500 nt:           9 (100.00 %)
# contigs > 1K nt:            9 (100.00 %)
# contigs > 10K nt:           9 (100.00 %)
# contigs > 100K nt:          9 (100.00 %)
# contigs > 1M nt:            8 (88.89 %)
# N50                4958227
# L50                4
# N80                4694329
# L80                6

busco -m genome -i necat_asm_v1.fasta -o necat_busco -l hypocreales_odb10 -f

#assembly with flye using corrected reads from canu
#getting ride of duplicate IDs
seqkit rename fuscarium_v0.0.1.correctedReads.fasta.gz > flye_renamed_fuscarium_v0.0.1.correctedReads.fasta

flye --nano-corr ./flye_renamed_fuscarium_v0.0.1.correctedReads.fasta --out-dir ./flye_asm/ --genome-size 46m -t 24 -i 3 --scaffold

gt seqstat assembly.fasta
# number of contigs:     10
# total contigs length:  42091861
# mean contig size:      4209186.10
# contig size first quartile: 4248396
# median contig size:         4835698
# contig size third quartile: 6506071
# longest contig:             7141792
# shortest contig:            50747
# contigs > 500 nt:           10 (100.00 %)
# contigs > 1K nt:            10 (100.00 %)
# contigs > 10K nt:           10 (100.00 %)
# contigs > 100K nt:          9 (90.00 %)
# contigs > 1M nt:            8 (80.00 %)
# N50                4957983
# L50                4
# N80                4727134
# L80                6

busco -m genome -i ./flye_asm/assembly.fasta -o flye_busco -l hypocreales_odb10

"""
C:94.1%[S:93.9%,D:0.2%],F:2.8%,M:3.1%,n:4494
	4227	Complete BUSCOs (C)
	4220	Complete and single-copy BUSCOs (S)
	7	Complete and duplicated BUSCOs (D)
	124	Fragmented BUSCOs (F)
	143	Missing BUSCOs (M)
	4494	Total BUSCO groups searched
"""

#moving ahead with the Flye asm as it had the best busco scores, longest assembly and shortest contig. N50 and contig number were also very good.

#polishing the consensus sequence using nanopore reads with Nanopolish
nanopolish index -d ./genomes/fuscarium/fast5/run1/ \
-s ./genomes/fuscarium/fast5/run1/hac_called/sequencing_summary.txt \
./genomes/fuscarium/filtered_fuscarium_run1.fastq.gz # for FAST5 inout

nanopolish index -d ./genomes/fuscarium/fast5/run1_restart/ \
-s ./genomes/fuscarium/fast5/run1_restart/hac_called/sequencing_summary.txt \
./genomes/fuscarium/filtered_fuscarium_run1_restart.fastq.gz

#indexing flye asm for polishing
bwa index ./genomes/fuscarium/polishing/flye_assembly.fa

#Align the basecalled reads to the draft sequence
cd ./genomes/fuscarium/polishing/

#left off at this spot and it was working just didn't have time to complete before shutting down desktop

bwa mem -x ont2d -t 32 ./flye_assembly.fasta ./genomes/fuscarium/filtered_fuscarium_run1.fastq.gz | samtools sort -o reads.sorted.bam -T reads.tmp

samtools index reads.sorted.bam

bwa mem -x ont2d -t 32 ./flye_assembly.fasta ./genomes/fuscarium/filtered_fuscarium_run1_restart.fastq.gz | samtools sort -o reads1.sorted.bam -T reads.tmp

samtools index reads1.sorted.bam

samtools merge reads_merged.sorted.bam reads.sorted.bam reads1.sorted.bam

samtools index reads_merged.sorted.bam

#polish assembly NO CHANGES FOUND
python3 nanopolish_makerange.py ./genomes/fuscarium/polishing/flye_assembly.fasta | parallel --results nanopolish.results -P 8 \
    nanopolish variants --consensus -o polished.{1}.vcf -w {1} -r ./genomes/fuscarium/filtered_fuscarium_run1.fastq.gz -b ./genomes/fuscarium/polishing/reads.sorted.bam -g ./genomes/fuscarium/polishing/flye_assembly.fasta -t 4 --min-candidate-frequency 0.1
    nanopolish vcf2fasta -g ./genomes/fuscarium/polishing/flye_assembly.fasta polished.*.vcf > polished_genome.fa
    #evalute polishing
    mkdir analysis
    dnadiff --prefix analysis/polished.dnadiff /./genomes/fuscarium/polishing/flye_assembly.fasta ./polished_genome.fa

##############
# annotation #
##############

#ran a R9.4.1 flowcell with cDNA (PCR free) and basecalled with guppy 6.1.5 using the super high accuracy DNA mode
#merge all the fastq pass files into a consensus transcript file
cd /var/lib/minknow/data/direct_cDNA_fusarium_aven/direct_cDNA_fusarium_aven/20220816_1603_MN37716_FAU08372_a913376f/pass
sudo cat *.gz > ./genomes/fuscarium/rna/LR_transcriptome.fastq.gz

#check for adapter contaimination
fastqc LR_transcriptome.fastq.gz

minimap2 -d flye_asm_sm ../flye_asm.softmask.fasta

paftools.js gff2bed ./genomes/fuscarium/flye_asm_fuscarium/complete_annotation.combined.gtf > anno.bed

minimap2 -ax splice --junc-bed ./anno.bed -t 28 ../flye_asm.softmask.fasta ../../rna/LR_transcriptome.fastq.gz > LR_mapping.sam

samtools view -bS LR_mapping.sam > LR_mapping.bam

samtools sort -o LR_mapping.sort.bam -@ 30 LR_mapping.bam

samtools index LR_mapping.bam

#mapping of short read SRA datasets
#copied mapping done during LoReAn ./genomes/fuscarium/flye_asm_fuscarium/annotation/LoReAn_annotation/run/STAR/short_reads_mapped.bam.sorted.bam

samtools index short_reads_mapped.bam.sorted.bam

"""#MAKER v3 annotation with short read only test
conda activate maker
export LIBDIR=./anaconda3/envs/maker/share/RepeatMasker/Libraries

cd ./genomes/fuscarium/flye_asm_fuscarium/annotation/maker

maker -CTL

maker #after editing the control files appropriately

fasta_merge -d assembly_master_datastore_index.log
gff3_merge -d assembly_master_datastore_index.log

conda deactivate

conda activate busco

busco -m transcriptome -i ./assembly.all.maker.transcripts.fasta -o busco_transcriptome -l hypocreales_odb10

--------------------------------------------------
|Results from dataset hypocreales_odb10           |
--------------------------------------------------
|C:97.9%[S:90.5%,D:7.4%],F:1.0%,M:1.1%,n:4494     |
|4396	Complete BUSCOs (C)                       |
|4065	Complete and single-copy BUSCOs (S)       |
|331	Complete and duplicated BUSCOs (D)        |
|44	Fragmented BUSCOs (F)                     |
|54	Missing BUSCOs (M)                        |
|4494	Total BUSCO groups searched               |
--------------------------------------------------
"""

#Maker using long reads

#map SRA reads to genome soft masked
STAR --runMode genomeGenerate --runThreadN 30 --genomeDir ./genomes/fuscarium/flye_asm_fuscarium/annotation/ --genomeFastaFiles ./genomes/fuscarium/flye_asm_fuscarium/flye_asm.softmask.fasta

STAR --genomeDir ./genomes/fuscarium/flye_asm_fuscarium/annotation/ --outFileNamePrefix SRA --readFilesType Fastx --runThreadN 28 --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --alignIntronMax 2000 --readFilesIn ./genomes/fuscarium/flye_asm_fuscarium/annotation/SRR7962513_1.fastq ./genomes/fuscarium/flye_asm_fuscarium/annotation/SRR7962513_2.fastq

samtools sort -o SRAAligned.out.sort.bam -@ 32 SRAAligned.out.bam

samtools index -@ 32 SRAAligned.out.sort.bam

#map long reads with splicing info
nanoQC -o ./ LR_transcriptome.fastq.gz
NanoPlot -t 20 --fastq LR_transcriptome.fastq
porechop -i LR_transcriptome.fastq.gz -o LR_transcriptome.trimmed.fastq.gz --discard_middle --threads 30

minimap2 -d flye_asm.softmask ./genomes/fuscarium/flye_asm_fuscarium/flye_asm.softmask.fasta

minimap2 -ax splice -k14 -t 32 -G 2000 -uf flye_asm.softmask ./genomes/fuscarium/flye_asm_fuscarium/annotation/LR_transcriptome.trimmed.fastq.gz > long_read_Align.out.sam

samtools view -b long_read_Align.out.sam > long_read_Align.out.bam

samtools sort -o long_read_Align.out.sort.bam -@ 32 long_read_Align.out.bam

samtools index -@ 32 long_read_Align.out.sort.bam

#create transcript models using stringtie
stringtie -o fusarium_stringtie_out.gtf -p 32 --mix SRAAligned.out.sort.bam long_read_Align.out.sort.bam
bascially
gffread -w fusarium_stringtie_out.fa -g ./genomes/fuscarium/flye_asm_fuscarium/flye_asm.softmask.fasta fusarium_stringtie_out.gtf

#rnaMMER fasta output from gensas to quickly reannotate rRNA
gffread -w rnammer.fa -g ./genomes/fuscarium/flye_asm_fuscarium/flye_asm.softmask.fasta ../rnammer.gensas.gff3

conda activate maker
export LIBDIR=/home/chris/anaconda3/envs/maker/share/RepeatMasker/Libraries

mpiexec -n 30 maker -base fusav_rd1 >& log1

mkdir augustus1
cd augustus1
gff3_merge -d ../fusav_rd1.maker.output/fusav_rd1_master_datastore_index.log

## filter gff file, only keep maker annotation in the filtered gff file
awk '{if ($2=="maker") print }' fusav_rd1.all.gff > maker_rnd1.gff

##convert the maker gff and fasta file into a Genbank formated file named pyu.gb
##We keep 2000 bp up- and down-stream of each gene for training the models
gff2gbSmallDNA.pl maker_rnd1.gff ./genomes/fuscarium/flye_asm_fuscarium/assembly.fasta 2000 fusav_rd1.gb

## check number of genes in training set
grep -c LOCUS fusav_rd1.gb

## train model
## first create a new Augustus species named
new_species.pl --species=fusav

## initial training
etraining --species=fusav fusav_rd1.gb
## the initial model should be in the directory
ls -ort $AUGUSTUS_CONFIG_PATH/species/fusav

##create a smaller test set for evaluation before and after optimization. Name the evaluation set pyu.gb.evaluation.
randomSplit.pl fusav_rd1.gb 200
mv fusav_rd1.gb.test fusav_rd1.gb.evaluation

# use the first model to predict the genes in the test set, and check the results
augustus --species=fusav fusav_rd1.gb.evaluation >& first_evaluate.out
grep -A 22 Evaluation first_evaluate.out

"---------------------------------------------\
                 | sensitivity | specificity |
---------------------------------------------|
nucleotide level |       0.944 |       0.803 |
---------------------------------------------/

----------------------------------------------------------------------------------------------------------\
           |  #pred |  #anno |      |    FP = false pos. |    FN = false neg. |             |             |
           | total/ | total/ |   TP |--------------------|--------------------| sensitivity | specificity |
           | unique | unique |      | part | ovlp | wrng | part | ovlp | wrng |             |             |
----------------------------------------------------------------------------------------------------------|
           |        |        |      |                251 |                163 |             |             |
exon level |    584 |    496 |  333 | ------------------ | ------------------ |       0.671 |        0.57 |
           |    584 |    496 |      |  108 |   17 |  126 |  113 |   17 |   33 |             |             |
----------------------------------------------------------------------------------------------------------/

----------------------------------------------------------------------------\
transcript | #pred | #anno |   TP |   FP |   FN | sensitivity | specificity |
----------------------------------------------------------------------------|
gene level |   215 |   200 |   89 |  126 |  111 |       0.445 |       0.414 |
----------------------------------------------------------------------------/
"

# optimize the model. this step is very time consuming. It could take days. To speed things up, you can create a smaller test set
# the following step will create a test and training sets. the test set has 1000
# genes. This test set will be splitted into 24 kfolds for optimization (the kfold
# can be set up to 48, with processed with one cpu core per kfold. Kfold must be
# same number as as cpus). The training, prediction and evaluation will be
# performed on each bucket in parallel (training on hh.gb.train+each bucket, then
# comparing each bucket with the union of the rest). By default, 5 rounds of
# optimization. As optimization for large genome could take days, I changed it to
# 3 here.

randomSplit.pl fusav_rd1.gb 1000
optimize_augustus.pl --species=fusav --kfold=32 --cpus=32 --rounds=5 --onlytrain=fusav_rd1.gb.train fusav_rd1.gb.test >& log &

#train again after optimization
etraining --species=fusav fusav_rd1.gb

augustus --species=fusav fusav_rd1.gb.evaluation >& second_evaluate.out
grep -A 22 Evaluation second_evaluate.out

'*******      Evaluation of gene prediction     *******

---------------------------------------------\
                 | sensitivity | specificity |
---------------------------------------------|
nucleotide level |       0.942 |       0.807 |
---------------------------------------------/

----------------------------------------------------------------------------------------------------------\
           |  #pred |  #anno |      |    FP = false pos. |    FN = false neg. |             |             |
           | total/ | total/ |   TP |--------------------|--------------------| sensitivity | specificity |
           | unique | unique |      | part | ovlp | wrng | part | ovlp | wrng |             |             |
----------------------------------------------------------------------------------------------------------|
           |        |        |      |                242 |                161 |             |             |
exon level |    577 |    496 |  335 | ------------------ | ------------------ |       0.675 |       0.581 |
           |    577 |    496 |      |  105 |   16 |  121 |  111 |   17 |   33 |             |             |
----------------------------------------------------------------------------------------------------------/

----------------------------------------------------------------------------\
transcript | #pred | #anno |   TP |   FP |   FN | sensitivity | specificity |
----------------------------------------------------------------------------|
gene level |   214 |   200 |   92 |  122 |  108 |        0.46 |        0.43 |
----------------------------------------------------------------------------/'

mpiexec -n 30 maker -base fusav_rd2 >& log2

mkdir augustus2
cd augustus2
gff3_merge -d ../fusav_rd2.maker.output/fusav_rd2_master_datastore_index.log

## filter gff file, only keep maker annotation in the filtered gff file
awk '{if ($2=="maker") print }' fusav_rd2.all.gff > maker_rnd2.gff

##convert the maker gff and fasta file into a Genbank formated file named pyu.gb
##We keep 2000 bp up- and down-stream of each gene for training the models
gff2gbSmallDNA.pl maker_rnd2.gff /media/chris/Drive_2/genomes/fuscarium/flye_asm_fuscarium/assembly.fasta 2000 fusav_rd2.gb

## check number of genes in training set
grep -c LOCUS fusav_rd2.gb

## train model
## first create a new Augustus species named
#new_species.pl --species=fusav

## initial training
etraining --species=fusav fusav_rd2.gb
## the initial model should be in the directory
ls -ort $AUGUSTUS_CONFIG_PATH/species/fusav

##create a smaller test set for evaluation before and after optimization. Name the evaluation set pyu.gb.evaluation.
randomSplit.pl fusav_rd2.gb 200
mv fusav_rd2.gb.test fusav_rd2.gb.evaluation

# use the first model to predict the genes in the test set, and check the results
augustus --species=fusav fusav_rd2.gb.evaluation >& first_evaluate.out
grep -A 22 Evaluation first_evaluate.out

'---------------------------------------------\
                 | sensitivity | specificity |
---------------------------------------------|
nucleotide level |       0.973 |        0.96 |
---------------------------------------------/

----------------------------------------------------------------------------------------------------------\
           |  #pred |  #anno |      |    FP = false pos. |    FN = false neg. |             |             |
           | total/ | total/ |   TP |--------------------|--------------------| sensitivity | specificity |
           | unique | unique |      | part | ovlp | wrng | part | ovlp | wrng |             |             |
----------------------------------------------------------------------------------------------------------|
           |        |        |      |                 77 |                 80 |             |             |
exon level |    583 |    586 |  506 | ------------------ | ------------------ |       0.863 |       0.868 |
           |    583 |    586 |      |   48 |    6 |   23 |   48 |    4 |   28 |             |             |
----------------------------------------------------------------------------------------------------------/

----------------------------------------------------------------------------\
transcript | #pred | #anno |   TP |   FP |   FN | sensitivity | specificity |
----------------------------------------------------------------------------|
gene level |   204 |   200 |  151 |   53 |   49 |       0.755 |        0.74 |
----------------------------------------------------------------------------/
'

# optimize the model. this step is very time consuming. It could take days. To speed things up, you can create a smaller test set
# the following step will create a test and training sets. the test set has 1000
# genes. This test set will be splitted into 24 kfolds for optimization (the kfold
# can be set up to 48, with processed with one cpu core per kfold. Kfold must be
# same number as as cpus). The training, prediction and evaluation will be
# performed on each bucket in parallel (training on hh.gb.train+each bucket, then
# comparing each bucket with the union of the rest). By default, 5 rounds of
# optimization. As optimization for large genome could take days, I changed it to
# 3 here.

randomSplit.pl fusav_rd2.gb 1000
optimize_augustus.pl --species=fusav --kfold=32 --cpus=32 --rounds=5 --onlytrain=fusav_rd2.gb.train fusav_rd2.gb.test >& log &

#train again after optimization
etraining --species=fusav fusav_rd2.gb

augustus --species=fusav fusav_rd2.gb.evaluation >& second_evaluate.out
grep -A 22 Evaluation second_evaluate.out

"*******      Evaluation of gene prediction     *******

---------------------------------------------\
                 | sensitivity | specificity |
---------------------------------------------|
nucleotide level |       0.984 |       0.961 |
---------------------------------------------/

----------------------------------------------------------------------------------------------------------\
           |  #pred |  #anno |      |    FP = false pos. |    FN = false neg. |             |             |
           | total/ | total/ |   TP |--------------------|--------------------| sensitivity | specificity |
           | unique | unique |      | part | ovlp | wrng | part | ovlp | wrng |             |             |
----------------------------------------------------------------------------------------------------------|
           |        |        |      |                 75 |                 71 |             |             |
exon level |    590 |    586 |  515 | ------------------ | ------------------ |       0.879 |       0.873 |
           |    590 |    586 |      |   49 |    6 |   20 |   49 |    4 |   18 |             |             |
----------------------------------------------------------------------------------------------------------/

----------------------------------------------------------------------------\
transcript | #pred | #anno |   TP |   FP |   FN | sensitivity | specificity |
----------------------------------------------------------------------------|
gene level |   205 |   200 |  153 |   52 |   47 |       0.765 |       0.746 |
----------------------------------------------------------------------------/"

cp fusav_rd2.all.gff ../

mpiexec -n 32 maker -base fusav_rd3 >& log3

cd ./fusav_rd3.maker.output

fasta_merge -d fusav_rd3_master_datastore_index.log
gff3_merge -n -g -d fusav_rd3_master_datastore_index.log

conda deactivate

conda activate busco

busco -m transcriptome -i ./fusav_rd3.all.maker.transcripts.fasta -o busco_transcriptome -l ./genomes/fuscarium/flye_asm_fuscarium/annotation/maker_final/busco_comparisons/hypocreales_odb10
--------------------------------------------------
|Results from dataset hypocreales_odb10           |
--------------------------------------------------
|C:96.4%[S:93.5%,D:2.9%],F:2.1%,M:1.5%,n:4494     |
|4333	Complete BUSCOs (C)                       |
|4203	Complete and single-copy BUSCOs (S)       |
|130	Complete and duplicated BUSCOs (D)        |
|94	Fragmented BUSCOs (F)                     |
|67	Missing BUSCOs (M)                        |
|4494	Total BUSCO groups searched               |
--------------------------------------------------


busco -m transcriptome -i ./fusav_rd3.all.maker.transcripts.fasta -o busco_transcriptome_sord -l ./fusarium/annnotation/busco_downloads/lineages/sordariomycetes_odb10

--------------------------------------------------
|Results from dataset sordariomycetes_odb10       |
--------------------------------------------------
|C:97.1%[S:94.0%,D:3.1%],F:1.5%,M:1.4%,n:3817     |
|3706	Complete BUSCOs (C)                       |
|3589	Complete and single-copy BUSCOs (S)       |
|117	Complete and duplicated BUSCOs (D)        |
|56	Fragmented BUSCOs (F)                     |
|55	Missing BUSCOs (M)                        |
|3817	Total BUSCO groups searched               |
--------------------------------------------------

agat_sp_statistics.pl --gff ./fusav_rd3.all.gff -o annotation.stats

conda deactivate
conda activate maker

AED_cdf_generator.pl -b 0.025 ./fusav_rd2.all.gff > AED_rnd2
AED_cdf_generator.pl -b 0.025 ./fusav_rd3.all.gff > AED_rnd3

mkdir final_anno
cp fusav_rd3.all.gff ../final_anno/Fusav_v1.0.gff
cp fusav_rd3.all.maker.transcripts.fasta ../final_anno/Fusav_v1.0.transcripts.fa
cp fusav_rd3.all.maker.proteins.fasta ../final_anno/Fusav_v1.0.proteins.fa

#running blastp against uniprot/swissprot
makeblastdb -in uniprot_sprot.fasta -input_type fasta -dbtype prot

blastp -query Fusav_v1.0.proteins.fa -db ./databases/uniprot_sprot.fasta -num_threads 30 -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out output.blastp

#running introproscn
interproscan.sh -cpu 30 -appl pfam -dp -f TSV -goterms -iprlookup -pa -t p -i ./Fusav_v1.0.proteins.fa -o output.iprscan

#correct names of genes
maker_map_ids --prefix Fusav_ --justify 8 --iterate 1 Fusav_v1.0.gff > Fusav_v1.0.id.map

cp Fusav_v1.0.gff Fusav_v1.0.renamed.gff
cp Fusav_v1.0.proteins.fa Fusav_v1.0.proteins.renamed.fasta
cp Fusav_v1.0.transcripts.fa Fusav_v1.0.transcripts.renamed.fasta
cp output.iprscan output.renamed.iprscan
cp output.blastp output.renamed.blastp

map_gff_ids Fusav_v1.0.id.map Fusav_v1.0.renamed.gff
map_fasta_ids Fusav_v1.0.id.map Fusav_v1.0.proteins.renamed.fasta
map_fasta_ids Fusav_v1.0.id.map Fusav_v1.0.transcripts.renamed.fasta
map_data_ids Fusav_v1.0.id.map output.renamed.iprscan
map_data_ids Fusav_v1.0.id.map output.renamed.blastp

maker_functional_gff ./uniprot_sprot.fasta output.renamed.blastp Fusav_v1.0.renamed.gff > Fusav_v1.0.putative_function.gff
maker_functional_fasta ./uniprot_sprot.fasta output.renamed.blastp Fusav_v1.0.proteins.renamed.fasta > Fusav_v1.0.proteins.putative_function.fasta
maker_functional_fasta ./uniprot_sprot.fasta output.renamed.blastp Fusav_v1.0.transcripts.renamed.fasta > Fusav_v1.0.transcripts.putative_function.fasta

ipr_update_gff Fusav_v1.0.putative_function.gff output.renamed.iprscan > Fusav_v1.0.putative_function.domain_added.gff
iprscan2gff3 output.renamed.iprscan Fusav_v1.0.renamed.gff > visible_iprscan_domains.gff

#########################################
# removal of contig 11 as it's circular #
#########################################

seqtk subseq assembly.fasta name.lst > assembly_filt.fa

'name.lst
contig_1
contig_10
contig_3
contig_4
contig_5
contig_6
contig_7
contig_8
contig_9
'

#removal of contig_11 from gff files
grep -v '^contig_11' assembly.all.renamed.putative_function.domain_added.gff > assembly.final_filt.gff
grep -v '^contig_11' rnammer.gensas.gff3 > rnammer_filt.gff3
grep -v '^contig_11' trnascan.gensas.gff3 > trnascan_filt.gff3

gt seqstat assembly_filt.fa
# number of contigs:     9
# total contigs length:  42041114
# mean contig size:      4671234.89
# contig size first quartile: 4727134
# median contig size:         4957983
# contig size third quartile: 6506071
# longest contig:             7141792
# shortest contig:            713051
# contigs > 500 nt:           9 (100.00 %)
# contigs > 1K nt:            9 (100.00 %)
# contigs > 10K nt:           9 (100.00 %)
# contigs > 100K nt:          9 (100.00 %)
# contigs > 1M nt:            8 (88.89 %)
# N50                4957983
# L50                4
# N80                4727134
# L80                6
