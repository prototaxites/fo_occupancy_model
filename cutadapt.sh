#!/bin/zsh

mamba activate cutadapt

mkdir -p outputs
mkdir -p amplicons/Leaf/ITS
mkdir -p amplicons/Stem/ITS
mkdir -p amplicons/Soil/ITS
mkdir -p amplicons/Leaf/16S
mkdir -p amplicons/Stem/16S
mkdir -p amplicons/Soil/16S
mkdir -p amplicons/D2E/16S
mkdir -p amplicons/D2E/ITS
mkdir -p amplicons/Agar/16S
mkdir -p amplicons/Agar/ITS

ulimit -S -n 4096

## Leaf
cutadapt \
    -e 0.15 --no-indels \
	--cores=8 \
    -g ^file:primers/novogene_16s_fwd.fasta \
    -G ^file:primers/novogene_16s_rev.fasta \
    -o amplicons/Leaf/16S/Leaf-16S-{name1}-{name2}.1.fastq.gz -p amplicons/Leaf/16S/Leaf-16S-{name1}-{name2}.2.fastq.gz \
    raw_data/FO16SL_FKDN230312727-1A_HFFM2DRX3_L1_1.fq.gz raw_data/FO16SL_FKDN230312727-1A_HFFM2DRX3_L1_2.fq.gz

echo "file,count" > outputs/leaf_16S_counts.txt
find amplicons/Leaf/ -name "*.1.fastq.gz" -exec sh -c 'echo {},$(echo $(zcat {} |wc -l)/4|bc)' \; >> outputs/leaf_its_counts.txts

cutadapt \
    -e 0.15 --no-indels \
	--cores=8 \
    -g ^file:primers/novogene_its_fwd.fasta \
    -G ^file:primers/novogene_its_rev.fasta \
    -o amplicons/Leaf/ITS/Leaf-ITS-{name1}-{name2}.1.fastq.gz -p amplicons/Leaf/ITS/Leaf-ITS-{name1}-{name2}.2.fastq.gz \
    raw_data/FO_ITS_L_FKDN230027390-1A_HGN7TDRX2_L2_1.fq.gz raw_data/FO_ITS_L_FKDN230027390-1A_HGN7TDRX2_L2_2.fq.gz

echo "file,count" > outputs/leaf_its_counts.txt
find amplicons/Leaf/ -name "*.1.fastq.gz" -exec sh -c 'echo {},$(echo $(zcat {} |wc -l)/4|bc)' \; >> outputs/leaf_its_counts.txt

## Stem
cutadapt \
    -e 0.15 --no-indels \
	--cores=8 \
    -g ^file:primers/novogene_16s_fwd.fasta \
    -G ^file:primers/novogene_16s_rev.fasta \
    -o amplicons/Stem/16S/Stem-16S-{name1}-{name2}.1.fastq.gz -p amplicons/Stem/16S/Stem-16S-{name1}-{name2}.2.fastq.gz \
    raw_data/FO16SS_FKDN230312728-1A_HFFM2DRX3_L1_1.fq.gz raw_data/FO16SS_FKDN230312728-1A_HFFM2DRX3_L1_2.fq.gz

echo "file,count" > outputs/stem_16S_counts.txt
find amplicons/Stem/ -name "*.1.fastq.gz" -exec sh -c 'echo {},$(echo $(zcat {} |wc -l)/4|bc)' \; >> outputs/stem_its_counts.txt

cutadapt \
    -e 0.15 --no-indels \
	--cores=8 \
    -g ^file:primers/novogene_its_fwd.fasta \
    -G ^file:primers/novogene_its_rev.fasta \
    -o amplicons/Stem/ITS/Stem-ITS-{name1}-{name2}.1.fastq.gz -p amplicons/Stem/ITS/Stem-ITS-{name1}-{name2}.2.fastq.gz \
    raw_data/FO_ITS_S_FKDN230027391-1A_HGN7TDRX2_L2_1.fq.gz raw_data/FO_ITS_S_FKDN230027391-1A_HGN7TDRX2_L2_2.fq.gz

echo "file,count" > outputs/stem_its_counts.txt
find amplicons/Stem/ -name "*.1.fastq.gz" -exec sh -c 'echo {},$(echo $(zcat {} |wc -l)/4|bc)' \; >> outputs/stem_its_counts.txt

## Soil
cutadapt \
    -e 0.15 --no-indels \
	--cores=8 \
    -g ^file:primers/novogene_16s_fwd.fasta \
    -G ^file:primers/novogene_16s_rev.fasta \
    -o amplicons/Soil/16S/Soil-16S-{name1}-{name2}.1.fastq.gz -p amplicons/Soil/16S/Soil-16S-{name1}-{name2}.2.fastq.gz \
    raw_data/FO16SR_FKDN230312729-1A_HFFM2DRX3_L1_1.fq.gz raw_data/FO16SR_FKDN230312729-1A_HFFM2DRX3_L1_2.fq.gz

cutadapt \
    -e 0.15 --no-indels \
	--cores=8 \
    -g ^file:primers/novogene_its_fwd.fasta \
    -G ^file:primers/novogene_its_rev.fasta \
    -o amplicons/Soil/ITS/Soil-ITS-{name1}-{name2}.1.fastq.gz -p amplicons/Soil/ITS/Soil-ITS-{name1}-{name2}.2.fastq.gz \
    raw_data/FO_ITS_R_FKDN230027392-1A_HGN7TDRX2_L2_1.fq.gz raw_data/FO_ITS_R_FKDN230027392-1A_HGN7TDRX2_L2_2.fq.gz

echo "file,count" > outputs/soil_its_counts.txt
find amplicons/Soil/ -name "*.1.fastq.gz" -exec sh -c 'echo {},$(echo $(zcat {} |wc -l)/4|bc)' \; >> outputs/soil_its_counts.txt

## D2E

cutadapt \
    -e 0.15 --no-indels \
	--cores=8 \
    -g ^file:primers/novogene_16s_fwd.fasta \
    -G ^file:primers/novogene_16s_rev.fasta \
    -o amplicons/D2E/16S/D2E-16S-{name1}-{name2}.1.fastq.gz -p amplicons/D2E/16S/D2E-16S-{name1}-{name2}.2.fastq.gz \
    raw_data/A1_FKDN220492785-1A_HGN3VDRX2_L2_1.fq.gz raw_data/A1_FKDN220492785-1A_HGN3VDRX2_L2_2.fq.gz

echo "file,count" > outputs/d2e_16s_counts.txt
find amplicons/D2E/ -name "*.1.fastq.gz" -exec sh -c 'echo {},$(echo $(zcat {} |wc -l)/4|bc)' \; >> outputs/d2e_16s_counts.txt

cutadapt \
    -e 0.15 --no-indels \
	--cores=8 \
    -g ^file:primers/novogene_its_fwd.fasta \
    -G ^file:primers/novogene_its_rev.fasta \
    -o amplicons/D2E/ITS/D2E-ITS-{name1}-{name2}.1.fastq.gz -p amplicons/D2E/ITS/D2E-ITS-{name1}-{name2}.2.fastq.gz \
    raw_data/A2_FKDN220492786-1A_HGMTHDRX2_L1_1.fq.gz raw_data/A2_FKDN220492786-1A_HGMTHDRX2_L1_2.fq.gz

echo "file,count" > outputs/d2e_its_counts.txt
find amplicons/D2E/ -name "*.1.fastq.gz" -exec sh -c 'echo {},$(echo $(zcat {} |wc -l)/4|bc)' \; >> outputs/d2e_its_counts.txt

## Agar

cutadapt \
    -e 0.15 --no-indels \
	--cores=8 \
	--action none \
    -g file:primers/16s_its_fwd.fasta \
    -G file:primers/16s_its_rev.fasta \
    -o raw_data/Agar/Agar-{name1}-{name2}.1.fastq.gz -p raw_data/Agar/Agar-{name1}-{name2}.2.fastq.gz \
    raw_data/Metagenome-FastQ_S1_L001_R2_001.fastq.gz raw_data/Metagenome-FastQ_S1_L001_R1_001.fastq.gz

echo "file,count" > outputs/agar_primer_counts.txt
find raw_data/Agar/ -name "*.1.fastq.gz" -exec sh -c 'echo {},$(echo $(zcat {} |wc -l)/4|bc)' \; >> outputs/agar_primer_counts.txt

cutadapt \
    -e 0.15 --no-indels \
	--cores=8 \
    -g ^file:primers/novogene_16s_fwd.fasta \
    -G ^file:primers/novogene_16s_rev.fasta \
    -o amplicons/Agar/16S/Agar-16S-{name1}-{name2}.1.fastq.gz -p amplicons/Agar/16S/Agar-16S-{name1}-{name2}.2.fastq.gz \
    raw_data/Agar/Agar-16S-16S.1.fastq.gz raw_data/Agar/Agar-16S-16S.2.fastq.gz

echo "file,count" > outputs/agar_16s_counts.txt
find amplicons/Agar/16S -name "*.1.fastq.gz" -exec sh -c 'echo {},$(echo $(zcat {} |wc -l)/4|bc)' \; >> outputs/agar_16s_counts.txt

cutadapt \
    -e 0.15 --no-indels \
	--cores=8 \
	-g ^file:primers/novogene_its_fwd.fasta \
	-G ^file:primers/novogene_its_rev.fasta \
    -o amplicons/Agar/ITS/Agar-ITS-{name1}-{name2}.1.fastq.gz -p amplicons/Agar/ITS/Agar-ITS-{name1}-{name2}.2.fastq.gz \
    raw_data/Agar/Agar-ITS-ITS.1.fastq.gz raw_data/Agar/Agar-ITS-ITS.2.fastq.gz

echo "file,count" > outputs/agar_its_counts.txt
find amplicons/Agar/ITS -name "*.1.fastq.gz" -exec sh -c 'echo {},$(echo $(zcat {} |wc -l)/4|bc)' \; >> outputs/agar_its_counts.txt