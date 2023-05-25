# PBMC uG single-cell RNA-seq MTD validation
# 1nd: Gammaretrovirus
# extract the Gammaretrovirus reads
mkdir ~/scRNAseq_analysis/Gammaretrovirus_reads
cd ~/scRNAseq_analysis/Gammaretrovirus_reads
# female
for i in $lsn; do
    python ~/MTD/Tools/KrakenTools/extract_kraken_reads.py \
        -k /bigrock/FurmanLab/Fei/scRNAseq_runs/scRNA_hPBMCmicrogravity_050622/MTD_F/output/temp/Report_non-host_raw_${i}.kraken \
        -s1 /bigrock/FurmanLab/Fei/scRNAseq_runs/scRNA_hPBMCmicrogravity_050622/MTD_F/output/temp/${i}_non-host_raw.fq \
        -o ~/scRNAseq_analysis/Gammaretrovirus_reads/${i}_153135.fq \
        -r /bigrock/FurmanLab/Fei/scRNAseq_runs/scRNA_hPBMCmicrogravity_050622/MTD_F/output/temp/Report_non-host.raw_${i}.txt \
        --fastq-output \
        --taxid 153135 --include-children
done
# male
for i in $lsn; do
    python ~/MTD/Tools/KrakenTools/extract_kraken_reads.py \
        -k /bigrock/FurmanLab/Fei/scRNAseq_runs/scRNA_hPBMCmicrogravity_052522/MTD_M/output/temp/Report_non-host_raw_${i}.kraken \
        -s1 /bigrock/FurmanLab/Fei/scRNAseq_runs/scRNA_hPBMCmicrogravity_052522/MTD_M/output/temp/${i}_non-host_raw.fq \
        -o ~/scRNAseq_analysis/Gammaretrovirus_reads/${i}_153135.fq \
        -r /bigrock/FurmanLab/Fei/scRNAseq_runs/scRNA_hPBMCmicrogravity_052522/MTD_M/output/temp/Report_non-host.raw_${i}.txt \
        --fastq-output \
        --taxid 153135 --include-children
done

# run Magic-blast alignment for extracted Gammaretrovirus reads
# Create a Gammaretrovirus index for Magic-BLAST
makeblastdb -in ~/MTD/hisat2_index_Gammaretrovirus/genome.fa -dbtype nucl -parse_seqids -out ~/MTD/Gammaretrovirus_blastdb/Gammaretrovirus_blastdb
# female
for i in $lsn; do
    magicblast -query ~/scRNAseq_analysis/Gammaretrovirus_reads/Female/${i}_153135.fq \
            -db ~/MTD/Gammaretrovirus_blastdb/Gammaretrovirus_blastdb \
            -infmt fastq \
            -out ~/scRNAseq_analysis/Gammaretrovirus_reads/Female/${i}_153135_blast.txt \
            -outfmt tabular \
            -num_threads 20
done
# count the unique mapped reads
for i in $lsn; do
    echo "mapped reads in sample: "$i
    awk 'BEGIN {min_identity=90} ($1 ~ /^[a-zA-Z]/ && $3 >= min_identity) {unique_reads[$1]=1} END {for (read_id in unique_reads) {count++} print count}' ${i}_153135_blast.txt
done

# male
for i in $lsn; do
    magicblast -query ~/scRNAseq_analysis/Gammaretrovirus_reads/Male/${i}_153135.fq \
            -db ~/MTD/Gammaretrovirus_blastdb/Gammaretrovirus_blastdb \
            -infmt fastq \
            -out ~/scRNAseq_analysis/Gammaretrovirus_reads/Male/${i}_153135_blast.txt \
            -outfmt tabular \
            -num_threads 20
done
# count the unique mapped reads
for i in $lsn; do
    echo "mapped reads in sample: "$i
    awk 'BEGIN {min_identity=90} ($1 ~ /^[a-zA-Z]/ && $3 >= min_identity) {unique_reads[$1]=1} END {for (read_id in unique_reads) {count++} print count}' ${i}_153135_blast.txt
done

# 2nd: Mycobacterium canettii
# concatinate fna.gz files of the Mycobacterium canettii to build a overall genome
for f in *.fna.gz; do zcat "$f"; done > genome.fa
# Create a Mycobacterium canettii index for Magic-BLAST
makeblastdb -in ~/MTD/Mycobacterium_canettii_blastdb/genome.fa -dbtype nucl -parse_seqids -out ~/MTD/Mycobacterium_canettii_blastdb/Mycobacterium_canettii_blastdb

# extract the Mycobacterium canettii reads
# female
for i in $lsn; do
    python ~/MTD/Tools/KrakenTools/extract_kraken_reads.py \
        -k /bigrock/FurmanLab/Fei/scRNAseq_runs/scRNA_hPBMCmicrogravity_050622/MTD_F/output/temp/Report_non-host_raw_${i}.kraken \
        -s1 /bigrock/FurmanLab/Fei/scRNAseq_runs/scRNA_hPBMCmicrogravity_050622/MTD_F/output/temp/${i}_non-host_raw.fq \
        -o ~/scRNAseq_analysis/Mycobacterium_canettii_reads/Female/${i}_78331.fq \
        -r /bigrock/FurmanLab/Fei/scRNAseq_runs/scRNA_hPBMCmicrogravity_050622/MTD_F/output/temp/Report_non-host.raw_${i}.txt \
        --fastq-output \
        --taxid 78331 --include-children
done
# male
for i in $lsn; do
    python ~/MTD/Tools/KrakenTools/extract_kraken_reads.py \
        -k /bigrock/FurmanLab/Fei/scRNAseq_runs/scRNA_hPBMCmicrogravity_052522/MTD_M/output/temp/Report_non-host_raw_${i}.kraken \
        -s1 /bigrock/FurmanLab/Fei/scRNAseq_runs/scRNA_hPBMCmicrogravity_052522/MTD_M/output/temp/${i}_non-host_raw.fq \
        -o ~/scRNAseq_analysis/Mycobacterium_canettii_reads/Male/${i}_78331.fq \
        -r /bigrock/FurmanLab/Fei/scRNAseq_runs/scRNA_hPBMCmicrogravity_052522/MTD_M/output/temp/Report_non-host.raw_${i}.txt \
        --fastq-output \
        --taxid 78331 --include-children
done
# run Magic-blast alignment for extracted Mycobacterium canettii reads
# female
cd ~/scRNAseq_analysis/Mycobacterium_canettii_reads/Female
for i in $lsn; do
    magicblast -query ~/scRNAseq_analysis/Mycobacterium_canettii_reads/Female/${i}_78331.fq \
            -db ~/MTD/Mycobacterium_canettii_blastdb/Mycobacterium_canettii_blastdb \
            -infmt fastq \
            -out ~/scRNAseq_analysis/Mycobacterium_canettii_reads/Female/${i}_78331_blast.txt \
            -outfmt tabular \
            -num_threads 20
done
# count the unique mapped reads
for i in $lsn; do
    echo "mapped reads in sample: "$i
    awk 'BEGIN {min_identity=90} ($1 ~ /^[a-zA-Z]/ && $3 >= min_identity) {unique_reads[$1]=1} END {for (read_id in unique_reads) {count++} print count}' ${i}_78331_blast.txt
done

# male
cd ~/scRNAseq_analysis/Mycobacterium_canettii_reads/Male
for i in $lsn; do
    magicblast -query ~/scRNAseq_analysis/Mycobacterium_canettii_reads/Male/${i}_78331.fq \
            -db ~/MTD/Mycobacterium_canettii_blastdb/Mycobacterium_canettii_blastdb \
            -infmt fastq \
            -out ~/scRNAseq_analysis/Mycobacterium_canettii_reads/Male/${i}_78331_blast.txt \
            -outfmt tabular \
            -num_threads 20
done
# count the unique mapped reads
for i in $lsn; do
    echo "mapped reads in sample: "$i
    awk 'BEGIN {min_identity=90} ($1 ~ /^[a-zA-Z]/ && $3 >= min_identity) {unique_reads[$1]=1} END {for (read_id in unique_reads) {count++} print count}' ${i}_78331_blast.txt
done

