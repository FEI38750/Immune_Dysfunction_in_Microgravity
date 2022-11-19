# --- 10x Cellranger pipeline --- #
# Get matrix for seprate samples for female
cd ~/scRNAseq_runs/scRNA_hPBMCmicrogravity_050622/
cellranger count --id=Control_1G_stimulated \
                   --transcriptome=/home/buckcenter.org/fwu/cellranger/refdata-gex-GRCh38-2020-A \
                   --fastqs=/bigrock/Winer_Furman/scRNA_hPBMCmicrogravity_050622/HC5GFDRX2/outs/fastq_path/HC5GFDRX2 \
                   --sample=Control_1G_stimulated \
                   --include-introns \
                   --localcores=48

cellranger count --id=Microgravity_uG_simulated \
                   --transcriptome=/home/buckcenter.org/fwu/cellranger/refdata-gex-GRCh38-2020-A \
                   --fastqs=/bigrock/Winer_Furman/scRNA_hPBMCmicrogravity_050622/HC5GFDRX2/outs/fastq_path/HC5GFDRX2 \
                   --sample=Microgravity_uG_simulated \
                   --include-introns \
                   --localcores=48

cellranger count --id=Control_1G_un-stimulated \
                   --transcriptome=/home/buckcenter.org/fwu/cellranger/refdata-gex-GRCh38-2020-A \
                   --fastqs=/bigrock/Winer_Furman/scRNA_hPBMCmicrogravity_050622/HC5GFDRX2/outs/fastq_path/HC5GFDRX2 \
                   --sample=Control_1G_un-stimulated \
                   --include-introns \
                   --localcores=48

cellranger count --id=Microgravity_uG_un-simulated \
                   --transcriptome=/home/buckcenter.org/fwu/cellranger/refdata-gex-GRCh38-2020-A \
                   --fastqs=/bigrock/Winer_Furman/scRNA_hPBMCmicrogravity_050622/HC5GFDRX2/outs/fastq_path/HC5GFDRX2 \
                   --sample=Microgravity_uG_un-simulated \
                   --include-introns \
                   --localcores=48

# Get matrix for seprate samples for male
cd ~/scRNAseq_runs/scRNA_hPBMCmicrogravity_052522/
cellranger count --id=Control_1G_stimulated \
                   --transcriptome=/home/buckcenter.org/fwu/cellranger/refdata-gex-GRCh38-2020-A \
                   --fastqs=/bigrock/Winer_Furman/scRNA_hPBMCmicrogravity_052522/HGW2KDRX2/outs/fastq_path/HGW2KDRX2 \
                   --sample=Control-1G-stimulated \
                   --include-introns \
                   --localcores=18

cellranger count --id=Control_1G_unstimulated \
                   --transcriptome=/home/buckcenter.org/fwu/cellranger/refdata-gex-GRCh38-2020-A \
                   --fastqs=/bigrock/Winer_Furman/scRNA_hPBMCmicrogravity_052522/HGW2KDRX2/outs/fastq_path/HGW2KDRX2 \
                   --sample=Control-1G-un-stimulated \
                   --include-introns \
                   --localcores=18

cellranger count --id=Microgravity_uG_stimulated \
                   --transcriptome=/home/buckcenter.org/fwu/cellranger/refdata-gex-GRCh38-2020-A \
                   --fastqs=/bigrock/Winer_Furman/scRNA_hPBMCmicrogravity_052522/HGW2KDRX2/outs/fastq_path/HGW2KDRX2 \
                   --sample=Microgravity-uG-simulated \
                   --include-introns \
                   --localcores=18

cellranger count --id=Microgravity_uG_unstimulated \
                   --transcriptome=/home/buckcenter.org/fwu/cellranger/refdata-gex-GRCh38-2020-A \
                   --fastqs=/bigrock/Winer_Furman/scRNA_hPBMCmicrogravity_052522/HGW2KDRX2/outs/fastq_path/HGW2KDRX2 \
                   --sample=Microgravity-uG-un-simulated \
                   --include-introns \
                   --localcores=18

# --- MTD piepleine --- #
bash ~/MTD/MTD_singleCell.sh -i ~/scRNAseq_analysis/MTD_F/samplesheet_SC.csv -o ~/scRNAseq_analysis/MTD_F/output -h 9606 -t 20 -p 3 -d 5
bash ~/MTD/MTD_singleCell.sh -i ~/scRNAseq_analysis/MTD_M/samplesheet_SC.csv -o ~/scRNAseq_analysis/MTD_M/output -h 9606 -t 20 -p 3 -d 5