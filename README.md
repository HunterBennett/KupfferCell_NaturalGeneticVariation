# Python and R Code from Kupffer Cell Natural Genetic Variation Paper (Bennett et al., 2023)
____
In this repository you can find all the code that was used for analysis and visualization <a href="https://www.biorxiv.org/content/10.1101/2022.09.22.509046v1">Bennett et al</a>.

Briefly, this paper explores the effect of natural genetic variation *in vivo*, using Kupffer cells as a model cell type. We collect and analyze transcriptional and epigenetic data (ATAC-seq, H3K27Ac ChIP-seq) to identify putative *trans* regulators driving differential gene expression across inbred strains of mice. Additionally, we provide evidence that *trans* effects control a majority of strain differential genes at homeostasis while *cis* effects dominate the transcriptional response to an external signal (lipopolysaccharide).

### Content of this repository
____
We have organized our scripts by the figures they produce within the manuscript.
* figure1
    * figure1_assembly.R: R script to assemble figure 1
    * fn_clear_session.R: R function to clear session
    * fn_figure1_assembly.R: R function called by assembly script to build fig 1
* Figure 2
    * Figure2_ATAC_00_Annotation.ipynb: Generates peak counts at IDR peaks and provides basic QC of ATAC-seq data.
    * Figure2_ATAC_01_UpSetPlot.ipynb: Generates UpSetPlot found in Figure 2.
    * Figure2_ATAC_02_AccessibleMotifs.ipynb: Perform motif analysis in accessible sites across the three strains.
    * Figure2_ATAC_03_Differential.ipynb: Identify differentially accessible peaks across the three strains.
    * Figure2_ATAC_04_Differential_Peak_Mutation_Burden.ipynb: Identify mutation rate across differentially accessible peaks.
    * Figure2_ATAC_05_Homer.ipynb: Perform motif analysis across differentially accessible peaks.
    * Figure2_ATAC_06_Maggie.ipynb: Perform motive *mutation* analysis across differentially accesssible peaks.
    * Figure2_H3K27Ac_00_Annotation.ipynb: Generates peak counts at acetylated regions and provides basic QC of H3K27Ac ChIP-seq data.
    * Figure2_H3K27Ac_01_UpSetPlot.ipynb: Generates UpSetPlot found in Figure 2.
    * Figure2_H3K27Ac_02_AccessibleMotifs.ipynb: Perform motif analysis in acetylated sites across the three strains.
    * Figure2_H3K27Ac_03_Differential.ipynb: Identify differentially acetylated peaks across the three strains.
    * Figure2_H3K27Ac_04_Differential_Peak_Mutation_Burden.ipynb: Identify mutation rate across differentially acetylated peaks.
    * Figure2_H3K27Ac_05_Homer.ipynb: Perform motif analysis across differentially acetylated peaks.
    * Figure2_H3K27Ac_06_Maggie.ipynb: Perform motive *mutation* analysis across differentially acetylated peaks.   
* Figure 3
    * Figure3_00_NicheNet_Data_Prep.ipynb: Download a prepare data used in NicheNet analysis.
    * Figure3_01_AJ.R: Perform NicheNet analysis for A/J differentially expressed genes.
    * Figure3_02_BALBcJ.R Perform NicheNet analysis for BALB/cJ differentially expressed genes.
    * Figure3_03_C57BL6J.R: Perform NicheNet analysis for C57BL/6J differentially expressed genes.
    * Figure3_04_Circos.R: Generate circos plot found in Figure 3.
* Figure 4
    * Figure4_00_CB6F1J_NSG_MA_Plots.ipynb: Generate MA plots found in Figure 4.
    * Figure4_01_CB6F1J_CisTrans_Plot.ipynb: Generate *cis*/*trans* gene expression plots for CB6F1/J mice found in Figure 4.
    * Figure4_02_CB6F1J_NSG_CisTrans_Plot.ipynb: Generate *cis*/*trans* gene expression plots for F1NSG mice found in Figure 4.
    * Figure4_03_NSG_CisTrans_Plot: Generate *cis*/*trans* gene expression plots for F0-NSG mice found in Figure 4.
    * Figure4_04_UpSetPlot.ipynb: Generate UpSetPlot found in Figure 4.
    * Figure4_05_CisTrans_Effectsize_CumulativeDist.ipynb: Generate plot of cumulative effect size of *cis* and *trans* acting mutations in CB6F1/J mice.
    * Figure4_06_Cell_Extrinsic_Intrinsic_Heatmaps.ipynb: Create heatmaps illustrating cell intrinsic vs. extrinsic patterns of gene expression.
    * Figure4_07_Promoter_Motifs_and_Gene_Ontology.ipynb: Perform gene ontology and promoter motif analysis for *cis* and *trans* gene sets.
* Figure 5
    * Figure5_ATAC_00_Annotation_Github.ipynb: Generates allellic peak counts at IDR peaks and performs basic principal component analysis.
    * Figure5_ATAC_01_Differential_Peaks.ipynb: Identifies differentially accessible ATAC-seq peaks in CB6F1/J mice.
    * Figure5_ATAC_02_CisTrans.ipynb: Identifies *cis* and *trans* regulated ATAC-seq peaks in parental strains using CB6F1/J ATAC-seq data.
    * Figure5_ATAC_03_CisTrans_KnownMotifs.ipynb: Motif analysis of *cis* and *trans* regulated ATAC-seq peaks.
    * Figure5_ATAC_04_CisTrans_Promoter_Motifs.ipynb: Motif analysis of accessible enhancers associated with *trans* regulated genes.
    * Figure5_H3K27Ac_00_Annotation_Github.ipynb: Generates allellic peak counts at IDR peaks and performs basic principal component analysis.
    * Figure5_H3K27Ac_01_Differential_Peaks.ipynb: Identifies differentially accessible ATAC-seq peaks in CB6F1/J mice.
    * Figure5_H3K27Ac_02_CisTrans.ipynb: Identifies *cis* and *trans* regulated ATAC-seq peaks in parental strains using CB6F1/J ATAC-seq data.
    * Figure5_H3K27Ac_03_CisTrans_KnownMotifs.ipynb: Motif analysis of *cis* and *trans* regulated ATAC-seq peaks.
    * Figure5_H3K27Ac_CisTrans_Promoter_Motifs.ipynb: Motif analysis of accessible enhancers associated with *trans* regulated genes.
* Figure 6
    * figure6_assembly.R: R script to build most of figure 6
    * fn_clear_session.R: R function to clear session
    * fn_process_peaks_strict.R: R function to process peak files
    * fn_process_plot.R: R function to process fig 6 plot
    * Figure6_00_Maggie_LowEqualHighBasal_Github.ipynb: MAGGIE analysis of low/equal/high basal ATAC-seq peaks from LPS treatment experiment.  
* Extended Data Figure 1
    * extended_figure1.R: script to build extended data figure 1
* Extended Data Figure 2
    * extended_figure2.R: script to build extended data figure 2
    * fn_clear_session.R: R function to clear session
    * pheatmap2.R: function to build heatmap for figure 2
    
If you encounter a problem with the provided code please post in the <a href="https://github.com/HunterBennett/KupfferCell_NaturalGeneticVariation/issues">issues section</a>.  

### Packages utlized in this work
____

This work utilizes the following published algorithms:
* <a href="https://github.com/zeyang-shen/maggie">MAGGIE</a> Shen et al., Bioinformatics 2020.
* <a href="https://github.com/saeyslab/nichenetr">NicheNet</a> Browaeys et al., Nature Methods 2019.
* <a href="https://github.com/vlink/marge">MMARGE</a> Link et al., Bioinformatics, 2018.
* <a href="http://homer.ucsd.edu/homer/">HOMER</a> Heinz et al., Molecular Cell, 2010.
* <a href="https://github.com/s-andrews/FastQC">FastQC</a> Andrews et al.
* <a href="https://github.com/alexdobin/STAR">STAR</a> Dobin et al., Bioinformatics 2013.
* <a href="https://github.com/BenLangmead/bowtie2">BowTie2</a> Langmead et al., Nature Methods 2012.
* <a href="https://bioconductor.org/packages/release/bioc/html/DESeq2.html">DESeq2</a> Love et al., Genome Biology 2014.
* <a href="https://github.com/z5ouyang/A.C.Rsuite">A.C.RSuite</a> Zhengyu Ouyang.

### Contributors
___
Code was written by Hunter Bennett and Ty Troutman.

### References
___
Hunter Bennett, Ty D. Troutman, Enchen Zhou, Nathanael J. Spann, Verena M. Link, Jason S. Seidman, Christian K. Nickl, Yohei Abe, Mashito Sakai, Martina P. Pasillas, Justin M. Marlman, Carlos Guzman, Mojgan Hosseini, Bernd Schnabl, Christopher K. Glass
bioRxiv 2022.09.22.509046; doi: https://doi.org/10.1101/2022.09.22.509046
