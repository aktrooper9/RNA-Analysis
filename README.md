
# RNA-Seq Analysis Pipeline

##  Folder Structure
Each folder represents a **sample set**. Inside each folder you will find:
- **Raw data**  
- **Matrix file**  
- **Different analyses** performed  
- **R scripts** used to run the analyses  

---

##  Running the R Scripts
To run the R scripts properly:

1. Download the raw data and matrix files into the **same directory** as your R project.  
2. Run the import_pandas_as_pd.py script in your terminal: 


## Special Instructions for Sample Set `135092`

For dataset **135092**, additional preprocessing steps are required:

1. Run **merge.py** to combine all the samples.
2. Use **135092\_combine\_geneid.py** to merge based on gene IDs.

---

##  Plot Summaries

### Cnetplot

A **network plot** that visualizes the relationships between enriched biological terms (e.g., pathways, functions) and the genes associated with them.

* Nodes represent genes and enriched terms.
* Edges (connections) show which genes contribute to which terms.
  This plot is helpful for identifying shared genes across multiple pathways.

---

### DEG Counts Plot

A **bar plot** that displays the number of **differentially expressed genes (DEGs)**.

* Separate bars usually represent upregulated vs. downregulated genes.
  This provides a quick summary of how many genes show significant expression changes under different conditions.

---

### Heatmap

A **heatmap** of gene expression values across samples.

* Rows represent genes.
* Columns represent samples.
* Colors correspond to expression levels (e.g., high expression = red, low expression = blue).
  Heatmaps allow clustering of genes and samples to reveal expression patterns or groups of samples with similar profiles.

---

### KEGG Dotplot

A **dotplot** that summarizes enrichment analysis results for KEGG pathways.

* Each dot corresponds to a pathway.
* Dot size represents the number of genes involved in that pathway.
* Dot color indicates the statistical significance (e.g., adjusted p-value).
  This makes it easy to see which pathways are most enriched and biologically relevant.

---

### MA Plot

A **scatter plot** comparing expression levels across conditions.

* X-axis = mean expression of each gene.
* Y-axis = log fold-change.
* Genes with significant changes are highlighted (e.g., red for upregulated, blue for downregulated).
  The MA plot helps distinguish overall gene expression distribution from significantly altered genes.

---

### PCA Plot

A **Principal Component Analysis (PCA) plot** that reduces data dimensionality and clusters samples.

* Each point represents a sample.
* Samples closer together have more similar gene expression profiles.
* PCA helps reveal natural groupings, batch effects, or treatment differences.

---

### Reactome Barplot

A **barplot** displaying enriched pathways from the Reactome database.

* Each bar corresponds to a pathway.
* Bar length/height represents either the number of genes in the pathway or the enrichment score.
  This plot highlights which biological processes or pathways are significantly involved.

---

### Volcano Plot

A **scatter plot** combining statistical significance and magnitude of change.

* X-axis = log fold-change (effect size).
* Y-axis = -log10(p-value) (significance).
* Genes with large fold-changes and strong statistical significance appear in the “volcano arms.”
  The volcano plot is a quick way to identify the most biologically important and statistically robust DEGs.

---





---

---
---

##  Key Results

### Major Findings in Males
The differential expression analysis for **male samples** (mgs level 4 vs 1) revealed that:  
- **Gene ID 4057** was **far more expressed** than any other gene.  
- The **top DEGs plot** showed a single dominant signal (gene 4057), while the rest of the genes showed minimal or no differential expression.  

This result is unusual because in typical RNA-seq analyses, we expect to see a distribution of multiple genes that are upregulated or downregulated.  

### Comparison with Females
In contrast, the **female samples** showed:  
- Multiple genes differentially expressed.  
- A more typical distribution, with several genes both upregulated and downregulated across conditions.  

### Why is This Strange?
The male-specific result suggests that:  
1. **One gene (4057) drives nearly all the observed differences** in expression.  
2. This could be due to:  
   - **Biological reasons** → 4057 may have a sex-specific regulatory role.  
   - **Technical artifact** → sample imbalance, mapping issues, or errors in raw counts.  

Because females showed multiple DEGs, the **male-specific dominance of one gene** stands out as biologically intriguing but also calls for careful QC and validation.  

---

##  Supporting Figures

### Top DEGs (Males)  
<img src="115828/Male 1 vs 4/results/plots/top_degs_plot_mgs_1_vs_4_males.png" alt="Top DEGs Males" width="600"/>  

### Top DEGs (Females)  
<img src="115828/Female%201%20vs%204/results/plots/top_degs_plot_mgs_1_vs_4_females.png?raw=true" alt="Top DEGs Females" width="600"/>  

These visualizations highlight the difference between males (single dominant DEG) and females (multiple DEGs).  
