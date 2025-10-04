# 🧬 Shiny Adjustable PCA Plot for Proteomics Data

![R Version](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue)
![Shiny](https://img.shields.io/badge/Shiny-1.7.0-green)
![Docker](https://img.shields.io/badge/Docker-Ready-brightgreen)
![License](https://img.shields.io/badge/License-MIT-yellow)

## 🎯 Overview

An interactive **Shiny web application** for real-time Principal Component Analysis (PCA) of proteomics datasets. This tool provides researchers with dynamic exploration of protein expression data through customizable PCA visualization, specifically designed for comparing Wild Type (WT) and Knockout (KO) samples.

![PCA Proteomics Analysis Interface](pca_proteomics.png)

*Figure: Interactive PCA visualization showing WT vs KO sample clustering with explained variance*

---

## 🔬 Application Purpose

This Shiny application is designed for:
- **Sample Quality Control**: Identify outliers and batch effects in proteomics experiments
- **Group Separation**: Visualize differences between experimental conditions (WT vs KO)
- **Data Exploration**: Interactively adjust normalization and filtering parameters
- **Reproducible Analysis**: Consistent PCA computation with adjustable parameters

---

## ✨ Key Features

### Interactive Controls

| Control | Description | Options/Range |
|---------|-------------|---------------|
| **Run PCA Button** | Trigger PCA analysis with current settings | Click to execute |
| **Normalization Method** | Choose data normalization approach | Voom or Log CPM |
| **Min Expression Threshold** | Filter low-expression proteins | Default: 0.5 (adjustable) |
| **Low Variance Filter** | Toggle filtering of low-variance genes | On/Off checkbox |
| **Variance Quantile Slider** | Set threshold for variance filtering | 0-1 (default: 0.25) |

### Visualization Output

- **PCA Biplot**: PC1 vs PC2 scatter plot with:
  - Color-coded groups (WT: blue, KO: red)
  - Semi-transparent polygons for group boundaries
  - Sample labels with automatic positioning (ggrepel)
  - Explained variance percentages on axes
  
- **PCA Summary**: Statistical summary of principal components

---

## 🚀 Quick Start

### Option 1: Docker Deployment

```bash
# Clone repository
git clone https://github.com/carlosbuss1/SHINY_Adjustable_PCA_plot_4_Proteomics.git
cd SHINY_Adjustable_PCA_plot_4_Proteomics

# Build Docker image
docker build -t shiny-pca-app .

# Run container
docker run -p 3838:3838 shiny-pca-app

# Access in browser
# http://localhost:3838
```

### Option 2: Local R Environment

```r
# Install required packages
install.packages(c("shiny", "ggplot2", "ggrepel", "limma", "edgeR"))

# Clone and run
setwd("path/to/SHINY_Adjustable_PCA_plot_4_Proteomics")
shiny::runApp("shiny_pca_plot.R")
```

---

## 📁 Repository Structure

```
SHINY_Adjustable_PCA_plot_4_Proteomics/
├── shiny_pca_plot.R          # Main Shiny application
├── Dockerfile                # Container specification
├── proteomics_counts.csv     # Input data file (tab-delimited)
├── pca_proteomics.png        # Example output screenshot
└── README.md                 # This documentation
```

---

## 📋 Data Format Requirements

### Input File: `proteomics_counts.csv`

- **Format**: Tab-delimited text file (`.csv` extension but tab-separated)
- **Structure**: 
  - First column: `Gene.Symbol` (protein/gene identifiers)
  - Subsequent columns: Sample data (KO1-KO5, WT1-WT5)
- **Values**: Raw count data or intensity values

### Example Structure

```
Gene.Symbol	KO1	KO2	KO3	KO4	KO5	WT1	WT2	WT3	WT4	WT5
ACTB	15420	18903	12456	20145	17890	22456	19872	25631	21098	23456
GAPDH	8934	7652	9871	6543	8765	10234	11456	9876	12345	10987
TP53	23456	19872	25631	21098	18765	15432	17890	14567	16789	15234
```

---

## ⚙️ Processing Pipeline

### 1. Data Loading
- Reads tab-delimited `proteomics_counts.csv`
- Maps sample names to standardized format (KO1-5, WT1-5)

### 2. Data Preprocessing
- Replaces missing values with column minimum
- Creates experimental design matrix (WT vs KO)

### 3. Normalization Options

| Method | Description | Use Case |
|--------|-------------|----------|
| **Voom** | Variance modeling at the observation level | RNA-seq style count data with mean-variance relationship |
| **Log CPM** | Log2 counts per million | General normalization for count data |

### 4. Filtering Steps
- **Expression filtering**: Removes genes below mean expression threshold
- **Variance filtering**: Removes bottom quantile of low-variance genes (optional)

### 5. PCA Computation
- Performs scaling and centering
- Calculates principal components
- Computes explained variance percentages

---

## 🎨 Visual Design

### Color Scheme
- **WT (Wild Type)**: 
  - Points: Royal Blue (`royalblue`)
  - Polygon fill: Light Blue (`lightblue`)
- **KO (Knockout)**: 
  - Points: Red (`red3`)
  - Polygon fill: Pink (`pink`)

### Plot Elements
- Point size: 4 with 0.8 alpha transparency
- Text labels: Size 3, black color
- Theme: Minimal with base size 15
- Title: Centered, size 18

---

## 🔧 Parameter Guidelines

### Recommended Settings

```r
# For high-quality proteomics data
normalization: "Voom"
min_expression: 0.5
filter_low: TRUE
low_variation_threshold: 0.25

# For noisy or sparse data
normalization: "Log CPM"
min_expression: 1.0
filter_low: TRUE
low_variation_threshold: 0.5
```

---

## 🐛 Troubleshooting

### Common Issues

#### Data Not Loading
```r
# Check file format - must be tab-delimited
data <- read.csv("proteomics_counts.csv", sep = "\t")
```

#### PCA Not Running
- Click "Run PCA" button after adjusting parameters
- Ensure data has sufficient non-zero values

#### Memory Issues
```r
# Increase Shiny upload size for large files
options(shiny.maxRequestSize = 30*1024^2)  # 30MB
```

---

## 📊 Sample Requirements

- **Minimum samples**: 3 per group (ideally 5+ for robust PCA)
- **Minimum features**: 100+ proteins/genes after filtering
- **Data completeness**: Handle missing values before upload or use built-in imputation

---

## 📚 Dependencies

### R Packages
- `shiny` - Web application framework
- `ggplot2` - Advanced plotting
- `ggrepel` - Text label positioning
- `limma` - Linear models for microarray/RNA-seq
- `edgeR` - Differential expression analysis

### Version Requirements
- R ≥ 4.0.0
- All packages should be latest CRAN/Bioconductor versions

---

## 🤝 Contributing

Contributions welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Test with multiple datasets
4. Submit a pull request

---

## 📜 License

MIT License - See `LICENSE` file for details

---

## 👨‍🔬 Author

**Carlos E. Buss, PhD**  
Bioinformatics Researcher  
Signal Transduction and Metabolism Laboratory  
Université libre de Bruxelles (ULB)  
Brussels, Belgium  

📧 Email: carlos.eduardo.buss@ulb.be  
🌐 Lab: [www.stmlaboratory.com](https://www.stmlaboratory.com)  
💻 GitHub: [@carlosbuss1](https://github.com/carlosbuss1)  

---

## 🙏 Acknowledgments

- [Bioconductor](https://bioconductor.org/) for `limma` and `edgeR` packages
- [RStudio](https://www.rstudio.com/) for Shiny framework
- Proteomics community for feedback

---

*Last updated: January 2025 | Version: 1.0.0*
