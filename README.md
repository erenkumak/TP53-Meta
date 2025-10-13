# TP53-META, a meta-analysis tool for comparative transcriptomics of TP53 dependency

## Background
TP53 is the most frequently mutated transcription factor (TF) in sporadic cancers; and its targets exhibit dysregulation at the level of expression in both cancer and non-cancer pathologies. However, there is not yet a web-based tool that enables meta-analysis and visualization of TP53-related gene expression datasets, although several databases exist to access and annotate TP53 target information.  

To address this gap, we developed **TP53-META**, an interactive **R Shiny-based web application** that allows users to upload and simultaneously analyze user-supplied or integrated public RNA-seq datasets for effects of TP53 depletion and/or induction on the transcriptome.

---

## Results and Functionality
**TP53-META** enables:
- Visualization of **significant expression clusters**
- Construction of **TF-TF**, **pathway-pathway**, **disease-gene**, and **treatment-gene** networks
- Exploration of **TP53-dependency** of selected treatment contrasts
- **User data upload** and **meta-analysis** of public datasets  

Example public instance: [http://konulabapps.bilkent.edu.tr:3838/TP53-Meta1.5/](http://konulabapps.bilkent.edu.tr:3838/TP53-Meta1.5/)

## 🐳 Running TP53-META with Docker

This repository includes a `Dockerfile` generated with shiny2docker.  
It allows anyone to reproduce the Shiny app environment without manual dependency installation.

### 1. Clone the repository
```bash
git clone https://github.com/erenkumak/TP53-Meta.git
cd TP53-Meta
docker build -t tp53-meta .
# Expose the container’s Shiny port (3838) to your host:
docker run -d -p 3838:3838 --name tp53_meta tp53-meta
# Prefer port 80?
docker run -d -p 80:3838 --name tp53_meta tp53-meta
# Verify it’s running
docker ps
docker logs -f tp53_meta
# To Stop
docker stop tp53_meta
