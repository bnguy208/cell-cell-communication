# **CellChat Overview**

CellChat is a computational tool for analyzing and visualizing cell-cell communication based on single-cell RNA sequencing (scRNA-seq) data. It leverages known ligand-receptor interactions to infer signaling networks, enabling researchers to study cellular communication patterns in health and disease.

---

## **Key Features**
1. **Database of Ligand-Receptor Pairs**
   - Curated from KEGG and primary literature.
   - Includes secreted signaling, ECM-receptor interactions, and cell-cell contact signaling.

2. **Modeling Cell-Cell Communication**
   - Quantifies communication probability using the **law of mass action**.
   - Identifies statistically and biologically significant interactions.

3. **Visualization Tools**
   - **Hierarchy Plot**: Directional communication flow between cell groups.
   - **Circle Plot**: Connectivity between cell groups.
   - **Bubble Plot**: Strength of ligand-receptor interactions across groups.

4. **Advanced Analysis**
   - **Network Centrality Analysis**: Identifies sender, receiver, mediator, and influencer roles.
   - **Signaling Pathway Clustering**: Groups pathways by functional similarity.
   - **Comparative Analysis**: Identifies shared and context-specific signaling patterns.

---

## **Workflow**

### **Step 1: Data Input and Preprocessing**
- Input: scRNA-seq gene expression matrix (genes vs. cells).
- Options:
  - Use predefined cell groups.
  - Infer relationships via **latent space** or **latent trajectory** analysis.

### **Step 2: Cellular Communication Modeling**
- Identify overexpressed genes in each cell group.
- Calculate communication probabilities based on ligand-receptor interactions.

### **Step 3: Visualization**
- Generate **hierarchy plots**, **circle plots**, and **bubble plots** for communication networks.

### **Step 4: Advanced Analysis**
- Perform **network centrality analysis** to determine cell roles.
- Compare signaling across conditions or datasets.

---

## **Installation**
Install CellChat using R:
```R
# Install CellChat
install.packages("devtools")
devtools::install_github("sqjin/CellChat")
