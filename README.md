# epiCausalCellNet

**Granger Causality-based Cell-Cell Communication Network Inference in Alzheimer's Disease**

[![Python](https://img.shields.io/badge/Python-3.12-blue)](https://python.org)
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)

## Overview

A novel computational framework that identifies **causal** cell-cell 
communication networks from single-nucleus RNA sequencing (snRNA-seq) 
data using Granger causality analysis.

Unlike existing tools (dsCellNet, CellPhoneDB, NATMI), this framework 
establishes **directional causality** between interacting cell types — 
identifying which cell is the causal sender and which is the receiver.

## Key Features

- Granger Causality testing for directional causal inference
- Permutation validation (1,000 iterations) for robustness
- Identifies AD-specific interactions absent in healthy controls
- Cross-validation in independent datasets
- Compatible with any snRNA-seq dataset

## Installation

```bash
pip install scanpy statsmodels gseapy scipy pandas numpy matplotlib
```

## Quick Start

```python
from epiCausalCellNet import load_data, quality_control, granger_causality

# Load data (GEO: GSE174367)
adata = load_data('GSE174367.h5', 'cell_meta.csv')

# Quality control
adata = quality_control(adata)

# Run Granger Causality
results = granger_causality(
    adata,
    lr_pairs  = [('CCL2','CCR2'), ('APOE','LRP1'), ('CX3CL1','CX3CR1')],
    ct_pairs  = [('ASC','MG'), ('ODC','MG'), ('MG','EX')],
    diagnosis_val = 'AD'
)
print(results[results['Causal'] == 'YES'])
```

## Key Results

| Sender | Receiver | Ligand | Receptor | Perm p |
|--------|----------|--------|----------|--------|
| ASC | MG | CCL2 | CCR2 | 0.001 |
| ODC | ASC | APOE | LRP1 | 0.002 |
| ASC | MG | NRG1 | ERBB4 | 0.008 |
| ODC | MG | PROS1 | AXL | 0.008 |
| MG | EX | CX3CL1 | CX3CR1 | 0.048 |

- **10** AD-specific causal interactions identified
- **0** causal interactions in healthy controls
- **9/10** absent in normal brain development (GSE153164)

## Data Availability

- Primary dataset: [GEO GSE174367](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174367)
- Validation dataset: [GEO GSE153164](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153164)

## Repository Structure# epiCausalCellNet
Granger Causality-based Cell-Cell Communication in Alzheimer's Disease
## Comparison with Existing Tools

| Feature | dsCellNet | CellPhoneDB | NATMI | epiCausalCellNet |
|---------|-----------|-------------|-------|-----------------|
| Causal inference | No | No | No | **Yes** |
| Time series | Yes | No | No | Yes |
| Permutation validation | No | Yes | No | **Yes** |
| AD-specific (validated) | N/A | N/A | N/A | **10** |

## Citation

If you use this tool, please cite:## License

MIT License — free to use and modify with attribution.

## Contact

jamilasma13122024@gmail.com
