[![PyPI-Server](https://img.shields.io/pypi/v/firmament.svg)](https://pypi.org/project/firmament/)
![Unit tests](https://github.com/genentech/firmament/actions/workflows/run-tests.yml/badge.svg)
# Firmament

Firmament is a Python package for performing high-throughput gene signature searches on large collections of single-cell RNA-seq datasets stored as H5Ad objects. It helps identify datasets containing sets of cells enriched with cells expressing a gene signature.

The package leverages the [signature score calculation methods](https://github.com/lilab-bcb/pegasus/blob/master/signature_score.pdf) from the [pegasus](https://github.com/lilab-bcb/pegasus) package to compute and search for gene signatures. Firmament's two-phase approach (offline computation + online search) enables efficient signature score searches across 100M+ single cells.

## Installation

```shell
pip install firmament
```

## Usage

### Computing Signature Scores

Easily compute Z-scores for gene signatures across your single-cell data:

```python
import anndata
import firmament

# Load your data
ad = anndata.read_h5ad("path/to/your/data.h5ad")

# Compute Z-scores
zscore_df = firmament.compute.calc_zscores(
    ad, 
    layer=None,  # Use .X matrix (or specify a layer name)
    cell_type_labels="level1class"  # Column in adata.obs containing cell type labels
)

print(zscore_df)
```

    ### output
                Tspan12     Tshz1    Fnbp1l  Adamts15    Cldn12     Rxfp1  \
    1772071015_C02 -0.996770 -0.162126 -0.171727 -0.454830 -0.370999 -0.424958   
    1772071017_G12 -0.773850 -0.321869 -0.658392 -0.516157 -0.128007 -0.273857   
    1772071017_A05 -0.830232 -0.829437  0.106758 -0.402295 -0.438613 -0.577686   
    1772071014_B06  0.046787 -0.065567 -0.399475 -0.278929 -0.757392 -0.670489   
    1772067065_H06 -1.004179 -0.033353 -0.709410 -0.278839 -0.722943 -0.405029   
    ...                  ...       ...       ...       ...       ...       ...   

One can imagine a system stores these computed signature Z-score matrices to rapidly identify cell types or individual cells enriched for the expression of specific genes across large data collections.

## Searching for Gene Signature Enrichment

Search for cells or cell types enriched for specific gene signatures:

```python
from firmament import signature_search

# Search for cells enriched for a set of genes
results = signature_search(
    "path/to/your/data.h5ad",
    genes=["Fnbp1l", "Tspan12", "Vipas39"],
    layer_name=None,  # Use .X matrix (or specify a layer name)
    cell_type_label_column_name="level1class",
    verbose=True
)

print(results)  # Can also be converted into a Pandas DataFrame
```

    ### output

        fdr_prop	total_count	fdr_count	zscore_hist	pval_hist	mean_zscore	set_zscore	set_pval	label
    0	0.001997	3005	6	{'bins': ['-3.20', '-3.14', '-3.09', '-3.03', ...	{'bins': ['0.01', '0.02', '0.03', '0.04', '0.0...	-0.932011	-51.090881	1.0	overall


### Batch Processing

For larger collections of files, you can use the `batch_signature_search` function to perform searches _on-the-fly_ across multiple datasets.

<!-- pyscaffold-notes -->

## Note

This project has been set up using [BiocSetup](https://github.com/biocpy/biocsetup)
and [PyScaffold](https://pyscaffold.org/).
