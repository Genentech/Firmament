import logging
from typing import Union

import numpy as np
import pandas as pd
import pegasus
from anndata import AnnData
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests as padjust

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def calc_zscores(
    data: AnnData,
    gene_symbols: str = None,
    celltype_annotation: str = "cellTypeOntologyID",
    verbose: bool = False,
) -> pd.DataFrame:
    """Calculate z-scores on an ``AnnData`` Object.

    Friendly wrapper around `calculate_z_score` method.

    Args:
        data:
            Input ``AnnData`` object.

        gene_symbols:
            Column name from `var` that contains gene symbols.
            Defaults to None.

        celltype_annotation:
            Column name from `obs` that specifies cell type labels. Defaults to None.

    Returns:
        A :py:class:`~pandas.DataFrame` of z-scores with each column as gene.
    """
    logging.info("preprocessing data...")

    if gene_symbols is not None:
        if verbose:
            logging.info(f"setting column {gene_symbols} as index for `data.var`")
        data.var = data.var.set_index(gene_symbols, drop=False)
    else:
        if "symbol" in data.var.columns:
            if verbose:
                logging.info("dataset contains gene symbols; setting this as index to `data.var`")
            data.var = data.var.set_index("symbol", drop=False)

    if verbose:
        logging.info("converting X to float32")
    data.X = data.X.astype(np.float32)

    if verbose:
        logging.info("calculating z-scores using pegasus")

    z_mat = pegasus.calculate_z_score(data)
    z_mat_df = pd.DataFrame(data=z_mat, columns=data.var.index.values)
    z_mat_df = z_mat_df.apply(pd.to_numeric, errors="coerce", downcast="float")

    if celltype_annotation is not None and celltype_annotation in data.obs.columns:
        if verbose:
            logging.info(f"using {celltype_annotation} as celltype annotation")
        z_mat_df["cell_types"] = data.obs[celltype_annotation].values
    else:
        if verbose:
            logging.info("Dataset does not contains cell type labels, using 'NA' as default.")
        z_mat_df["cell_types"] = "NA"

    z_mat_df.index = pd.Series(data.obs.index.values)
    # welp somethings gotta go because symbols are not unique
    z_mat_df = z_mat_df.loc[:, ~z_mat_df.columns.duplicated()]
    return z_mat_df


def compute_hist(
    x: Union[list, pd.Series],
    num_bins: int = 100,
    verbose: bool = False,
) -> dict:
    """Compute a histogram for a vector.

    Args:
        x:
            Input list or pandas Series object.

        num_bins:
            Number of bins.
            Defaults to 100.

    Returns:
        A dict containing the bins and the frequencies.
    """
    if verbose:
        logging.info("computing histogram bins")

    hist = np.histogram(x, bins=num_bins)
    res = {}
    res["bins"] = [f"{x:.2f}" for x in hist[1]]
    res["freq"] = hist[0]
    return res


def calc_stats(
    zscores: Union[list, pd.Series],
    pvals: Union[list, pd.Series],
    num_hist_bins: int = 100,
    alpha: int = 0.1,
    verbose: bool = False,
) -> dict:
    """Calculate stats to for the dataset.

    Args:
        zscores:
            Mean zscore vector for genes of interest.

        pvals:
            Calculated pvalues.

        num_hist_bins:
            Number of histogram bins to generate. Defaults to 100.

        alpha:
            Alpha for error rate. Defaults to 0.1.

    Returns:
        A Dictionary containing the stats.
    """
    if verbose:
        logging.info("computing stats..")

    fdr_res = padjust(pvals, method="fdr_bh", alpha=alpha)
    mean_zscore = np.mean(zscores)
    set_zscore = mean_zscore * np.sqrt(len(zscores))
    set_pval = 1.0 - norm.cdf(set_zscore)
    zscoreGt2 = [z for z in zscores if z > 2]

    res = {
        "fdr_prop": len(zscoreGt2) / len(zscores),  # np.mean(fdr_res[0]) # % of  zscore > 2
        "total_count": len(fdr_res[0]),
        "fdr_count": len(zscoreGt2),  # np.sum(fdr_res[0]), # of  zscore > 2
        "zscore_hist": compute_hist(zscores, num_bins=num_hist_bins, verbose=verbose),
        "pval_hist": compute_hist(pvals, num_bins=num_hist_bins, verbose=verbose),
        "mean_zscore": mean_zscore,
        "set_zscore": set_zscore,
        "set_pval": set_pval,
    }

    return res
