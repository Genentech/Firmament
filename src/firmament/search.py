import logging
from typing import List, Union

import anndata
import pandas
from scipy.stats import norm

from .compute import calculate_stats, compute_zscores

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def signature_search(
    adata_or_path: Union[str, anndata.AnnData],
    genes: list,
    layer_name: str,
    gene_symbols_column_name: str = None,
    cell_type_label_column_name: str = None,
    num_hist_bins: int = 100,
    alpha: int = 0.1,
    include_celltype_enrichment: bool = False,
    verbose: bool = False,
) -> List[dict]:
    """Compute signature scores on a local dataset.

    Args:
        adata_or_path:
            Path to a H5ad file or :py:class:`~anndata.AnnData` object.

        genes:
            Input genes or features of interest.

        layer_name:
            Name of the layer containing the count matrix.

            Use None to use `X` slot of the AnnData object.

            This is expected to be consistent across all the anndata objects.

        gene_symbols_column_name:
            Column in the ``var`` slot of ``AnnData`` containing gene symbols.
            Defaults to the index of the ``var`` slot.

        cell_type_label_column_name:
            Column in the ``obs`` slof of ``AnnData`` containing cell type labels.
            Defaults to the index of the ``obs`` slot.

        num_hist_bins:
            Number of histrogram bins.
            Defaults to 100.

        alpha:
            Alpha for family wise error rate.
            Defaults to 0.1.

        include_celltype_enrichment:
            Wehether to group cells by celltype and compute enrichment scores.
            Must provide a ``celltype_label_column_name`` containing celltype labels.
            Defaults to False.

        verbose:
            Whether to display logs.
            Defaults to False.

    Returns:
        list containing results; each element in the list is a dictionary containng
        the cell type label, enrichment scores and other computed statistics.
    """
    results = []

    try:
        data = adata_or_path
        if isinstance(adata_or_path, str):
            data = anndata.read_h5ad(adata_or_path)

        z_mats = compute_zscores(
            data,
            gene_symbols=gene_symbols_column_name,
            cell_type_labels=cell_type_label_column_name,
            layer_name=layer_name,
            verbose=verbose,
        )

        missing_genes = list(set(genes).difference(set(z_mats.columns)))

        if len(missing_genes) > 0:
            msg = f"Input Dataset does not contain full set of input genes: {missing_genes}."
            if verbose:
                logging.error(msg)
            raise ValueError(msg)

        z_mat = z_mats[genes]
        zscores = z_mat.mean(axis=1)
        pvals = 1.0 - norm.cdf(zscores)
        res = calculate_stats(zscores, pvals, num_hist_bins=num_hist_bins, alpha=alpha)
        res["label"] = "overall"
        results.append(res)

        if include_celltype_enrichment:
            celltypes = pandas.DataFrame({"cell_types": z_mats["cell_types"], "zscores": zscores, "pvals": pvals})

            groups = celltypes.groupby("cell_types")

            if len(groups) > 1:
                for cname, group in groups:
                    res = calculate_stats(
                        group["zscores"].values,
                        group["pvals"].values,
                        num_hist_bins=num_hist_bins,
                        alpha=alpha,
                    )
                    res["label"] = cname
                    results.append(res)
    except Exception as e:
        logging.error(str(e))
        results = None
    finally:
        return results


def batch_signature_search(
    adatas_or_paths: List[Union[str, anndata.AnnData]],
    genes: List[str],
    gene_symbols_column_name: str = None,
    celltype_label_column_name: str = None,
    num_hist_bins: int = 100,
    alpha: float = 0.1,
    include_celltype_enrichment: bool = False,
    verbose: bool = False,
):
    """Compute signature scores on a collection of datasets (Batch mode).

    Args:
        adatas_or_paths:
            List of paths to a H5ad file or :py:class:`~anndata.AnnData` objects.

        genes:
            Input genes or features of interest.

        gene_symbols_column_name:
            Column in the ``var`` slot of ``AnnData`` containing gene symbols.
            Defaults to the index of the ``var`` slot.

        celltype_label_column_name:
            Column in the ``obs`` slof of ``AnnData`` containing cell type labels.
            Defaults to the index of the ``obs`` slot.

        num_hist_bins:
            Number of histrogram bins.
            Defaults to 100.

        alpha:
            Alpha for family wise error rate.
            Defaults to 0.1.

        include_celltype_enrichment:
            Wehether to group cells by celltype and compute enrichment scores.
            Must provide a ``celltype_label_column_name`` containing celltype labels.
            Defaults to False.

        verbose:
            Whether to display logs.
            Defaults to False.

    Returns:
        list containing results; each element in the list is a dictionary containng
        the cell type label, enrichment scores and other computed statistics.
    """
    results = []

    for path in adatas_or_paths:
        if verbose:
            logging.info(f"processing dataset: {path}")
        res = signature_search(
            path,
            genes=genes,
            gene_symbols_column_name=gene_symbols_column_name,
            celltype_label_column_name=celltype_label_column_name,
            num_hist_bins=num_hist_bins,
            alpha=alpha,
            include_celltype_enrichment=include_celltype_enrichment,
            verbose=verbose,
        )
        if res is not None:
            results.extend(res)

    return results
