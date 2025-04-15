import os

from firmament.search import signature_search

__author__ = "kancherj"
__copyright__ = "kancherj"
__license__ = "MIT"


def test_sigsearch():
    """Test signature search on local H5AD."""
    results = signature_search(
        os.getcwd() + "/tests/data/test.h5ad",
        genes=["Tspan12", "Vipas39", "Cldn12"],
        layer_name=None,
        cell_type_label_column_name="level1class",
        include_celltype_enrichment=True,
        verbose=True,
    )

    assert isinstance(results, list)
    assert isinstance(results[0], dict)
