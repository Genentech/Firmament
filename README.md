<!-- These are examples of badges you might want to add to your README:
     please update the URLs accordingly

[![Built Status](https://api.cirrus-ci.com/github/<USER>/firmament.svg?branch=main)](https://cirrus-ci.com/github/<USER>/firmament)
[![ReadTheDocs](https://readthedocs.org/projects/firmament/badge/?version=latest)](https://firmament.readthedocs.io/en/stable/)
[![Coveralls](https://img.shields.io/coveralls/github/<USER>/firmament/main.svg)](https://coveralls.io/r/<USER>/firmament)
[![PyPI-Server](https://img.shields.io/pypi/v/firmament.svg)](https://pypi.org/project/firmament/)
[![Conda-Forge](https://img.shields.io/conda/vn/conda-forge/firmament.svg)](https://anaconda.org/conda-forge/firmament)
[![Monthly Downloads](https://pepy.tech/badge/firmament/month)](https://pepy.tech/project/firmament)
[![Twitter](https://img.shields.io/twitter/url/http/shields.io.svg?style=social&label=Twitter)](https://twitter.com/firmament)
-->

[![Project generated with PyScaffold](https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold)](https://pyscaffold.org/)

# firmament

Provides method to perform gene signature search on a collection of anndata objects.

## Installation

```shell
pip install firmament
```

## Usage

```python
from firmament import signature_search
import os

results = signature_search(os.getcwd() +  "/tests/data/test.h5ad",
    genes=["HLA-F"], verbose=True
)

print(results) # can also be converted into a Pandas DataFrame
```

<!-- pyscaffold-notes -->

## Note

This project has been set up using [BiocSetup](https://github.com/biocpy/biocsetup)
and [PyScaffold](https://pyscaffold.org/).
