# Cluster Alignment Tool (CAT)

[![PyPI][pypi-badge]][pypi-link]
[![build][build-badge]][build-link]
[![Documentation][docs-badge]][docs-link]

## Installation

```console
$ pip install cat-python
```

## Installation from source

```console
$ pip install git+https://github.com/brickmanlab/cat-python.git@master
```

## Running CAT

```console
$ catcli \
    --ds1 ./tests/datasets/mock.h5ad \
    --ds1_name DS1 \
    --ds1_cluster Condition_E+D \
    --ds2 ./tests/datasets/mock.h5ad \
    --ds2_name DS2 \
    --ds2_cluster Condition_E+D \
    --output ./res
```

## Build documentation

```console
$ sphinx-build -M html docs docs/_build
```

## Citation

Please consider citing scANVI Explainer if you use in your research.

> Rothová, M.M., Nielsen, A.V., Proks, M. et al. <br>
> Identification of the central intermediate in the extra-embryonic to embryonic endoderm transition through single-cell transcriptomics. <br>
> Nat Cell Biol 24, 833–844 (2022). [10.1038/s41556-022-00923-x]

```BibTeX
@article{rothova2022,
  title = {Identification of the Central Intermediate in the Extra-Embryonic to Embryonic Endoderm Transition through Single-Cell Transcriptomics},
  author = {Rothov{\'a}, Michaela Mrugala and Nielsen, Alexander Valentin and Proks, Martin and Wong, Yan Fung and Riveiro, Alba Redo and {Linneberg-Agerholm}, Madeleine and David, Eyal and Amit, Ido and Trusina, Ala and Brickman, Joshua Mark},
  year = {2022},
  month = jun,
  journal = {Nature Cell Biology},
  volume = {24},
  number = {6},
  pages = {833--844},
  publisher = {Nature Publishing Group},
  issn = {1476-4679},
  doi = {10.1038/s41556-022-00923-x}
}

```

[pypi-badge]: https://img.shields.io/pypi/v/cat-python.svg
[pypi-link]: https://pypi.org/project/cat-python
[docs-badge]: https://readthedocs.org/projects/brickmanlabcat/badge/?version=latest
[docs-link]: https://brickmanlabcat.readthedocs.io/en/latest/
[build-badge]: https://github.com/brickmanlab/cat-python/actions/workflows/build.yml/badge.svg
[build-link]: https://github.com/brickmanlab/cat-python/actions/workflows/build.yml
[10.1038/s41556-022-00923-x]: https://doi.org/10.1038/s41556-022-00923-x
