# CAT

[![Tests][badge-tests]][link-tests]
[![Documentation][badge-docs]][link-docs]

[badge-tests]: https://img.shields.io/github/actions/workflow/status/matq007/CAT/test.yaml?branch=main
[link-tests]: https://github.com/brickmanlab/CAT/actions/workflows/test.yml
[badge-docs]: https://img.shields.io/readthedocs/CAT

Cluster Alignment Tool (CAT)

## Testing

```bash
catcli \
    --ds1 ./tests/datasets/mock.h5ad \
    --ds1_name DS1 \
    --ds1_cluster Condition_E+D \
    --ds2 ./tests/datasets/mock.h5ad \
    --ds2_name DS2 \
    --ds2_cluster Condition_E+D \
    --output ./res
```

## Getting started

Please refer to the [documentation][link-docs]. In particular, the

-   [API documentation][link-api].

## Installation

You need to have Python 3.10 or newer installed on your system. If you don't have
Python installed, we recommend installing [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge).

There are several alternative options to install CAT:

1. Install the latest release of `cat-python` from [PyPI][link-pypi]:

```bash
pip install cat-python
```

1. Install the latest development version:

```bash
pip install git+https://github.com/brickmanlab/CAT.git@main
```

## Release notes

See the [changelog][changelog].

## Contact

If you found a bug, please use the [issue tracker][issue-tracker].

## Citation

> Rothová, M.M., Nielsen, A.V., Proks, M. et al. Identification of the central intermediate in the extra-embryonic to embryonic endoderm transition through single-cell transcriptomics. Nat Cell Biol 24, 833–844 (2022). https://doi.org/10.1038/s41556-022-00923-x

[issue-tracker]: https://github.com/matq007/CAT/issues
[changelog]: https://CAT.readthedocs.io/latest/changelog.html
[link-docs]: https://CAT.readthedocs.io
[link-api]: https://CAT.readthedocs.io/latest/api.html
[link-pypi]: https://pypi.org/project/CAT
