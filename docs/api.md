# API

Import CAT as:

```python
import cat
```

```{eval-rst}
.. currentmodule:: cat
```

## CLI

```{eval-rst}
.. currentmodule:: cat
.. autosummary::
    :toctree: generated

    cli.init_logger
    cli.parse_args
    cli.main
```

## CAT

```{eval-rst}
.. currentmodule:: cat
.. autosummary::
    :toctree: generated

    cat.internal_preprocessing
    cat.compare
    cat.run
```

## Dataset

```{eval-rst}
.. currentmodule:: cat
.. autosummary::
    :toctree: generated

    Dataset._fix_metadata
    Dataset._fix_genes
    Dataset._fix_genes
    Dataset._filter_genes
    Dataset._save
    Dataset.prepare
```

## DatasetDiff

```{eval-rst}
.. currentmodule:: cat
.. autosummary::
    :toctree: generated
```

## Report

```{eval-rst}
.. currentmodule:: cat
.. autosummary::
    :toctree: generated

    report.generate_tables
    report.to_excel
    report.save_tables
```

## Utils

```{eval-rst}
.. currentmodule:: cat
.. autosummary::
    :toctree: generated

    utils.get_nz_mean
    utils.get_nz_median
    utils.normalize
    utils.rename_ds
    utils.read_features
```
