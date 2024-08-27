# Contributing guide

## Installing dev dependencies

In addition to the packages needed to _use_ this package, you need additional python packages to _run tests_ and _build
the documentation_. It's easy to install them using `pip`:

```bash
cd CAT
pip install -e ".[dev,test,doc]"
```

## Code-style

This package uses [pre-commit][] to enforce consistent code-styles.
On every commit, pre-commit checks will either automatically fix issues with the code, or raise an error message.

To enable pre-commit locally, simply run

```bash
pre-commit install
```

in the root of the repository. Pre-commit will automatically download all dependencies when it is run for the first time.

Alternatively, you can rely on the [pre-commit.ci][] service enabled on GitHub. If you didn't run `pre-commit` before
pushing changes to GitHub it will automatically commit fixes to your pull request, or show an error message.

If pre-commit.ci added a commit on a branch you still have been working on locally, simply use

```bash
git pull --rebase
```

to integrate the changes into yours.
While the [pre-commit.ci][] is useful, we strongly encourage installing and running pre-commit locally first to understand its usage.

Finally, most editors have an _autoformat on save_ feature. Consider enabling this option for [ruff][ruff-editors]
and [prettier][prettier-editors].

[ruff-editors]: https://docs.astral.sh/ruff/integrations/
[prettier-editors]: https://prettier.io/docs/en/editors.html

## Publishing a release

### Updating the version number

Before making a release, you need to update the version number in the `pyproject.toml` file. Please adhere to [Semantic Versioning][semver], in brief

> Given a version number MAJOR.MINOR.PATCH, increment the:
>
> 1. MAJOR version when you make incompatible API changes,
> 2. MINOR version when you add functionality in a backwards compatible manner, and
> 3. PATCH version when you make backwards compatible bug fixes.
>
> Additional labels for pre-release and build metadata are available as extensions to the MAJOR.MINOR.PATCH format.

Once you are done, commit and push your changes and navigate to the "Releases" page of this project on GitHub.
Specify `vX.X.X` as a tag name and create a release. For more information, see [managing GitHub releases][]. This will automatically create a git tag and trigger a Github workflow that creates a release on PyPI.

## Writing documentation

Please write documentation for new or changed features and use-cases. This project uses [sphinx][] with the following features:

-   the [myst][] extension allows to write documentation in markdown/Markedly Structured Text
-   [Numpy-style docstrings][numpydoc] (through the [napoloen][numpydoc-napoleon] extension).
-   Jupyter notebooks as tutorials through [myst-nb][] (See [Tutorials with myst-nb](#tutorials-with-myst-nb-and-jupyter-notebooks))
-   [Sphinx autodoc typehints][], to automatically reference annotated input and output types
-   Citations (like {cite:p}`Rothova2022`) can be included with [sphinxcontrib-bibtex](https://sphinxcontrib-bibtex.readthedocs.io/)

See the [scanpy developer docs](https://scanpy.readthedocs.io/en/latest/dev/documentation.html) for more information
on how to write documentation.

### Tutorials with myst-nb and jupyter notebooks

The documentation is set-up to render jupyter notebooks stored in the `docs/notebooks` directory using [myst-nb][].
Currently, only notebooks in `.ipynb` format are supported that will be included with both their input and output cells.
It is your responsibility to update and re-run the notebook whenever necessary.

If you are interested in automatically running notebooks as part of the continuous integration, please check
out [this feature request](https://github.com/scverse/cookiecutter-scverse/issues/40) in the `cookiecutter-scverse`
repository.

#### Hints

-   If you refer to objects from other packages, please add an entry to `intersphinx_mapping` in `docs/conf.py`. Only
    if you do so can sphinx automatically create a link to the external documentation.
-   If building the documentation fails because of a missing link that is outside your control, you can add an entry to
    the `nitpick_ignore` list in `docs/conf.py`

#### Building the docs locally

```bash
cd docs
make html
open _build/html/index.html
```

<!-- Links -->

[pre-commit.ci]: https://pre-commit.ci/
[myst-nb]: https://myst-nb.readthedocs.io/en/latest/
[pre-commit]: https://pre-commit.com/
[semver]: https://semver.org/
[sphinx]: https://www.sphinx-doc.org/en/master/
[myst]: https://myst-parser.readthedocs.io/en/latest/intro.html
[numpydoc-napoleon]: https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html
[numpydoc]: https://numpydoc.readthedocs.io/en/latest/format.html
[sphinx autodoc typehints]: https://github.com/tox-dev/sphinx-autodoc-typehints
[managing GitHub releases]: https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository
