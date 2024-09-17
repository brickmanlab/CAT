# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog][],
and this project adheres to [Semantic Versioning][].

[keep a changelog]: https://keepachangelog.com/en/1.0.0/
[semantic versioning]: https://semver.org/spec/v2.0.0.html

## [Unreleased]

## [v1.1.0]

> [!WARNING]
> **Breaking changes**
> New changes to the API, please see the documentation

### Added

- Build check with automatic Pypi release

### Changed

- New API documentation
- Switch to `uv` as package manager
- Sankey plots are generated with new Python code, removed the old R script
- Silence `RuntimeWarning` message when performing median normalization

## [v1.0.1]

### Fixed

- Error in saving XlsxWriter [#2](https://github.com/brickmanlab/cat-python/issues/2)

## [v1.0]

Initial release for publication
