# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
### ...
- ...



## [2.0.1]
### Changed
- Re-packed software for pypi.

### Added
- Default parameter values.



## [1.4.5]
### Fixed
- Fixed issue with new interaction type labeling system.

### Added
- Automatic R packages installation (argparser and ggplot2).

### Changed
- Separated curve calculation and fraction calculation to different functions.
- Harmonized `melt_X.py` scripts output.
- Better fasta reading management: not kept in memory and considering sequences with same header as separate fasta items.

### Added
- `plot_melt_curves_coupled.R` for single-oligo coupled melting curve plot.
- Allowed for initial formamide and sodium concentration for melting temperature correction.



## [1.4.4]
### Changed
- Renamed `oligomeltlib.py` to `meltlib.py`.
- Renamed `oligomelt.py` to `melt_duplex.py`.

### Added
- `melt_second.py` to perform formamide correction and melting curve calculation for secondary structures predicted with OligoArrayAux.

## [1.4.3]
### Added
- Melting curve plotting script (from fish-conditions repo).

### Changed
- Moved functions to separate library

### Fixed
- Now melting curves show the proper (inverted) dissociation curve.



## [1.4.2]
### Changed
- Now using user-defined formamide m-value for wright correction.



## [1.4.0]
### Added
- New feature: formamide correction.



## [1.3.0]
### Added
- New feature: temperature curve calculation. Proper fasta input.



## [1.2.2]
### Fixed
- Allawi and Freier thermodynamic tables.



## [1.2.1]
### Fixed
- Introduced Sugimoto (DNA:RNA) thermodynamic table.



## [1.2.0]
### Added
- New feature: DNA/RNA and RNA/DNA duplex calculation.



## [1.1.0]
### Added
- New input file mode.



### Fixed
- Mg2+ correction now skips Na+ correction.



## [1.0.0]



* [unreleased] https://github.com/ggirelli/oligo-melting
* [2.0.1] https://github.com/ggirelli/oligo-melting/releases/tag/v2.0.1
* [2.0.0] https://github.com/ggirelli/oligo-melting/releases/tag/v2.0.0
* [1.4.5] https://github.com/ggirelli/oligo-melting/releases/tag/v1.4.5
* [1.4.4] https://github.com/ggirelli/oligo-melting/releases/tag/v1.4.4
* [1.4.3] https://github.com/ggirelli/oligo-melting/releases/tag/v1.4.3
* [1.4.2] https://github.com/ggirelli/oligo-melting/releases/tag/v1.4.2
* [1.4.0] https://github.com/ggirelli/oligo-melting
* [1.3.0] https://github.com/ggirelli/oligo-melting
* [1.2.2] https://github.com/ggirelli/oligo-melting
* [1.2.1] https://github.com/ggirelli/oligo-melting
* [1.2.0] https://github.com/ggirelli/oligo-melting
* [1.1.0] https://github.com/ggirelli/oligo-melting
* [1.0.0] https://github.com/ggirelli/oligo-melting
