# Change Log

## [1.2.10] [WIP] - 2024-05-xx
 
### Added

### Changed
- Updated user guide

### Fixed


## [1.2.9] - 2024-05-02
 
### Added
- miRNA specific preprocessing, QC and plots both in backend and GUI.
- Possibility to generate svg files and HQ image resolution.
- miRNA example dataset
- Evaluation visualization option to choose plot PCA by batch or by group.
- Access to relevant directories, files and user guide from GUI
   
### Changed

### Fixed
- Log file can be accessed from changed output folder via symlink

## [1.2.8] - 2024-02-16
 
### Added
   
### Changed
- Minor functional changes, mainly "Black" formatting
### Fixed

## [1.2.7] - 2024-01-19
Stable version

### Changed
 - GUI text to become more clear

### Fixed
- Install in Windows/Mac/Linux working properly in Python >3.9 (3.12.1 still might arise some dependencies problems)
- Proper template for normalization report
- 50 counts filter applies now to housekeeping and endogenous genes
- Ponderated mean calculation from refgenes for content normalization scaling factor

## [1.2.6] - 2024-01-19
Stable version
### Added
 - Changelog.md
### Changed
 - Proper timing


## [1.2.5] - 2024-01-17
 
### Added
   
### Changed
- Improved data handling and visualization
### Fixed
 - Alternative negative controls selection
 - Kruskal-Wallis and Wilcoxon filtering

## [1.2.4] - 2023-12-18
### Added
   
### Changed
- Redesigned data management between steps
- Improved PCA visualization
- Improved recursive back and forth analysis
### Fixed
- Options for pipeline1 normalization and connection to gui

## [1.2.3] - 2023-11-27
Major update
### Added
   - pydeseq2 mean of ratios implementation
   - RUVgnorm normalization

### Changed
- Tab system to include both scaling factors and RUVg normalization
### Fixed
- Path managing
- Group management in output files
- Code readability


The format is based on [Keep a Changelog](http://keepachangelog.com/).