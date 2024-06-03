# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v0.2.1]
### Added
- System tests based on single chromosome subsample

### Changed
- Match default parameters with wf-human-variation
- Copy CN from INFO field to sample column in vcf output
- Metadata and blacklist can be parsed from package resources
- Metadata and blacklist duplicated with omitting "chr" from chromosome names

### Fixed
- Missing dependencies
- Parsing mosdepth with number/string chromosome names

## [v0.2.0]
### Changed
- Forked from upstread spectre project.
- Reworked package to be packaged.
