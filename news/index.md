# Changelog

## Changes in version 1.5.3

### New Features

- Exported the previously internal function
  [`overdispFromBoots()`](https://smped.github.io/transmogR/reference/digestSalmon.md)
  for standalone importing of bootstraps

## Changes in version 1.5.2

### Bug Fixes

- Switched from ComplexUpset to SimpleUpset as dependency

## Changes in version 1.3.8

### New Features

- Added
  [`shiftByVar()`](https://smped.github.io/transmogR/reference/shiftByVar.md)
  to produce shifted coordinates which match those after incorporation
  of variants
- Added
  [`cleanVariants()`](https://smped.github.io/transmogR/reference/cleanVariants-methods.md)
  to identify and resolve overlapping variants

## Changes in version 1.3.1

### Improvements

- Changed default behaviour of
  [`digestSalmon()`](https://smped.github.io/transmogR/reference/digestSalmon.md)
  to exclude ‘TPM’ and ‘effectiveLength’ assays, which are now optional
  via the `extra_assays` argument

## Changes in version 1.1.1

### New Features

- Added
  [`digestSalmon()`](https://smped.github.io/transmogR/reference/digestSalmon.md)

## Changes in version 1.0.0

### Major Changes

- Initial Bioconductor release
