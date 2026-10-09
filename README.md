# Exact Raster Boundary — Version 2.0.1

## Global-ready architecture

Version 2 removes the dependency on Rasterio, GeoPandas and Shapely.

The plugin uses:

- QGIS / PyQGIS
- GDAL / OGR
- QGIS geometry engine
- NumPy for block-wise raster mask construction

These are part of the normal QGIS Python ecosystem, so users do not need to
install a separate GIS Python environment.

## What it creates

For each raster:

`<raster_name>_boundary.shp`

The boundary is derived from valid raster pixels, not from the rectangular
dataset extent.

## Valid-pixel logic

The plugin considers a pixel valid when:

1. GDAL's dataset mask considers it valid;
2. it is not transparent through an alpha band;
3. it is not the raster's explicit NoData value;
4. if enabled, it is not RGB `(0,0,0)` background;
5. for a single-band raster, zero is treated as background.

The black RGB option is intended for rotated imagery with black corners or
black outside-background.

## Memory-safe processing

RGB and mask calculations are processed in blocks (default 1024 x 1024)
rather than loading the entire raster into memory.

## Installation

QGIS:

Plugins -> Manage and Install Plugins -> Install from ZIP

Select:

`ExactRasterBoundary_V2_Global.zip`

## Recommended publishing process

For the official QGIS Plugin Repository, test the plugin on the QGIS versions
you intend to support, then submit the ZIP through the QGIS plugin repository.

## Important

For very large rasters, polygonization can still take time because the valid
mask must be converted into vector geometry. This is normal for an exact
pixel-derived footprint.

## Supported formats

- TIF
- TIFF
- IMG
- ECW
- JP2
- JPG
- JPEG
- PNG


## QGIS version compatibility

This package declares compatibility with QGIS 3.28 through the QGIS 3.x series (metadata maximum: 3.99), including QGIS 3.44.10-Solothurn. It targets QGIS 3.x APIs; QGIS 4.x compatibility is not claimed by this build.
