# -*- coding: utf-8 -*-
"""
Exact Raster Boundary - Version 2.0.1

Global QGIS plugin implementation using:
- QGIS / PyQGIS
- GDAL (osgeo.gdal / osgeo.ogr)
- NumPy (used only for block-wise raster mask construction)

No Rasterio, GeoPandas or Shapely are required.

The algorithm:
1. Opens the raster through GDAL.
2. Builds a Byte validity mask in blocks.
3. Uses the GDAL dataset mask / alpha / NoData information.
4. Optionally removes RGB (0,0,0) pixels.
5. Polygonizes the valid mask using GDAL.
6. Dissolves the polygonized valid pixels using QGIS geometry.
7. Writes one exact boundary Shapefile per input raster.

The output follows the actual valid pixels, including rotated raster
geotransforms, instead of simply using the rectangular raster extent.
"""

import os
import traceback
import tempfile

from qgis.PyQt.QtCore import Qt
from qgis.PyQt.QtWidgets import (
    QAction, QDialog, QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QLineEdit, QPushButton, QCheckBox, QProgressBar,
    QPlainTextEdit, QMessageBox, QFileDialog, QGroupBox,
    QApplication
)

from qgis.core import (
    QgsProject,
    QgsVectorLayer,
    QgsFeature,
    QgsGeometry,
    QgsField,
    QgsFields,
    QgsVectorFileWriter,
    QgsWkbTypes,
    QgsCoordinateReferenceSystem
)

from qgis.PyQt.QtCore import QVariant

try:
    from osgeo import gdal, ogr, osr
except Exception as e:
    gdal = None
    ogr = None
    osr = None
    GDAL_IMPORT_ERROR = str(e)

try:
    import numpy as np
    NUMPY_IMPORT_ERROR = None
except Exception as e:
    np = None
    NUMPY_IMPORT_ERROR = str(e)


SUPPORTED_EXTENSIONS = {
    ".tif", ".tiff", ".img", ".ecw", ".jp2", ".jpg", ".jpeg", ".png"
}

GDAL_USE_EXCEPTIONS = getattr(gdal, "UseExceptions", None)
if GDAL_USE_EXCEPTIONS:
    gdal.UseExceptions()


def is_alpha_band(band):
    """Return True when a GDAL band is an alpha band."""
    try:
        color_interp = band.GetColorInterpretation()
        return color_interp == gdal.GCI_AlphaBand
    except Exception:
        return False


def nodata_equal(arr, nodata):
    """Return a NumPy boolean array identifying NoData pixels."""
    if nodata is None:
        return np.zeros(arr.shape, dtype=bool)

    if np.issubdtype(arr.dtype, np.floating) and np.isnan(nodata):
        return np.isnan(arr)

    return arr == nodata


def build_valid_mask_raster(
    src_ds,
    mask_path,
    black_rgb_invalid=True,
    block_size=1024,
    progress_callback=None
):
    """
    Build a Byte GDAL raster:
      1 = valid pixel
      0 = invalid/background

    The mask is created block-by-block to avoid loading a huge raster into
    memory at once.
    """

    if np is None:
        raise RuntimeError(
            "NumPy is required by QGIS's Python environment for block-wise "
            "RGB mask creation."
        )

    width = src_ds.RasterXSize
    height = src_ds.RasterYSize
    bands = src_ds.RasterCount

    driver = gdal.GetDriverByName("GTiff")
    if driver is None:
        raise RuntimeError("GDAL GTiff driver is unavailable.")

    if os.path.exists(mask_path):
        driver.Delete(mask_path)

    mask_ds = driver.Create(
        mask_path,
        width,
        height,
        1,
        gdal.GDT_Byte,
        options=[
            "TILED=YES",
            "COMPRESS=LZW",
            "BIGTIFF=IF_SAFER"
        ]
    )

    if mask_ds is None:
        raise RuntimeError("Could not create temporary mask raster.")

    mask_ds.SetGeoTransform(src_ds.GetGeoTransform())
    projection = src_ds.GetProjection()
    if projection:
        mask_ds.SetProjection(projection)

    out_band = mask_ds.GetRasterBand(1)
    out_band.SetNoDataValue(0)

    # Dataset mask is the preferred GDAL validity mask.
    dataset_mask_band = None
    try:
        if bands > 0:
            dataset_mask_band = src_ds.GetRasterBand(1).GetMaskBand()
    except Exception:
        dataset_mask_band = None

    # Identify alpha band, if one exists.
    alpha_band = None
    for bidx in range(1, bands + 1):
        band = src_ds.GetRasterBand(bidx)
        if is_alpha_band(band):
            alpha_band = band
            break

    # For RGB background detection, use the first three bands.
    rgb_available = bands >= 3
    r_band = src_ds.GetRasterBand(1) if rgb_available else None
    g_band = src_ds.GetRasterBand(2) if rgb_available else None
    b_band = src_ds.GetRasterBand(3) if rgb_available else None

    # NoData values are checked band-by-band. Usually the first band is
    # enough, but checking all bands is safer for multiband imagery.
    nodata_values = []
    for bidx in range(1, bands + 1):
        band = src_ds.GetRasterBand(bidx)
        try:
            nodata_values.append(band.GetNoDataValue())
        except Exception:
            nodata_values.append(None)

    total_blocks = ((height + block_size - 1) // block_size) * (
        (width + block_size - 1) // block_size
    )
    block_number = 0

    for yoff in range(0, height, block_size):
        ysize = min(block_size, height - yoff)

        for xoff in range(0, width, block_size):
            xsize = min(block_size, width - xoff)

            # Start with GDAL's dataset mask where available.
            if dataset_mask_band is not None:
                base_mask = dataset_mask_band.ReadAsArray(
                    xoff, yoff, xsize, ysize
                )
                if base_mask is None:
                    valid = np.ones((ysize, xsize), dtype=bool)
                else:
                    valid = base_mask > 0
            else:
                valid = np.ones((ysize, xsize), dtype=bool)

            # Apply alpha transparency.
            if alpha_band is not None:
                alpha = alpha_band.ReadAsArray(
                    xoff, yoff, xsize, ysize
                )
                if alpha is not None:
                    valid &= alpha > 0

            # Apply explicit NoData values.
            for bidx, nodata in enumerate(nodata_values, start=1):
                if nodata is None:
                    continue

                # Do not re-read alpha as a data band for this purpose.
                if alpha_band is not None and bidx == alpha_band.GetBand():
                    continue

                band = src_ds.GetRasterBand(bidx)
                arr = band.ReadAsArray(xoff, yoff, xsize, ysize)
                if arr is None:
                    continue

                invalid_nd = nodata_equal(arr, nodata)
                valid &= ~invalid_nd

                # NaN is always invalid for floating point data.
                if np.issubdtype(arr.dtype, np.floating):
                    valid &= ~np.isnan(arr)

            # Optional black RGB background removal.
            if black_rgb_invalid and rgb_available:
                r = r_band.ReadAsArray(xoff, yoff, xsize, ysize)
                g = g_band.ReadAsArray(xoff, yoff, xsize, ysize)
                b = b_band.ReadAsArray(xoff, yoff, xsize, ysize)

                if r is not None and g is not None and b is not None:
                    black = (r == 0) & (g == 0) & (b == 0)
                    valid &= ~black

            # Single-band zero background is treated as invalid.
            elif bands == 1:
                arr = src_ds.GetRasterBand(1).ReadAsArray(
                    xoff, yoff, xsize, ysize
                )
                if arr is not None:
                    valid &= arr != 0
                    if np.issubdtype(arr.dtype, np.floating):
                        valid &= ~np.isnan(arr)

            out = np.where(valid, 1, 0).astype(np.uint8)
            out_band.WriteArray(out, xoff, yoff)

            block_number += 1
            if progress_callback:
                progress_callback(
                    block_number / float(max(total_blocks, 1))
                )

    out_band.FlushCache()
    mask_ds.FlushCache()

    return mask_ds


def polygonize_mask(mask_path, polygon_path):
    """Polygonize the Byte mask using GDAL."""

    src = gdal.Open(mask_path, gdal.GA_ReadOnly)
    if src is None:
        raise RuntimeError("Could not reopen temporary mask raster.")

    band = src.GetRasterBand(1)
    if band is None:
        src = None
        raise RuntimeError("Temporary mask has no raster band.")

    shp_driver = ogr.GetDriverByName("ESRI Shapefile")
    if shp_driver is None:
        src = None
        raise RuntimeError("GDAL ESRI Shapefile driver is unavailable.")

    if os.path.exists(polygon_path):
        shp_driver.DeleteDataSource(polygon_path)

    out_ds = shp_driver.CreateDataSource(polygon_path)
    if out_ds is None:
        src = None
        raise RuntimeError("Could not create polygon output.")

    projection = src.GetProjection()
    srs = None
    if projection:
        srs = osr.SpatialReference()
        srs.ImportFromWkt(projection)

    layer = out_ds.CreateLayer(
        "valid_pixels",
        srs=srs,
        geom_type=ogr.wkbPolygon
    )

    field = ogr.FieldDefn("DN", ogr.OFTInteger)
    layer.CreateField(field)

    # GDAL polygonize creates polygons for connected equal-value regions.
    result = gdal.Polygonize(
        band,
        None,
        layer,
        0,
        [],
        callback=None
    )

    layer.SyncToDisk()
    out_ds.FlushCache()
    out_ds = None
    src = None

    if result != 0:
        raise RuntimeError(
            f"GDAL Polygonize returned error code {result}."
        )


def union_valid_polygons(polygon_path, output_shp, source_name):
    """
    Read polygonized valid regions with QGIS, keep DN=1, union all valid
    geometries, and save one exact boundary feature as an ESRI Shapefile.
    """

    layer = QgsVectorLayer(polygon_path, "valid_pixels", "ogr")
    if not layer.isValid():
        raise RuntimeError("Polygonized vector layer is invalid.")

    dn_index = layer.fields().indexOf("DN")

    geometries = []
    valid_feature_count = 0

    for feature in layer.getFeatures():
        if dn_index >= 0:
            dn_value = feature["DN"]
            if dn_value is None or int(dn_value) != 1:
                continue

        geom = feature.geometry()
        if geom is not None and not geom.isEmpty():
            geometries.append(QgsGeometry(geom))
            valid_feature_count += 1

    if not geometries:
        raise RuntimeError("No valid-pixel polygons were generated.")

    # QGIS performs the geometry union; no Shapely is needed.
    union_geom = QgsGeometry.unaryUnion(geometries)

    if union_geom is None or union_geom.isEmpty():
        raise RuntimeError("The union of valid pixels is empty.")

    # Repair minor geometry issues if possible.
    if not union_geom.isGeosValid():
        fixed = union_geom.makeValid()
        if fixed is not None and not fixed.isEmpty():
            union_geom = fixed

    crs = layer.crs()

    fields = QgsFields()
    fields.append(QgsField("RASTER", QVariant.String, "String", 254))
    fields.append(QgsField("VALID_POLY", QVariant.Int))
    fields.append(QgsField("AREA", QVariant.Double, "Double", 20, 3))

    geometry_type = QgsWkbTypes.multiType(QgsWkbTypes.Polygon)

    mem_uri = (
        f"MultiPolygon?crs={crs.authid()}"
        if crs.isValid()
        else "MultiPolygon"
    )

    mem = QgsVectorLayer(mem_uri, "exact_boundary", "memory")
    if not mem.isValid():
        raise RuntimeError("Could not create temporary QGIS vector layer.")

    provider = mem.dataProvider()
    provider.addAttributes(fields.toList())
    mem.updateFields()

    feature = QgsFeature(mem.fields())
    feature.setGeometry(union_geom)
    feature["RASTER"] = source_name
    feature["VALID_POLY"] = valid_feature_count
    feature["AREA"] = union_geom.area()
    provider.addFeature(feature)
    mem.updateExtents()

    # Ensure destination directory exists.
    os.makedirs(os.path.dirname(output_shp), exist_ok=True)

    # Remove old shapefile sidecars if present.
    for ext in (".shp", ".shx", ".dbf", ".prj", ".cpg", ".qpj"):
        candidate = os.path.splitext(output_shp)[0] + ext
        if os.path.exists(candidate):
            try:
                os.remove(candidate)
            except OSError as exc:
                raise RuntimeError(
                    "Could not remove existing output sidecar '{}': {}".format(
                        candidate, exc
                    )
                ) from exc

    transform_context = QgsProject.instance().transformContext()

    options = QgsVectorFileWriter.SaveVectorOptions()
    options.driverName = "ESRI Shapefile"
    options.fileEncoding = "UTF-8"
    options.layerName = os.path.splitext(
        os.path.basename(output_shp)
    )[0][:30]

    result = QgsVectorFileWriter.writeAsVectorFormatV3(
        mem,
        output_shp,
        transform_context,
        options
    )

    # V3 returns (WriterError, error_message, new_filename, new_layer)
    if isinstance(result, tuple):
        error_code = result[0]
        error_message = result[1] if len(result) > 1 else ""
    else:
        error_code = result
        error_message = ""

    if error_code != QgsVectorFileWriter.NoError:
        raise RuntimeError(
            f"Could not write Shapefile: {error_message}"
        )

    return {
        "area": union_geom.area(),
        "valid_polygons": valid_feature_count,
        "crs": crs.authid() if crs.isValid() else "Unknown"
    }


def process_one_raster(
    input_raster,
    output_shp,
    black_rgb_invalid=True,
    block_size=1024,
    progress_callback=None
):
    """Complete exact-boundary workflow for one raster."""

    src = gdal.Open(input_raster, gdal.GA_ReadOnly)
    if src is None:
        raise RuntimeError("GDAL could not open the raster.")

    if not src.GetProjection():
        src = None
        raise RuntimeError(
            "Raster has no CRS/projection. Assign a CRS before processing."
        )

    temp_dir = tempfile.mkdtemp(prefix="exact_boundary_")

    try:
        mask_path = os.path.join(temp_dir, "valid_mask.tif")
        polygon_path = os.path.join(temp_dir, "valid_pixels.shp")

        def mask_progress(p):
            if progress_callback:
                progress_callback(p * 0.70)

        build_valid_mask_raster(
            src,
            mask_path,
            black_rgb_invalid=black_rgb_invalid,
            block_size=block_size,
            progress_callback=mask_progress
        )

        src = None

        if progress_callback:
            progress_callback(0.70)

        polygonize_mask(mask_path, polygon_path)

        if progress_callback:
            progress_callback(0.85)

        result = union_valid_polygons(
            polygon_path,
            output_shp,
            os.path.basename(input_raster)
        )

        if progress_callback:
            progress_callback(1.0)

        return result

    finally:
        src = None

        # Ignore cleanup errors because shutil_rmtree uses ignore_errors=True.
        shutil_rmtree(temp_dir)


def shutil_rmtree(path):
    import shutil
    shutil.rmtree(path, ignore_errors=True)


class ExactRasterBoundaryDialog(QDialog):

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Exact Raster Boundary — QGIS 3.x")
        self.resize(820, 680)
        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout(self)

        title = QLabel("<h2>Exact Raster Boundary</h2>")
        subtitle = QLabel(
            "Create the actual valid-data footprint of raster imagery. "
            "The rectangular raster extent is not used as the boundary."
        )
        subtitle.setWordWrap(True)
        layout.addWidget(title)
        layout.addWidget(subtitle)

        io_box = QGroupBox("Input / Output")
        grid = QGridLayout(io_box)

        grid.addWidget(QLabel("Input raster folder:"), 0, 0)
        self.input_edit = QLineEdit()
        grid.addWidget(self.input_edit, 0, 1)

        in_btn = QPushButton("Browse...")
        in_btn.clicked.connect(self.browse_input)
        grid.addWidget(in_btn, 0, 2)

        grid.addWidget(QLabel("Output boundary folder:"), 1, 0)
        self.output_edit = QLineEdit()
        grid.addWidget(self.output_edit, 1, 1)

        out_btn = QPushButton("Browse...")
        out_btn.clicked.connect(self.browse_output)
        grid.addWidget(out_btn, 1, 2)

        layout.addWidget(io_box)

        options_box = QGroupBox("Boundary Options")
        options_layout = QVBoxLayout(options_box)

        self.black_check = QCheckBox(
            "Treat black RGB pixels (0, 0, 0) as invalid background"
        )
        self.black_check.setChecked(True)
        options_layout.addWidget(self.black_check)

        self.subfolder_check = QCheckBox(
            "Process raster files in subfolders"
        )
        self.subfolder_check.setChecked(False)
        options_layout.addWidget(self.subfolder_check)

        self.load_check = QCheckBox(
            "Automatically load generated boundaries into QGIS"
        )
        self.load_check.setChecked(True)
        options_layout.addWidget(self.load_check)

        layout.addWidget(options_box)

        info = QLabel(
            "<b>Supported:</b> TIF/TIFF, IMG, ECW, JP2, JPG/JPEG, PNG"
            "<br><b>Engine:</b> QGIS + GDAL + QGIS geometry engine"
            "<br><b>External GIS packages:</b> None required"
        )
        info.setWordWrap(True)
        layout.addWidget(info)

        self.status_label = QLabel("Ready.")
        layout.addWidget(self.status_label)

        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        self.progress.setValue(0)
        layout.addWidget(self.progress)

        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        layout.addWidget(self.log, 1)

        buttons = QHBoxLayout()

        self.run_button = QPushButton("Create Exact Boundaries")
        self.run_button.setDefault(True)
        self.run_button.clicked.connect(self.run_processing)
        buttons.addWidget(self.run_button)

        clear_btn = QPushButton("Clear Log")
        clear_btn.clicked.connect(self.log.clear)
        buttons.addWidget(clear_btn)

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        buttons.addWidget(close_btn)

        layout.addLayout(buttons)

    def browse_input(self):
        folder = QFileDialog.getExistingDirectory(
            self,
            "Select Input Raster Folder"
        )
        if folder:
            self.input_edit.setText(folder)

            if not self.output_edit.text().strip():
                self.output_edit.setText(
                    os.path.join(folder, "exact_raster_boundary")
                )

    def browse_output(self):
        folder = QFileDialog.getExistingDirectory(
            self,
            "Select Output Boundary Folder"
        )
        if folder:
            self.output_edit.setText(folder)

    def write_log(self, message):
        self.log.appendPlainText(message)
        bar = self.log.verticalScrollBar()
        bar.setValue(bar.maximum())

    def collect_rasters(self, root, recursive=False):
        rasters = []

        if recursive:
            for current_root, dirs, files in os.walk(root):
                for name in files:
                    ext = os.path.splitext(name)[1].lower()
                    if ext in SUPPORTED_EXTENSIONS:
                        rasters.append(os.path.join(current_root, name))
        else:
            for name in os.listdir(root):
                path = os.path.join(root, name)
                if os.path.isfile(path):
                    ext = os.path.splitext(name)[1].lower()
                    if ext in SUPPORTED_EXTENSIONS:
                        rasters.append(path)

        return sorted(rasters, key=lambda p: p.lower())

    def run_processing(self):
        if gdal is None or ogr is None:
            QMessageBox.critical(
                self,
                "GDAL unavailable",
                "The QGIS Python environment could not import GDAL/OGR.\n\n"
                + str(GDAL_IMPORT_ERROR)
            )
            return

        if np is None:
            QMessageBox.critical(
                self,
                "NumPy unavailable",
                "The QGIS Python environment could not import NumPy.\n\n"
                + str(NUMPY_IMPORT_ERROR)
            )
            return

        input_folder = self.input_edit.text().strip()
        output_folder = self.output_edit.text().strip()

        if not input_folder or not os.path.isdir(input_folder):
            QMessageBox.warning(
                self,
                "Input folder",
                "Please select a valid input raster folder."
            )
            return

        if not output_folder:
            QMessageBox.warning(
                self,
                "Output folder",
                "Please select an output folder."
            )
            return

        os.makedirs(output_folder, exist_ok=True)

        rasters = self.collect_rasters(
            input_folder,
            self.subfolder_check.isChecked()
        )

        if not rasters:
            QMessageBox.information(
                self,
                "No rasters found",
                "No supported raster files were found in the selected folder."
            )
            return

        self.run_button.setEnabled(False)
        self.log.clear()
        self.progress.setValue(0)

        success = 0
        failed = 0

        self.write_log("=" * 78)
        self.write_log("EXACT RASTER BOUNDARY — VERSION 2")
        self.write_log("=" * 78)
        self.write_log("Engine: QGIS + GDAL + QGIS geometry")
        self.write_log(f"Input : {input_folder}")
        self.write_log(f"Output: {output_folder}")
        self.write_log(f"Files : {len(rasters)}")
        self.write_log("")

        for i, raster in enumerate(rasters, start=1):
            base = os.path.splitext(os.path.basename(raster))[0]
            output_shp = os.path.join(
                output_folder,
                base + "_boundary.shp"
            )

            self.status_label.setText(
                f"Processing {i}/{len(rasters)}: {os.path.basename(raster)}"
            )

            self.write_log(
                f"[{i}/{len(rasters)}] {os.path.basename(raster)}"
            )

            try:
                def one_progress(p):
                    overall = ((i - 1) + p) / float(len(rasters))
                    self.progress.setValue(int(overall * 100))
                    QApplication.processEvents()

                result = process_one_raster(
                    raster,
                    output_shp,
                    black_rgb_invalid=self.black_check.isChecked(),
                    block_size=1024,
                    progress_callback=one_progress
                )

                success += 1

                self.write_log(
                    f"  OK -> {output_shp}\n"
                    f"     Valid regions: {result['valid_polygons']}\n"
                    f"     CRS: {result['crs']}\n"
                    f"     Area: {result['area']:.3f}"
                )

                if self.load_check.isChecked():
                    layer = QgsVectorLayer(
                        output_shp,
                        base + "_boundary",
                        "ogr"
                    )

                    if layer.isValid():
                        QgsProject.instance().addMapLayer(layer)
                    else:
                        self.write_log(
                            "     Warning: output created but could "
                            "not be loaded into QGIS."
                        )

            except Exception as e:
                failed += 1
                self.write_log(
                    f"  FAILED: {e}\n"
                    f"{traceback.format_exc()}"
                )

            self.progress.setValue(
                int((i / float(len(rasters))) * 100)
            )
            QApplication.processEvents()

        self.write_log("")
        self.write_log("=" * 78)
        self.write_log("PROCESS COMPLETED")
        self.write_log("=" * 78)
        self.write_log(f"Total   : {len(rasters)}")
        self.write_log(f"Success : {success}")
        self.write_log(f"Failed  : {failed}")

        self.status_label.setText(
            f"Completed: {success} successful, {failed} failed."
        )
        self.run_button.setEnabled(True)

        if failed == 0:
            QMessageBox.information(
                self,
                "Completed",
                f"Successfully created {success} exact raster boundary layer(s)."
            )
        else:
            QMessageBox.warning(
                self,
                "Completed with errors",
                f"Created {success} boundary layer(s); {failed} failed.\n"
                "Check the processing log for details."
            )


class ExactRasterBoundary:

    def __init__(self, iface):
        self.iface = iface
        self.action = None
        self.menu = "Exact Raster Boundary"

    def initGui(self):
        self.action = QAction(
            "Exact Raster Boundary",
            self.iface.mainWindow()
        )
        self.action.setToolTip(
            "Create exact valid-data boundaries from raster files"
        )
        self.action.triggered.connect(self.run)

        self.iface.addPluginToMenu(self.menu, self.action)
        self.iface.addToolBarIcon(self.action)

    def unload(self):
        if self.action:
            self.iface.removePluginMenu(
                self.menu,
                self.action
            )
            self.iface.removeToolBarIcon(self.action)

    def run(self):
        dialog = ExactRasterBoundaryDialog(
            self.iface.mainWindow()
        )
        dialog.exec()
