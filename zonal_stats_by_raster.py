"""Calculate stats per landcover code type."""
import argparse
import datetime
import glob
import os
import logging
import shutil
import sys
import tempfile

from osgeo import gdal
import pygeoprocessing
import numpy

# set a 1GB limit for the cache
gdal.SetCacheMax(2**30)


logging.basicConfig(
    level=logging.DEBUG,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(pathname)s.%(funcName)s:%(lineno)d] %(message)s'),
    stream=sys.stdout)
LOGGER = logging.getLogger(__name__)


def get_unique_values(raster_path):
    """Return a list of non-nodata unique values from `raster_path`."""
    nodata = pygeoprocessing.get_raster_info(raster_path)['nodata'][0]
    unique_set = set()
    for offset_data, array in pygeoprocessing.iterblocks((raster_path, 1)):
        unique_set |= set(numpy.unique(array[~numpy.isclose(array, nodata)]))
    return unique_set


def mask_out_op(mask_data, base_data, mask_code, base_nodata):
    """Return 1 where base data == mask_code, 0 or nodata othewise."""
    result = numpy.empty_like(base_data)
    result[:] = base_nodata
    valid_mask = mask_data == mask_code
    result[valid_mask] = base_data[valid_mask]
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Landcover zonal stats')
    parser.add_argument(
        'landcover_raster_path',
        help='Path to landcover raster.')
    parser.add_argument(
        'raster_pattern_list', metavar='sample_raster_path', nargs='+',
        help='One or more paths or wildcard patterns to calculate stats over.')
    args = parser.parse_args()

    if not os.path.isfile(args.landcover_raster_path):
        raise ValueError('%s is not a valid file' % args.landcover_raster_path)

    sample_raster_path_list = []
    for file_pattern in args.raster_pattern_list:
        sample_raster_path_list.extend(glob.glob(file_pattern))

    if not sample_raster_path_list:
        raise ValueError(
            'No files found matching: %s' % args.raster_pattern_list)

    LOGGER.info('calculating unique landcode values')
    lulc_nodata = pygeoprocessing.get_raster_info(
        args.landcover_raster_path)['nodata']
    unique_values = get_unique_values(args.landcover_raster_path)
    LOGGER.debug('unique landcode values: %s', unique_values)

    working_dir = tempfile.mkdtemp(
        prefix="zonal_stats_temp_workspace_", dir='.')
    time_as_str = ''.join(
        [x if x not in [' ', ':', '.', '-'] else '_'
         for x in str(datetime.datetime.now())])
    stats_table = open('stats_table_%s.csv' % time_as_str, 'w')
    stats_table.write('raster_filename,lucode,min,max,mean,stdev\n')

    for sample_raster_path in sample_raster_path_list:
        LOGGER.info('processing %s', sample_raster_path)
        base_raster_path_list = [
            args.landcover_raster_path, sample_raster_path]
        aligned_raster_path_list = [
            os.path.join(working_dir, os.path.basename(path))
            for path in base_raster_path_list]
        other_raster_info = pygeoprocessing.get_raster_info(sample_raster_path)
        pygeoprocessing.align_and_resize_raster_stack(
            base_raster_path_list, aligned_raster_path_list, ['mode', 'near'],
            other_raster_info['pixel_size'], 'intersection',
            target_sr_wkt=other_raster_info['projection'])
        for index, mask_code in enumerate(sorted(unique_values)):
            LOGGER.info(
                'processing %d of %d (lucode %d) for %s', index+1,
                len(unique_values), mask_code,
                os.path.basename(sample_raster_path))
            mask_raster_path = os.path.join(working_dir, '%d.tif' % mask_code)
            pygeoprocessing.raster_calculator(
                [(aligned_raster_path_list[0], 1),
                 (aligned_raster_path_list[1], 1),
                 (mask_code, 'raw'),
                 (other_raster_info['nodata'][0], 'raw')],
                mask_out_op, mask_raster_path, gdal.GDT_Float32,
                other_raster_info['nodata'][0])
            raster = gdal.OpenEx(mask_raster_path, gdal.OF_RASTER)
            band = raster.GetRasterBand(1)
            (raster_min, raster_max, raster_mean, raster_stdev) = (
                band.GetStatistics(0, 1))
            band = None
            raster = None
            stats_table.write(
                '%s,%d,%f,%f,%f,%f\n' % (
                    os.path.basename(sample_raster_path), mask_code,
                    raster_min, raster_max, raster_mean, raster_stdev))
    stats_table.close()
    try:
        shutil.rmtree(working_dir)
    except OSError:
        LOGGER.warning('unable to remove temp working dir: %s' % working_dir)
