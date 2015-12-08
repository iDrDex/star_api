import os

SERIES_MATRIX_URL = 'ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/'

def configure(directory):
    """Configure starapi to create cache locations in the specified directory"""

    global SERIES_MATRIX_MIRROR
    global CSV_CACHE

    SERIES_MATRIX_MIRROR = os.path.join(directory, 'geo_mirror/DATA/SeriesMatrix/')
    CSV_CACHE = os.path.join(directory, 'csv/')

    if not os.path.exists(SERIES_MATRIX_MIRROR):
        print('creating SERIES_MATRIX_MIRROR directory')
        os.makedirs(SERIES_MATRIX_MIRROR)

    if not os.path.exists(CSV_CACHE):
        print('creating CSV_CACHE directory')
        os.makedirs(CSV_CACHE)
