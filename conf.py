import os

SERIES_MATRIX_URL = 'ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/'

SERIES_MATRIX_MIRROR = "geo_mirror/DATA/SeriesMatrix/"
if not os.path.exists(SERIES_MATRIX_MIRROR): os.makedirs(SERIES_MATRIX_MIRROR)

CSV_CACHE = "csv"
if not os.path.exists(CSV_CACHE): os.makedirs(CSV_CACHE)