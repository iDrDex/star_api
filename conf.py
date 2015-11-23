import os

SERIES_MATRIX_URL = 'ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/'
SERIES_MATRIX_MIRROR = "geo_mirror/DATA/SeriesMatrix/"
os.makedirs(SERIES_MATRIX_MIRROR, exist_ok=True)
CSV_CACHE = "csv"
os.makedirs(CSV_CACHE, exist_ok=True)