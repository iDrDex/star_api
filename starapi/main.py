__author__ = 'dex', 'dhimmel'

import sys
import os
import shutil
import urllib2
import StringIO
import gzip
import re

from funcy import cat, first, re_all
import psycopg2
import psycopg2.extras
import pandas as pd
import numpy as np
import conf

###connect to DB###
import db_conf #PRIVATE
conn = psycopg2.connect(db_conf.DB_PARAMATERS)
cursor = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

def __getMatrixNumHeaderLines(inStream):
    import re

    p = re.compile(r'^"ID_REF"')
    for i, line in enumerate(inStream):
        if p.search(line):
            return i


def matrix_filenames(series_id, platform_id):
    gse_name = query_record(series_id, "series")['gse_name']
    yield "%s/%s_series_matrix.txt.gz" % (gse_name, gse_name)

    gpl_name = query_record(platform_id, "platform")['gpl_name']
    yield "%s/%s-%s_series_matrix.txt.gz" % (gse_name, gse_name, gpl_name)


def get_matrix_filename(series_id, platform_id):
    filenames = list(matrix_filenames(series_id, platform_id))
    mirror_filenames = (os.path.join(conf.SERIES_MATRIX_MIRROR, filename) for filename in filenames)
    mirror_filename = first(filename for filename in mirror_filenames if os.path.isfile(filename))
    if mirror_filename:
        return mirror_filename

    for filename in filenames:
        print 'Loading URL', conf.SERIES_MATRIX_URL + filename, '...'
        try:
            res = urllib2.urlopen(conf.SERIES_MATRIX_URL + filename)
        except urllib2.URLError:
            pass
        else:
            mirror_filename = os.path.join(conf.SERIES_MATRIX_MIRROR, filename)
            print 'Cache to', mirror_filename

            directory = os.path.dirname(mirror_filename)
            if not os.path.exists(directory):
                os.makedirs(directory)
            with open(mirror_filename, 'wb') as f:
                shutil.copyfileobj(res, f)

            return mirror_filename

    raise LookupError("Can't find matrix file for series %s, platform %s"
                      % (series_id, platform_id))

def get_data(series_id, platform_id, impute = False):
    matrixFilename = get_matrix_filename(series_id, platform_id)
    # setup data for specific platform
    for attempt in (0, 1):
        try:
            headerRows = __getMatrixNumHeaderLines(gzip.open(matrixFilename))
            na_values = ["null", "NA", "NaN", "N/A", "na", "n/a", ""]
            data = pd.io.parsers.read_table(gzip.open(matrixFilename),
                                            skiprows=headerRows,
                                            index_col=["ID_REF"],
                                            na_values=na_values,
                                            skipfooter=1,
                                            engine='python')
        except IOError as e:
            # In case we have corrupt file
            print "Failed loading %s: %s" % (matrixFilename, e)
            os.remove(matrixFilename)
            if attempt:
                raise
            matrixFilename = get_matrix_filename(series_id, platform_id)
    data = clean_data(data) #drop samples
    if len(data.columns) == 1:
        data = data.dropna()
    elif impute:
        data = impute_data(data)
    data = log_data(data) #logc

    data.index = data.index.astype(str)
    data.index.name = "probe"
    data.columns.name = 'gsm_name'
    for column in data.columns:
        data[column] = data[column].astype(np.float64)

    # data.to_csv("float64.data.csv")
    return data


def get_platform_probes(platform_id):
    sql = "select * from platform_probe where platform_id = %s"
    return pd.read_sql(sql, conn, "probe", params=(platform_id,))


def query_platform_probes(gpl_name):
    platform_id = query_record(gpl_name, "platform", "gpl_name")['id']
    return get_platform_probes(platform_id)


def get_samples(series_id, platform_id):
    sql = "select * from sample where series_id = %s and platform_id = %s"
    return pd.read_sql(sql, conn, "id", params=(series_id, platform_id,))


def query_samples(gse_name, gpl_name):
    series_id = query_record(gse_name, "series", "gse_name")['id']
    platform_id = query_record(gpl_name, "platform", "gpl_name")['id']
    return get_samples(series_id, platform_id)


def get_gene_data(series_id, platform_id):
    data = get_data(series_id, platform_id)
    platform_probes = get_platform_probes(platform_id)
    gene_data = platform_probes[['mygene_sym', 'mygene_entrez']] \
        .join(data) \
        .set_index(['mygene_sym', 'mygene_entrez'])
    gene_data.columns.name = 'gsm_name'
    return gene_data


def query_record(id, table, id_field="id"):
    sql = """select * from %s where %s """ % (table, id_field) + """= %s"""
    # print sql
    cursor.execute(sql, (id,))
    return cursor.fetchone()


def query_gene_data(gse_name, gpl_name):
    series_id = query_record(gse_name, "series", "gse_name")['id']
    platform_id = query_record(gpl_name, "platform", "gpl_name")['id']
    gene_data = get_gene_data(series_id, platform_id)
    gene_data.columns = gene_data.columns + "_" + gpl_name + "_" + gse_name
    return gene_data

def query_data(gse_name, gpl_name, impute=False):
    series_id = query_record(gse_name, "series", "gse_name")['id']
    platform_id = query_record(gpl_name, "platform", "gpl_name")['id']
    data = get_data(series_id, platform_id, impute)
    # data.columns = data.columns + "_" + gpl_name + "_" + gse_name
    return data

def query_tags_annotations(tokens):
    df = pd.read_sql('''
        SELECT
            sample_id,
            sample.gsm_name,
            annotation,
            series_annotation.series_id,
            series.gse_name,
            series_annotation.platform_id,
            platform.gpl_name,
            tag.tag_name
        FROM
            sample_annotation
            JOIN sample ON (sample_annotation.sample_id = sample.id)
            JOIN series_annotation ON (sample_annotation.serie_annotation_id = series_annotation.id)
            JOIN platform ON (series_annotation.platform_id = platform.id)
            JOIN tag ON (series_annotation.tag_id = tag.id)
            JOIN series ON (series_annotation.series_id = series.id)

        WHERE
            tag.tag_name ~* %(tags)s
    ''', conn, params={'tags': '^(%s)$' % '|'.join(map(re.escape, tokens))})
    # wide = get_wide_annotations(df, tokens)
    return df

def get_unique_annotations(df):
    df= df.query(""" annotation != ''""")\
            .groupby(['sample_id', 'series_id', 'platform_id', 'gsm_name', 'gpl_name'],
                            as_index=False)\
            .filter(lambda x: len(x) == 1)
    return get_wide_annotations(df)

def get_wide_annotations(df):
    tokens = df.tag_name.unique()
    # df = df.groupby(['sample_id', 'series_id', 'platform_id', 'gsm_name', 'gpl_name'],
    #                 as_index=False).filter(lambda x: len(x) == 1) #extracts unique
    # Make tag columns
    df.tag_name = df.tag_name.str.lower()
    df.annotation = df.annotation.str.lower()
    # create outcome column
    # df['outcome'] = None

    for tag in tokens:
        tag_name = tag.lower()
        df[tag_name] = df[df.tag_name == tag_name].annotation

    return df

def get_annotations(case_query, control_query, modifier_query=""):
    # Fetch all relevant data
    queries = [case_query, control_query, modifier_query]
    tokens = set(cat(re_all('[a-zA-Z]\w*', query) for query in queries))
    df = query_tags_annotations(tokens)

    # Make tag columns
    df.tag_name = df.tag_name.str.lower()
    df.annotation = df.annotation.str.lower()
    for tag in tokens:
        tag_name = tag.lower()
        df[tag_name] = df[df.tag_name == tag_name].annotation

    # Select only cells with filled annotations
    df = df.drop(['tag_name', 'annotation'], axis=1)
    df = df.groupby(['sample_id', 'series_id', 'platform_id', 'gsm_name', 'gpl_name'],
                    as_index=False).first()

    df = df.convert_objects(convert_numeric=True)

    # Apply case/control/modifier
    if modifier_query:
        df = df.query(modifier_query.lower())
    case_df = df.query(case_query.lower())
    control_df = df.query(control_query.lower())

    # Set 0 and 1 for analysis
    overlap_df = df.ix[set(case_df.index).intersection(set(control_df.index))]

    df['sample_class'] = None
    df['sample_class'].ix[case_df.index] = 1
    df['sample_class'].ix[control_df.index] = 0
    df['sample_class'].ix[overlap_df.index] = -1

    return df.dropna(subset=["sample_class"])

import numexpr as ne

def log_data(df):
    if is_logged(df):
        return df

    data = df.values
    floor = np.abs(np.nanmin(data, axis=0))
    res = ne.evaluate('log(data + floor + 1) / log(2)')
    return pd.DataFrame(res, index=df.index, columns=df.columns)

def is_logged(df):
    return np.max(df.values) < 10

def impute_data(data):
    import rpy2.robjects as robjects
    r = robjects.r
    import pandas.rpy.common as com
    r.library("impute")
    r_data = com.convert_to_r_matrix(data)
    r_imputedData = r['impute.knn'](r_data)
    npImputedData = np.asarray(r_imputedData[0])
    imputedData = pd.DataFrame(npImputedData)
    imputedData.index = data.index
    imputedData.columns = data.columns
    return imputedData

def drop_missing_genes(data, naLimit=0.5):
    """Filters a data frame to weed out cols with missing data"""
    thresh = len(data.columns) * (1 - naLimit)
    return data.dropna(thresh=thresh, axis="rows")

def drop_missing_samples(data, naLimit=0.8):
    """Filters a data frame to weed out cols with missing data"""
    thresh = len(data.index) * (1 - naLimit)
    return data.dropna(thresh=thresh, axis="columns")

def translate_negative_cols(data):
    """Translate the minimum value of each col to 1"""
    data = data.replace([np.inf, -np.inf], np.nan) #replace infinities
    return data + np.abs(np.min(data)) + 1

def clean_data(data):
    """convenience function to trannslate the data before analysis"""
    if not data.empty:
        # data = log_data(translate_negative_cols(data))
        data = drop_missing_samples(data)
    return data
