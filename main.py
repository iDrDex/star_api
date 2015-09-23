__author__ = 'dex'
from funcy import str_join, first, log_durations
import conf, urllib2, os, shutil, gzip, psycopg2, psycopg2.extras, pandas as pd, numpy as np
# get a connection, if a connect cannot be made an exception will be raised here
conn = psycopg2.connect(conf.DB_PARAMATERS)
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

def get_data(series_id, platform_id):
    matrixFilename = get_matrix_filename(series_id, platform_id)
    # setup data for specific platform
    for attempt in (0, 1):
        try:
            headerRows = __getMatrixNumHeaderLines(gzip.open(matrixFilename))
            na_values = ["null", "NA", "NaN", "N/A", "na", "n/a"]
            data = pd.io.parsers.read_table(gzip.open(matrixFilename),
                                            skiprows=headerRows,
                                            index_col=["ID_REF"],
                                            na_values=na_values,
                                            skipfooter=1,
                                            engine='python')
        except IOError as e:
            # In case we have cirrupt file
            print "Failed loading %s: %s" % (matrixFilename, e)
            os.remove(matrixFilename)
            if attempt:
                raise
            matrixFilename = get_matrix_filename(series_id, platform_id)

    data.index = data.index.astype(str)
    data.index.name = "probe"
    for column in data.columns:
        data[column] = data[column].astype(np.float64)
    # return data.head(100)
    return data

def get_platform_probes(platform_id):
    sql = "select * from platform_probe where platform_id = %s"
    return pd.read_sql(sql, conn, "probe", params=(platform_id,))

def query_platform_probes(gpl_name):
    platform_id = query_record(gpl_name, "platform", "gpl_name")['id']
    return get_platform_probes(platform_id)

def get_gene_data(series_id, platform_id):
    sample_data = get_data(series_id, platform_id)
    platform_probes = get_platform_probes(platform_id)
    gene_data = platform_probes[['mygene_sym', 'mygene_entrez']]\
                .join(sample_data)\
                .set_index(['mygene_sym', 'mygene_entrez'])
    return gene_data

def query_record(id, table, id_field = "id"):
    sql = """select * from %s where %s """%(table, id_field) + """= %s"""
    cursor.execute(sql, (id,))
    return cursor.fetchone()

def query_gene_data(gse_name, gpl_name):
    series_id = query_record(gse_name, "series", "gse_name")['id']
    platform_id = query_record(gpl_name, "platform", "gpl_name")['id']
    return get_gene_data(series_id, platform_id)

def query_data(gse_name, gpl_name):
    series_id = query_record(gse_name, "series", "gse_name")['id']
    platform_id = query_record(gpl_name, "platform", "gpl_name")['id']
    return get_data(series_id, platform_id)

import rpy2.robjects as robjects

r = robjects.r
import pandas.rpy.common as com

def dropMissingSamples(data, naLimit=0.8):
    """Filters a data frame to weed out cols with missing data"""
    thresh = len(data.index) * (1 - naLimit)
    return data.dropna(thresh=thresh, axis="columns")


def dropMissingGenes(data, naLimit=0.5):
    """Filters a data frame to weed out cols with missing data"""
    thresh = len(data.columns) * (1 - naLimit)
    return data.dropna(thresh=thresh, axis="rows")

def query_median_gene_data(gse_name, gpl_name):
    gene_data = query_gene_data(gse_name, gpl_name)
    gene_data_median = gene_data\
        .reset_index()\
        .groupby(['mygene_sym', 'mygene_entrez'])\
        .median()
    return gene_data_median

def getCombinedMatrix(names):
    """returns an averaged matrix of expression values over all supplid gses"""
    gse_name, gpl_name = names[0]
    m = query_median_gene_data(gse_name, gpl_name)
    for (gse_name, gpl_name) in names[1:]:
        print gse_name, gpl_name,
        median_gene_data = query_median_gene_data(gse_name, gpl_name)
        print median_gene_data.shape
        m = median_gene_data.join(m, how="outer")
    return m

def getCombat(m, samples, impute=False):
    if impute:
        m = getImputed(m)
    else: #drop genes with missing data
        m = m.dropna(axis=1)
    samples = samples.ix[m.columns]
    edata = com.convert_to_r_matrix(m)
    batch = robjects.StrVector(samples['GSE'] + '_' + samples['GPL'])
    pheno = robjects.FactorVector(samples.sample_class)
    r.library("sva")
    fmla = robjects.Formula('~pheno')
    # fmla.environment['pheno'] = r['as.factor'](pheno)
    fmla.environment['pheno'] = pheno
    mod = r['model.matrix'](fmla)
    r_combat_edata = r.ComBat(dat=edata, batch=batch, mod=mod, par_prior=True, prior_plots=False)
    combat = pd.DataFrame(np.asmatrix(r_combat_edata))
    combat.index = m.index
    combat.columns = m.columns
    return combat

def getImputed(data):
    r.library("impute")
    r_data = com.convert_to_r_matrix(data)
    r_imputedData = r['impute.knn'](r_data)
    npImputedData = np.asarray(r_imputedData[0])
    imputedData = pd.DataFrame(npImputedData)
    imputedData.index = data.index
    imputedData.columns = data.columns
    return imputedData

if __name__ == "__main__":
    # print query_data("GSE1", "GPL7")
    # print query_data("GSE3", "GPL9")
    names = ("GSE1", "GPL7"),("GSE3", "GPL9")
    getCombinedMatrix(names)
    # print get_platform_probes(1)[['mygene_sym', 'mygene_entrez']]

