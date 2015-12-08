import logging

from main import *
from easydict import EasyDict
from funcy import first, log_durations, imap, memoize, cat, re_all
import numpy as np
import pandas as pd
import scipy.stats as stats

#### create logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

def sanitize(filename):
    return "".join([c for c in filename if c.isalpha() or c.isdigit() or c==' ']).rstrip()

@log_durations(logger.debug)
def get_balanced_permutations(balanced, permutations):
    balanced_permutations = balanced
    if not permutations.empty:
        logger.info('Estimating significance by permutation for %s' % analysis.analysis_name)
        with log_durations(logger.debug, 'Estimating significance by permutation for %s' % analysis.analysis_name):
            recs = []
            permutations.sort_index(inplace=True)#to speed lookups
            # permutations.to_csv("permutations.csv")
            for mygene in balanced.index:
                perms = permutations.ix[mygene]
                random_TE = balanced.ix[mygene].random_TE
                random_side = 'right' if random_TE > 0 else "left"
                random_perms = perms.random_TE.order()
                random_nperms = random_perms.count()
                random_rank = random_perms.searchsorted(random_TE, side=random_side)[0]
                if random_side == "right":
                    random_rank  = random_nperms - random_rank
                random_pval_perm = float(random_rank)/random_nperms
                fixed_TE = balanced.ix[mygene].fixed_TE
                fixed_side = 'right' if fixed_TE > 0 else "left"
                fixed_perms = perms.fixed_TE.order()
                fixed_nperms = fixed_perms.count()
                fixed_rank = fixed_perms.searchsorted(fixed_TE, side=fixed_side)[0]
                if fixed_side == "right":
                    fixed_rank  = fixed_nperms - fixed_rank
                fixed_pval_perm = float(fixed_rank)/fixed_nperms
                rec = dict(random_rank = random_rank,
                           random_nperms = random_nperms,
                           random_pval_perm = random_pval_perm,
                           fixed_rank = fixed_rank,
                           fixed_nperms = fixed_nperms,
                           fixed_pval_perm = fixed_pval_perm)
                recs.append(rec)
            df = pd.DataFrame(recs)
            df.index = balanced.index
            balanced_permutations = balanced.join(df)
    return  balanced_permutations


@log_durations(logger.debug)
def perform_analysis(analysis, debug=False, impute = False, nperm = 0, mygene_filter = None):
    """
    Returns a tuple of sample_df, fold_change, balanced_permutations, permutations
    """
    logger.info('Started %s analysis', analysis.analysis_name)
    # from multiprocessing import Pool
    # pool = Pool(processes=4)

    with log_durations(logger.debug, 'Loading dataframe for %s' % analysis.analysis_name):
        df = get_analysis_df(analysis.case_query, analysis.control_query, analysis.modifier_query)
    debug and df.to_csv("%s.analysis_df.csv" % analysis.analysis_name)

    logger.info('Matching sources: %d' % df.groupby(['series_id', 'platform_id']).ngroups)

    # Remove single-class sources
    query = df.groupby(['series_id', 'platform_id']).sample_class.agg(lambda x: set(x)).map(lambda x: x >= {0, 1})
    df = filter_sources(df, query, 'as single-class')

    # Check for minimum number of samples
    if analysis.min_samples:
        counts = df.groupby(['series_id', 'platform_id']).sample_class.value_counts().unstack()
        query = (counts[0] >= analysis.min_samples) & (counts[1] >= analysis.min_samples)
        df = filter_sources(df, query, 'by min samples')

    # Check number of sources
    sources = df.groupby(['series_id', 'platform_id']).ngroups
    if sources <= 1:
        logger.error("FAIL Can't perform meta-analysis on %s"
                     % ('single source' if sources else 'no data'))
        return df, None, None, None

    # Calculating stats
    analysis.series_count = len(df.series_id.unique())
    analysis.platform_count = len(df.platform_id.unique())
    analysis.sample_count = len(df.sample_id.unique())
    analysis.series_ids = df.series_id.unique().tolist()
    analysis.platform_ids = df.platform_id.unique().tolist()
    analysis.sample_ids = df.sample_id.unique().tolist()
    # analysis.save(update_fields=['series_count', 'platform_count', 'sample_count',
    #                              'series_ids', 'platform_ids', 'sample_ids'])
    logger.info('Stats: %d sources, %d series, %d platforms, %d samples'
                % (sources, analysis.series_count, analysis.platform_count, analysis.sample_count))

    # Load GSE data, make and concat all fold change analyses results.
    # NOTE: we are doing load_gse() lazily here to avoid loading all matrices at once.
    logger.info('Loading data and calculating fold change for %s', analysis.analysis_name)
    with log_durations(logger.debug, 'Load/fold for %s' % analysis.analysis_name):
        gses = (load_gse(df, series_id, impute) for series_id in sorted(df.series_id.unique()))
        debugs = [debug]*df.series_id.nunique()
        nperms = [nperm]*df.series_id.nunique()
        mygene_filters = [mygene_filter]*df.series_id.nunique()

        # start a pool with 4 processes
        fold_change = pd.concat(imap(get_gene_fold_change, gses, debugs, nperms, mygene_filters))
        # fold_change = pd.concat(pool.imap(multi_run_wrapper, zip(gses, debugs, nperms)))
        debug and fold_change.to_csv("%s.fc.csv" % debug)

    #Start metaanalysis
    logger.info('Meta-Analyzing %s', analysis.analysis_name)
    with log_durations(logger.debug, 'Meta analysis for %s' % analysis.analysis_name):
        # logger.info('Meta analysis of real data for %s' % analysis.analysis_name)
        with log_durations(logger.debug, 'meta analysis of real data for %s' % analysis.analysis_name):
            balanced = get_full_meta(fold_change.query("""perm == 0"""), debug=debug)
            debug and balanced.to_csv("%s.meta.csv" % debug)
        # logger.info('Meta-Analyzing of permutations for %s', analysis.analysis_name)
        with log_durations(logger.debug, 'meta analysis of permutations for %s' % analysis.analysis_name):
            permutations = pd.DataFrame()
            fold_change = fold_change.reset_index().sort('perm').set_index('perm')
            for i in range(nperm):
                perm = i + 1
                # logger.info('Meta analysis of permutation %s for %s' % (perm, analysis.analysis_name))
                with log_durations(logger.debug, 'meta analysis of permutation %s / %s for %s' % (perm, nperm, analysis.analysis_name)):
                    # balanced_perm = get_full_meta(fold_change.query("""perm == %s"""%perm), debug=debug)
                    balanced_perm = get_full_meta(fold_change.ix[perm], debug=debug)
                    permutation = balanced_perm[['random_TE', 'fixed_TE']]
                    permutation['perm'] = perm
                    permutations = pd.concat([permutations, permutation])
        balanced_permutations = get_balanced_permutations(balanced, permutations)

    logger.info('DONE %s analysis', analysis.analysis_name)
    return df, fold_change, balanced_permutations, permutations

def filter_sources(df, query, reason):
    start_sources = df.groupby(['series_id', 'platform_id']).ngroups
    new_df = df.set_index(['series_id', 'platform_id']).loc[query].reset_index()
    sources = new_df.groupby(['series_id', 'platform_id']).ngroups
    excluded = start_sources - sources
    if excluded:
        logger.info('Excluded %d source%s %s'% (excluded, 's' if excluded > 1 else '', reason))
    return new_df

@log_durations(logger.debug)
def get_gene_fold_change(gse, debug=False, nperm=0, mygene_filter = None):
    samples = gse.samples

    if 'subset' not in samples.columns:
        samples['subset'] = "NA"

    groups = samples.ix[samples.sample_class >= 0] \
        .groupby(['subset', 'gpl_name'])

    results_list = []

    for group, df in groups:
        subset, gpl = group

        # NOTE: if data has changed then sample ids could be different
        if not set(df["gsm_name"]) <= set(gse.gpl2data[gpl].columns):
            logger.info("skipping %s: sample ids mismatch" % gpl)
            continue

        df = df.set_index("gsm_name")
        data = gse.gpl2data[gpl][df.index]
        probes = gse.gpl2probes[gpl]

        if mygene_filter is None: #filter probes for all genes
            mygene_filter = probes.dropna(subset=['mygene_sym', 'mygene_entrez'])\
                .set_index(['mygene_sym', 'mygene_entrez'])\
                .index.unique()

        mygene_probes = probes.reset_index()\
            .sort(['mygene_sym', 'mygene_entrez'])\
            .set_index(['mygene_sym', 'mygene_entrez'])

        mygene_filter = mygene_probes.index.intersection(mygene_filter).unique()

        data = data.ix[mygene_probes.ix[mygene_filter].probe]
        sample_class = df.ix[data.columns].sample_class

        debug = debug and debug + ".%s_%s_%s" % (gse.name, gpl, subset)
        for perm in range(nperm + 1):
            if perm:
                sample_class = np.random.permutation(sample_class)

            fold_change = get_fold_change(data, sample_class,
                                          debug=debug)
            fold_change = fold_change \
                    .join(probes[['mygene_entrez', 'mygene_sym']]) \
                    .dropna(subset=['mygene_entrez', 'mygene_sym'])

            debug and fold_change.to_csv("%s.table.csv" % debug)

            if not fold_change.empty:
                fold_change['direction'] = fold_change.log2foldChange.map(
                    lambda x: "up" if x > 0 else 'down')
                fold_change['subset'] = subset
                fold_change['gpl'] = gpl
                fold_change['gse'] = gse.name
                fold_change['perm'] = perm
                results_list.append(fold_change.reset_index())

    return pd.concat(results_list)


# COLUMNS = {
#     'sample__id': 'sample_id',
#     'sample__gsm_name': 'gsm_name',
#     'annotation': 'annotation',
#     'serie_annotation__series__id': 'series_id',
#     'serie_annotation__platform__id': 'platform_id',
#     'serie_annotation__platform__gpl_name': 'gpl_name',
#     'serie_annotation__tag__tag_name': 'tag_name',
# }

@log_durations(logger.debug)
def get_analysis_df(case_query, control_query, modifier_query=""):
    # Fetch all relevant data
    queries = [case_query, control_query, modifier_query]
    tokens = set(cat(re_all('[a-zA-Z]\w*', query) for query in queries))
    df = pd.read_sql_query('''
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

@log_durations(logger.debug)
def load_gse(df, series_id, impute = False):
    gse_name = series_gse_name(series_id)
    logger.debug('Loading data for %s, id = %d', gse_name, series_id)
    gpl2data = {}
    gpl2probes = {}

    for platform_id in df.query("""series_id == %s""" % series_id).platform_id.unique():
        gpl_name = platform_gpl_name(platform_id)
        gpl2data[gpl_name] = get_data(series_id, platform_id, impute=impute)
        # gpl2data[gpl_name].to_csv("%s.%s.csv"%(gse_name, gpl_name))
        gpl2probes[gpl_name] = get_probes(platform_id)
    samples = df.query('series_id == %s' % series_id)
    return Gse(gse_name, samples, gpl2data, gpl2probes)

def query_record(id, table, id_field="id"):
    sql = """select * from %s where %s """ % (table, id_field) + """= %s"""
    cursor.execute(sql, (id,))
    return cursor.fetchone()


# @make_lookuper
def series_gse_name(series_id):
    return query_record(series_id, "series")['gse_name']

    # return Series.objects.values_list('id', 'gse_name')

# @make_lookuper
def platform_gpl_name(platform_id):
    # return Platform.objects.values_list('id', 'gpl_name')
    return query_record(platform_id, "platform")['gpl_name']


def query_gsm_names(gse_name):
    sql = "select * from series inner join sample on series.id = series_id where gse_name = '%s'"%gse_name
    return list(pd.read_sql(sql, conn)['gsm_name'])
# def __getMatrixNumHeaderLines(inStream):
#     p = re.compile(r'^"ID_REF"')
#     for i, line in enumerate(inStream):
#         if p.search(line):
#             return i
#

# def matrix_filenames(series_id, platform_id):
#     gse_name = series_gse_name(series_id)
#     yield "%s/%s_series_matrix.txt.gz" % (gse_name, gse_name)
#
#     gpl_name = platform_gpl_name(platform_id)
#     yield "%s/%s-%s_series_matrix.txt.gz" % (gse_name, gse_name, gpl_name)


# def get_matrix_filename(series_id, platform_id):
#     filenames = list(matrix_filenames(series_id, platform_id))
#     mirror_filenames = (os.path.join(conf.SERIES_MATRIX_MIRROR, filename) for filename in filenames)
#     mirror_filename = first(filename for filename in mirror_filenames if os.path.isfile(filename))
#     if mirror_filename:
#         return mirror_filename
#
#     for filename in filenames:
#         logger.info('Loading URL %s...' % (conf.SERIES_MATRIX_URL + filename))
#         try:
#             res = urllib2.urlopen(conf.SERIES_MATRIX_URL + filename)
#         except urllib2.URLError:
#             pass
#         else:
#             mirror_filename = os.path.join(conf.SERIES_MATRIX_MIRROR, filename)
#             logger.info('Cache to %s' % mirror_filename)
#
#             directory = os.path.dirname(mirror_filename)
#             if not os.path.exists(directory):
#                 os.makedirs(directory)
#             with open(mirror_filename, 'wb') as f:
#                 shutil.copyfileobj(res, f)
#
#             return mirror_filename
#
#     raise LookupError("Can't find matrix file for series %s, platform %s"
#                       % (series_id, platform_id))


# @log_durations(logger.debug)
# def get_data(series_id, platform_id):
#     matrixFilename = get_matrix_filename(series_id, platform_id)
#     # setup data for specific platform
#     for attempt in (0, 1):
#         try:
#             headerRows = __getMatrixNumHeaderLines(gzip.open(matrixFilename))
#             na_values = ["null", "NA", "NaN", "N/A", "na", "n/a"]
#             data = pd.io.parsers.read_table(gzip.open(matrixFilename),
#                                             skiprows=headerRows,
#                                             index_col=["ID_REF"],
#                                             na_values=na_values,
#                                             lineterminator='\n',
#                                             engine='c')
#             # Drop last line
#             data.to_csv("data.1.csv")
#             data = data.drop(data.index[-1]).dropna()
#             break
#         except IOError as e:
#             # In case we have cirrupt file
#             logger.error("Failed loading %s: %s" % (matrixFilename, e))
#             os.remove(matrixFilename)
#             if attempt:
#                 raise
#             matrixFilename = get_matrix_filename(series_id, platform_id)
#
#     data.index = data.index.astype(str)
#     data.index.name = "probe"
#     data.to_csv("data.2.csv")
#
#     for column in data.columns:
#         data[column] = data[column].astype(np.float64)
#
#
#     return data


@log_durations(logger.debug)
def get_probes(platform_id):
    sql = "select * from platform_probe where platform_id = %s"
    return pd.read_sql(sql, conn, "probe", params=(platform_id,))
#     df = PlatformProbe.objects.filter(platform=platform_id).order_by('id').to_dataframe()
#     # df = db(Platform_Probe.platform_id == platform_id).select(processor=pandas_processor)
#     df.columns = [col.lower().replace("platform_probe.", "") for col in df.columns]
#     df.probe = df.probe.astype(str)  # must cast probes as str
#     df = df.set_index('probe')
#     return df


class Gse:
    def __init__(self, name, samples, gpl2data, gpl2probes):
        self.name = name
        self.samples = samples
        self.gpl2data = gpl2data
        self.gpl2probes = gpl2probes

    def __str__(self):
        return '<Gse %s>' % self.name


def get_full_meta(fold_change,  debug=False):
    debug and fold_change.to_csv("%s.fc.csv" % debug)
    all = []
    # i = 0
    all_gene_fold_change = fold_change.sort('p') \
        .drop_duplicates(['gse', 'gpl', 'mygene_sym', 'mygene_entrez', 'subset']) \
        .dropna()
    debug and all_gene_fold_change.to_csv("%s.geneestimates.csv" % debug)
    for group, gene_fold_change in all_gene_fold_change.groupby(['mygene_sym', 'mygene_entrez']):
        mygene_sym, mygene_entrez = group
        if len(gene_fold_change) == 1:
            continue
        metaAnalysis = get_gene_meta(gene_fold_change)
        metaAnalysis['caseDataCount'] = gene_fold_change['caseDataCount'].sum()
        metaAnalysis['controlDataCount'] = gene_fold_change['controlDataCount'].sum()
        metaAnalysis['mygene_sym'] = mygene_sym
        metaAnalysis['mygene_entrez'] = mygene_entrez
        all.append(metaAnalysis)

    full_meta_analysis = pd.DataFrame(all).set_index(['mygene_sym', 'mygene_entrez'])
    debug and full_meta_analysis.to_csv('full_meta_analysis.csv')
    full_meta_analysis['direction'] = full_meta_analysis['random_TE'].map(
        lambda x: "up" if x >= 0 else "down")

    return full_meta_analysis


# @dcache_new.cached
def get_gene_meta(gene_fold_change):
    return MetaAnalyser(gene_fold_change).get_results()



class MetaAnalyser:
    def isquared(self, H):
        """
        Calculate I-Squared.
        Higgins & Thompson (2002), Statistics in Medicine, 21, 1539-58.
        """
        if H.TE:
            t = lambda x: (x ** 2 - 1) / x ** 2
            return EasyDict(TE=t(H.TE), lower=t(H.lower), upper=t(H.upper))
        else:
            return EasyDict(TE=None, lower=None, upper=None)

    def calcH(self, Q, df, level):
        """
        Calculate H.
        Higgins & Thompson (2002), Statistics in Medicine, 21, 1539-58.
        """
        k = df + 1
        H = np.sqrt(Q / (k - 1))

        result = EasyDict(TE=None, lower=None, upper=None)
        if k > 2:
            if Q <= k:
                selogH = np.sqrt(1 / (2 * (k - 2)) * (1 - 1 / (3 * (k - 2) ** 2)))
            else:
                selogH = 0.5 * (np.log(Q) - np.log(k - 1)) / (np.sqrt(2 * Q) - np.sqrt(2 * k - 3))

            tres = self.getConfidenceIntervals(np.log(H), selogH, level)
            result = EasyDict(TE=1 if np.exp(tres.TE) < 1 else np.exp(tres.TE),
                              lower=1 if np.exp(tres.lower) < 1 else np.exp(tres.lower),
                              upper=1 if np.exp(tres.upper) < 1 else np.exp(tres.upper))
        return result


    norm_ppf = staticmethod(memoize(stats.norm.ppf))
    t_ppf = staticmethod(memoize(stats.t.ppf))

    def getConfidenceIntervals(self, TE, TE_se, level=.95, df=None):
        alpha = 1 - level
        zscore = TE / TE_se
        if not df:
            delta = self.norm_ppf(1 - alpha / 2) * TE_se
            lower = TE - delta
            upper = TE + delta
            pval = 2 * (1 - stats.norm.cdf(abs(zscore)))
        else:
            delta = self.t_ppf(1 - alpha / 2, df=df) * TE_se
            lower = TE - delta
            upper = TE + delta
            pval = 2 * (1 - stats.t.cdf(abs(zscore), df=df))

        result = EasyDict(TE=TE,
                          se=TE_se,
                          level=level,
                          df=df,
                          pval=pval,
                          zscore=zscore,
                          upper=upper,
                          lower=lower)

        return result

    @staticmethod
    def get_TE_se(gene_stats):
        # Convert to numpy arrays for speed
        caseSigma = gene_stats['caseDataSigma'].values
        caseCount = gene_stats['caseDataCount'].values
        controlSigma = gene_stats['controlDataSigma'].values
        controlCount = gene_stats['controlDataCount'].values

        # MD method
        na = np.sqrt(caseSigma ** 2 / caseCount + controlSigma ** 2 / controlCount)

        # Studies with non-positive variance get zero weight in meta-analysis
        for i in range(len(na)):
            if caseSigma[i] <= 0 or controlSigma[i] <= 0:
                na[i] = float('nan')

        return pd.Series(na, index=gene_stats.index)

    def __init__(self, gene_stats):
        TE = gene_stats.caseDataMu.values - gene_stats.controlDataMu.values

        # (7) Calculate results for individual studies
        TE_se = self.get_TE_se(gene_stats)
        w_fixed = (1 / TE_se ** 2).fillna(0)

        TE_fixed = np.average(TE, weights=w_fixed)
        TE_fixed_se = np.sqrt(1 / sum(w_fixed))
        self.fixed = self.getConfidenceIntervals(TE_fixed, TE_fixed_se)

        self.Q = sum(w_fixed * (TE - TE_fixed) ** 2)
        self.Q_df = TE_se.dropna().count() - 1

        self.cVal = sum(w_fixed) - sum(w_fixed ** 2) / sum(w_fixed)
        if self.Q <= self.Q_df:
            self.tau2 = 0
        else:
            self.tau2 = (self.Q - self.Q_df) / self.cVal
        self.tau = np.sqrt(self.tau2)
        self.tau2_se = None
        w_random = (1 / (TE_se ** 2 + self.tau2)).fillna(0)
        TE_random = np.average(TE, weights=w_random)
        TE_random_se = np.sqrt(1 / sum(w_random))
        self.random = self.getConfidenceIntervals(TE_random, TE_random_se)

        # Prediction interval
        self.level_predict = 0.95
        self.k = TE_se.count()
        self.predict = EasyDict(TE=None,
                                se=None,
                                level=None,
                                df=None,
                                pval=None,
                                zscore=None,
                                upper=None,
                                lower=None)
        if self.k >= 3:
            TE_predict_se = np.sqrt(TE_random_se ** 2 + self.tau2)
            self.predict = self.getConfidenceIntervals(TE_random, TE_predict_se, self.level_predict,
                                                       self.k - 2)

        # Calculate H and I-Squared
        self.level_comb = 0.95
        self.H = self.calcH(self.Q, self.Q_df, self.level_comb)
        self.I2 = self.isquared(self.H)

    def get_results(self):
        return dict(
            k=self.k,
            fixed_TE=self.fixed.TE,
            fixed_se=self.fixed.se,
            fixed_lower=self.fixed.lower,
            fixed_upper=self.fixed.upper,
            fixed_pval=self.fixed.pval,
            fixed_zscore=self.fixed.zscore,

            random_TE=self.random.TE,
            random_se=self.random.se,
            random_lower=self.random.lower,
            random_upper=self.random.upper,
            random_pval=self.random.pval,
            random_zscore=self.random.zscore,


            predict_TE=self.predict.TE,
            predict_se=self.predict.se,
            predict_lower=self.predict.lower,
            predict_upper=self.predict.upper,
            predict_pval=self.predict.pval,
            predict_zscore=self.predict.zscore,

            tau2=self.tau2,
            tau2_se=self.tau2_se,

            C=self.cVal,

            H=self.H.TE,
            H_lower=self.H.lower,
            H_upper=self.H.upper,

            I2=self.I2.TE,
            I2_lower=self.I2.lower,
            I2_upper=self.I2.upper,

            Q=self.Q,
            Q_df=self.Q_df,
        )

# @log_durations(logger.debug)
def get_fold_change(data, sample_class, debug=False):
    from scipy.stats import ttest_ind

    data = normalize_quantiles(log_data(data))

    summary = pd.DataFrame(index=data.index)

    summary['dataMu'] = data.mean(axis="columns")
    summary['dataSigma'] = data.std(axis="columns")
    summary['dataCount'] = data.count(axis="columns")

    caseData = data.T[sample_class == 1].T
    debug and caseData.to_csv("%s.case.data.csv" % debug)
    summary['caseDataMu'] = caseData.mean(axis="columns")
    summary['caseDataSigma'] = caseData.std(axis="columns") if len(caseData.columns) > 1 else 0
    summary['caseDataCount'] = caseData.count(axis="columns")

    controlData = data.T[sample_class == 0].T
    debug and controlData.to_csv("%s.control.data.csv" % debug)

    summary['controlDataMu'] = controlData.mean(axis="columns")
    summary['controlDataSigma'] = controlData.std(axis="columns") \
        if len(controlData.columns) > 1 else 0
    summary['controlDataCount'] = controlData.count(axis="columns")

    summary['fc'] = summary['caseDataMu'] - summary['controlDataMu']
    summary['log2foldChange'] = summary['fc']
    # else:
    # summary['fc'] = np.log2(summary['caseDataMu']/summary['controlDataMu'])

    summary['effect_size'] = summary['fc'] / summary['dataSigma']

    ttest, prob = ttest_ind(caseData, controlData, axis=1)
    summary['ttest'] = ttest
    summary['p'] = prob
    summary['direction'] = summary['effect_size'].map(lambda x: "up" if x >= 0 else "down")

    return summary


def normalize_quantiles(df):
    """
    df with samples in the columns and probes across the rows
    """
    # http://biopython.org/pipermail/biopython/2010-March/006319.html
    A = df.values
    AA = np.empty_like(A)
    I = np.argsort(A, axis=0)
    AA[I, np.arange(A.shape[1])] = np.mean(A[I, np.arange(A.shape[1])], axis=1)[:, np.newaxis]
    return pd.DataFrame(AA, index=df.index, columns=df.columns)

def query_median_gene_data(gse_name, gpl_name):
    import glob
    """returns the median intensity"""
    filename = os.path.join(conf.CSV_CACHE, "%s.%s.median.csv"%(gse_name, gpl_name))
    if glob.glob(filename):
        median_gene_data = pd.read_csv(filename)\
            .set_index(['mygene_sym', 'mygene_entrez'])
    else:
        gene_data = query_gene_data(gse_name, gpl_name)
        median_gene_data = gene_data \
            .reset_index() \
            .groupby(['mygene_sym', 'mygene_entrez']) \
            .median()
        median_gene_data.to_csv(filename)
    return median_gene_data

    # gene_index = gene_data_median.index
    # samples = gene_data_median.reset_index().T
    # samples['gpl_name'] = gpl_name
    # samples['gse_name'] = gse_name
    # gene_data = samples.reset_index().set_index(['gsm_name', 'gpl_name', 'gse_name']).T
    # gene_data.index = gene_index
    # return gene_data


def combine_matrix(names):
    """returns an combined matrix of genes vs samples expression values over all supplied gses"""
    i = 0
    gse_name, gpl_name = names[0]
    m = query_median_gene_data(gse_name, gpl_name)
    for (gse_name, gpl_name) in names[1:]:
        i += 1
        print "%s/%s"%(i, len(names)), gse_name, gpl_name
        median_gene_data = query_median_gene_data(gse_name, gpl_name)
        if median_gene_data.empty:
            continue
        m = m.join(median_gene_data,
                   how="outer")
    return m

def combine_samples(names):
    gse_name, gpl_name = names[0]
    combined_samples = query_samples(gse_name, gpl_name)
    combined_samples['gse_name'] = gse_name
    combined_samples['gpl_name'] = gpl_name
    for (gse_name, gpl_name) in names[1:]:
        print gse_name, gpl_name,
        samples = query_samples(gse_name, gpl_name)
        samples['gpl_name'] = gpl_name
        combined_samples = pd.concat([combined_samples, samples])
    return combined_samples


def save_upcs(gse_name):
    import random, glob
    for gsm_name in random.shuffle(gsm_names):
        upc = query_upc(gsm_name)


def query_upc(gsm_name):
    import glob
    import rpy2.rinterface
    filename = os.path.join(conf.CSV_CACHE, "%s.upc.csv"%(gsm_name))
    upc = None
    print "working on", filename
    if glob.glob(filename):
        upc = pd.read_csv(filename, index_col=0)
    else:
        try:
            print "STARTING R UPC"
            r.library("SCAN.UPC")
            r_matrix = r.exprs(r.UPC(gsm_name))
            upc = pd.DataFrame(np.asmatrix(r_matrix))
            upc.index = list(r_matrix.rownames)
            upc.columns = list(r_matrix.colnames)
            print "SAVING R UPC"
            upc.to_csv(filename)
        except rpy2.rinterface.RRuntimeError as e:
            print e.message

    return upc

def save_upcs(gsm_names):
    import random
    gsm_names = list(gsm_names)
    random.shuffle(gsm_names)
    for gsm_name in gsm_names:
        print gsm_name
        query_upc(gsm_name)

def combat(df):
    names = df[['gse_name', 'gpl_name']].drop_duplicates().to_records(index=False)
    # drop genes with missing data
    df['code'] = df.gsm_name + "_" + df.gpl_name + "_" + df.gse_name
    df = df.set_index('code')

    combined_matrix = combine_matrix(names)
    # combined_matrix.to_csv("combined_matrix.csv")
    m = drop_missing_samples(combined_matrix).dropna()
    # drop_missing_genes = drop_missing_genes(dropMissingSamples(combined_matrix)).dropna() #UNNECESSARY
    samples_m = df.index.intersection(m.columns)
    m = m[samples_m]
    m.to_csv("m.csv")
    samples = df \
        .ix[m.columns] \
        .reset_index()
    samples.to_csv("samples.csv")
    edata = com.convert_to_r_matrix(m)
    batch = robjects.StrVector(samples.gse_name + '_' +  samples.gpl_name)
    # pheno = robjects.FactorVector(samples.sample_class)
    pheno = robjects.FactorVector(samples.annotation)
    r.library("sva")
    fmla = robjects.Formula('~pheno')
    # fmla.environment['pheno'] = r['as.factor'](pheno)
    fmla.environment['pheno'] = pheno
    mod = r['model.matrix'](fmla)
    r_combat_edata = r.ComBat(dat=edata, batch=batch, mod=mod)
    combat_matrix = pd.DataFrame(np.asmatrix(r_combat_edata))
    combat_matrix.index = m.index
    combat_matrix.columns = m.columns
    return combat_matrix, samples
