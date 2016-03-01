from easydict import EasyDict
import pandas as pd

# queries = pd.read_table("https://raw.githubusercontent.com/dhimmel/stargeo/master/data/files.tsv")
# for i, row in queries.iterrows():
#     if  row['balanced_permutation.tsv.gz']:
#        continue
#     analysis = EasyDict(
#         analysis_name = row['slim_name'],
#             case_query = row['case_query'],
#             control_query = row['control_query'],
#             modifier_query = "",
#             min_samples = 3
#         )
#     print analysis
#     fc, results, permutations = perform_analysis(analysis=analysis, nperm=0)
#     results.to_csv("%s.csv"%row['slim_name'])

# 1/0

# tokens = "MB_Group4","MB_Group3", "MB_SHH", "MB_WNT", "MB_unlabeled"
# labels = query_tags_annotations(tokens)
# save_upcs(labels.gsm_name.unique().tolist())
# 1/0
# print query_upc("GSM555237")
# 1/0
# tokens = "MB_Group4","MB_Group3", "MB_SHH", "MB_WNT", "MB_unlabeled"
# labels = query_tags_annotations(tokens)
# labels = get_unique_annotations(labels)
# combat_matrix, samples  = combat(labels.tail(300))
# 1/0

from starapi.main import get_annotations
from starapi.analysis import combat, perform_analysis
from starapi import conf
conf.configure('./')


# tokens = "MB_Group4","MB_Cerebellum_Control", "GBM_samples", "GBM_controls"
# labels = query_tags_annotations(tokens)
# labels = get_annotations("""DHF=='DHF' or DSS=='DSS'""",
#                 """DF=='DF'""",
#                          """Dengue_Acute=="Dengue_Acute" or Dengue_Early_Acute=='Dengue_Early_Acute' or Dengue_Late_Acute == 'Dengue_Late_Acute' or Dengue_DOF < 10""")



# labels['annotation'] = None
# for label in "df", "dhf", "dss":
#     labels.annotation[labels[label]==label]=label

# labels['code'] = labels.gsm_name + "_" + labels.gpl_name + "_" + labels.gse_name
# labels = labels.set_index('code')

# combat_matrix, samples  = combat(labels)
# combat_matrix.to_csv("combat_dengue.csv")
# # combat_matrix.to_csv("combat_brain.csv")
# samples.to_csv("combat_dengue_samples.csv")

# 1/0

# analysis = EasyDict(
#     analysis_name = "test",
#     case_query = """ Senescent=='Senescent'""",
#     control_query = """Senescent_control=='Senescent_control'""",
#     modifier_query = "",
#     min_samples = 3
# )


# analysis = EasyDict(
#     analysis_name = "test",
#     case_query = """Smoker == 'Smoker'""",# MS == 'MS'""",
#     control_query = """Nonsmoker == 'Nonsmoker'""",#"""MS_control == 'MS_control'""",
#     modifier_query = "",
#     min_samples = 3
# )

# analysis = EasyDict(
#     analysis_name = "test",
#     case_query = """ PHT == 'PHT' or hypertension == 'hypertension' """,
#     control_query = """PHT_Control == 'PHT_Control' or hypertension_control == 'hypertension_control'""",
#     modifier_query = "",
#     min_samples = 3
# )
#
# analysis = EasyDict(
#     analysis_name = "test",
#     case_query = """ T1D == 'T1D' """,
#     control_query = """T1D_Control == 'T1D_Control'""",
#     modifier_query = "",
#     min_samples = 3
# )

analysis = EasyDict(
    analysis_name = "severe_dengue_top",
    case_query = """DHF=='DHF' or DSS=='DSS'""",
    control_query = """DF=='DF'""",
    modifier_query = """Dengue_Acute=="Dengue_Acute" or Dengue_Early_Acute=='Dengue_Early_Acute' or Dengue_Late_Acute == 'Dengue_Late_Acute' or Dengue_DOF <= 7""",
    min_samples = 3
)

# dengue_100_perm = pd.read_csv("severe_dengue_top.1000_perm.results.csv", dtype = dict(mygene_entrez=int))
# dengue_100_perm['random_TE_abs'] = dengue_100_perm.random_TE.abs()
# dengue_100_perm['fixed_TE_abs'] = dengue_100_perm.fixed_TE.abs()
# mygene_filter = dengue_100_perm\
#     .query("""random_rank == 0 and fixed_rank == 0 """)\
#     .set_index(['mygene_sym', 'mygene_entrez'])\
#     .index.unique()
# print "mygene_filter of", len(mygene_filter), "genes"
mygene_filter = None

nperm = 3
basename = "%s.%s_perm"%(analysis.analysis_name, nperm)

fc, results, permutations, _ = perform_analysis(analysis=analysis,
                                             impute=False,
                                             nperm=nperm,
                                             mygene_filter=mygene_filter,
                                             debug=basename)
if results is not None:
    results.to_csv(basename + ".results.csv")
    # fc.to_csv(basename + ".fc.csv")
if permutations is not None:
    permutations.to_csv(basename + ".perm.csv")
