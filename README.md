## Installation

For development, clone or download-and-extract the repository. Then run `pip install -e .` from the repository's root directory. The `-e` flag specifies [editable](https://pythonhosted.org/setuptools/setuptools.html#development-mode) mode, so updating the source updates your installation.

The current master branch can also be installed directly from GitHub using:

```
pip install git+git://github.com/idrdex/star_api.git@master
```

## Database configuration

Users must manually create a file at `starapi/db_conf.py` with the following:

```python
DB_PARAMATERS = """host='star.cwxuiazqnkkn.us-west-2.rds.amazonaws.com'
                   dbname='star'
                   user='XXXXXXXXXX'
                   password='XXXXXXXXXX'"""
```

Contact @idrdex for this information.


## Usage

```python
from starapi.main import get_annotations
from starapi.analysis import combat, perform_analysis
from starapi import conf
conf.configure('./')

analysis = EasyDict(
    analysis_name="severe_dengue_top",
    case_query="DHF=='DHF' or DSS=='DSS'",
    control_query="DF=='DF'",
    modifier_query="Dengue_Acute=='Dengue_Acute'",
    min_samples=3
)

nperm = 3
basename = "%s.%s_perm"%(analysis.analysis_name, nperm)

df, fc, results, permutations = perform_analysis(
    analysis=analysis,
    impute=False,
    nperm=nperm,
    mygene_filter=None,
    debug=basename
)
```
