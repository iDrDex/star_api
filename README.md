## installation

For development, clone or download-and-extract the repository. Then run `pip install -e .` from the repository's root directory. The `-e` flag specifies [editable](https://pythonhosted.org/setuptools/setuptools.html#development-mode) mode, so updating the source updates your installation.

## database configuration

Users must manually create a file at `starapi/db_conf.py` with the following:

```python
DB_PARAMATERS = """host='star.cwxuiazqnkkn.us-west-2.rds.amazonaws.com'
                   dbname='star'
                   user='XXXXXXXXXX'
                   password='XXXXXXXXXX'"""
```

Contact @idrdex for this information.
