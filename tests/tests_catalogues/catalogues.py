import numpy as np
from astropy.table import MaskedColumn
from pyvo import dal
from pyvo.dal import DALQueryError, DALFormatError

# =============================================================================
# Constants
# =============================================================================
TAP_SERVICE_URL = "https://archive.eso.org/tap_cat"
TAP_QUERY_TYPES = ["sync", "async"]
MAXREC = 1000

# =============================================================================
# Public API Functions
# =============================================================================

def list_catalogues(all_versions=False, collections=None, tables=None, verbose=False):
    """
    Retrieve a table with metadata about ESO catalogues.
    
    The returned table includes information such as collection, title, version,
    table_name, instrument, telescope, publication_date, etc. In addition, the 
    RA, Dec, and Source ID column names are added to each catalogue.
    
    Args:
        all_versions (bool): If True, include obsolete catalogue versions.
        collections (str or list): Filter results by collection name(s).
        tables (str or list): Filter results by table name(s).
        verbose (bool): If True, print additional query info.
    
    Returns:
        astropy.table.Table: Catalogue metadata table.
    """
    clean_collections = _is_collection_list_at_eso(collections)
    clean_tables = _is_table_list_at_eso(tables)
    query = _create_query_all_catalogues(all_versions, clean_collections, clean_tables)
    
    if verbose:
        _print_query(query)
    
    qobj = _ESOCatalogues(query=query)
    if collections is not None and tables is not None:
        print("Warning: Both `collections` and `tables` are set. Ensure this is the intended behavior.")
    
    qobj.run_query(to_string=True)
    qobj.result.sort(["collection", "table_name", "version"])
    qobj.set_last_version(update=True)
    catalogues_table = qobj.get_result()
    
    # Add RA, Dec, and Source ID columns based on UCD tokens.
    id_ra_dec_table = _get_id_ra_dec_from_columns(clean_collections)
    source_id, ra_id, dec_id = [], [], []
    for t_name in catalogues_table["table_name"]:
        source_id_table = id_ra_dec_table[
            (id_ra_dec_table["table_name"] == t_name) &
            (id_ra_dec_table["ucd"] == "meta.id;meta.main")
        ]["column_name"].tolist()
        ra_id_table = id_ra_dec_table[
            (id_ra_dec_table["table_name"] == t_name) &
            (id_ra_dec_table["ucd"] == "pos.eq.ra;meta.main")
        ]["column_name"].tolist()
        dec_id_table = id_ra_dec_table[
            (id_ra_dec_table["table_name"] == t_name) &
            (id_ra_dec_table["ucd"] == "pos.eq.dec;meta.main")
        ]["column_name"].tolist()
        
        source_id.append(source_id_table[0] if len(source_id_table) == 1 else None)
        ra_id.append(ra_id_table[0] if len(ra_id_table) == 1 else None)
        dec_id.append(dec_id_table[0] if len(dec_id_table) == 1 else None)
    
    catalogues_table.add_column(
        MaskedColumn(data=np.asarray(ra_id), name="table_RA", dtype=str,
                     description="Identifier for RA in the catalog")
    )
    catalogues_table.add_column(
        MaskedColumn(data=np.asarray(dec_id), name="table_Dec", dtype=str,
                     description="Identifier for Dec in the catalog")
    )
    catalogues_table.add_column(
        MaskedColumn(data=np.asarray(source_id), name="table_ID", dtype=str,
                     description="Identifier for Source ID in the catalog")
    )
    return catalogues_table


def all_list_catalogues(all_versions=False, verbose=False):
    """
    Retrieve a master table with metadata on all ESO catalogues.
    
    This is equivalent to calling:
        list_catalogues(all_versions=all_versions, collections=None, tables=None)
    
    Args:
        all_versions (bool): If True, include obsolete catalogue versions.
        verbose (bool): If True, print additional query information.
    
    Returns:
        astropy.table.Table: Table containing metadata for all catalogues.
    """
    return list_catalogues(all_versions=all_versions, collections=None, tables=None, verbose=verbose)


def list_catalogues_info(collections=None, tables=None, verbose=False):
    """
    Retrieve column metadata for the ESO catalogues.
    
    The returned table includes column name, UCD, datatype, description, and unit.
    
    Args:
        collections (str or list): Filter by collection name(s).
        tables (str or list): Filter by table name(s).
        verbose (bool): If True, print additional query information.
    
    Returns:
        astropy.table.Table: Table with column metadata.
    """
    clean_collections = _is_collection_list_at_eso(collections)
    clean_tables = _is_table_list_at_eso(tables)
    query = _create_query_all_columns(clean_collections, clean_tables)
    
    if verbose:
        _print_query(query)
    
    qobj = _ESOCatalogues(query=query)
    if collections is not None and tables is not None:
        print("Warning: Both `collections` and `tables` are set. Ensure this is the intended behavior.")
    qobj.run_query(to_string=True)
    return qobj.get_result()


def query_catalogues(collections=None, tables=None, columns=None, type_of_query='sync',
                     all_versions=False, maxrec=None, verbose=False,
                     conditions_dict=None, top=None, order_by=None, order='ascending'):
    """
    Query specific ESO catalogues from the TAP service.
    
    You can either supply a collection (or list of collections) or specific table names.
    If both are provided, the conditions are combined (AND).
    
    Args:
        collections (str or list): Collection name(s) to filter catalogues.
        tables (str or list): Specific table name(s) to query.
        columns (str or list): Column name(s) to retrieve.
        type_of_query (str): 'sync' or 'async' query mode.
        all_versions (bool): If True, include obsolete catalogue versions.
        maxrec (int): Maximum number of rows to retrieve per query.
        verbose (bool): If True, print query details.
        conditions_dict (dict): Additional query conditions.
        top (int): Return only the top N rows.
        order_by (str): Column name for ordering the result.
        order (str): Order direction ('ascending' or 'descending').
    
    Returns:
        astropy.table.Table or list of Tables: The queried catalogue(s).
    """
    clean_tables = _is_collection_and_table_list_at_eso(collections, tables, all_versions=all_versions)
    totrec_list = _get_catalogue_length_from_tables(clean_tables, maxrec=None, all_versions=all_versions)
    maxrec_list = [maxrec] * len(totrec_list) if maxrec is not None else [MAXREC] * len(totrec_list)
    
    list_of_catalogues = []
    for table_name, totrec, maxrec_val in zip(clean_tables, totrec_list, maxrec_list):
        valid_columns = _is_column_list_in_catalogues(columns, tables=table_name)
        query = _create_query_catalogues(table_name, valid_columns, conditions_dict, order_by, order, top)
        
        if verbose:
            _print_query(query)
        
        qobj = _ESOCatalogues(query=query, type_of_query=type_of_query, maxrec=maxrec_val)
        qobj.run_query(to_string=True)
        catalogue = qobj.get_result()
        list_of_catalogues.append(catalogue)
        print(f"The query to {table_name} returned {len(catalogue)} entries out of {totrec} "
              f"(with a limit set to maxrec={maxrec_val})")
    
    if len(list_of_catalogues) == 0:
        return None
    elif len(list_of_catalogues) == 1:
        return list_of_catalogues[0]
    else:
        return list_of_catalogues

# =============================================================================
# Internal Implementation (hidden from the user)
# =============================================================================

class _ESOCatalogues:
    """
    Internal class to manage ESO TAP queries.
    """
    def __init__(self, query=None, type_of_query="sync", maxrec=None):
        self.tap_service = _define_tap_service()
        self.query = query
        self.type_of_query = type_of_query if type_of_query in TAP_QUERY_TYPES else "sync"
        self.maxrec = maxrec or MAXREC
        self.result = None

    def run_query(self, to_string=True):
        """Execute the query and (optionally) convert byte columns to strings."""
        self.result = _run_query(self.tap_service, self.query, self.type_of_query, self.maxrec)
        if to_string and self.result is not None:
            for col in self.result.colnames:
                self.result[col].data[:] = _from_bytes_to_string(self.result[col].data.data)

    def get_result(self):
        """Return a copy of the query result."""
        return self.result.copy() if self.result is not None else None

    def set_last_version(self, update=True):
        """
        Add a `last_version` column to the result table indicating whether 
        a catalogue is the most recent version.
        """
        required_cols = ["title", "version"]
        if self.result is None:
            print("No results to update.")
            return
        for col in required_cols:
            if col not in self.result.colnames:
                print(f"Column '{col}' missing. Cannot create 'last_version'.")
                return
        if "last_version" in self.result.colnames and not update:
            print("'last_version' already exists; skipping update.")
            return

        unique_titles = np.unique(self.result["title"].data).tolist()
        last_version_flags = np.zeros_like(self.result["version"].data, dtype=bool)

        for title in unique_titles:
            version_data = self.result["version"].data[self.result["title"].data == title]
            latest_version = np.nanmax(version_data)
            flag = (self.result["title"].data == title) & (self.result["version"].data == latest_version)
            last_version_flags[flag] = True

        self.result.add_column(
            MaskedColumn(
                data=last_version_flags,
                name="last_version",
                dtype=bool,
                description="True if this is the latest version of the catalogue"
            )
        )

# -----------------------------------------------------------------------------
# Internal helper functions
# -----------------------------------------------------------------------------

def _define_tap_service():
    """Instantiate and return the TAP service."""
    return dal.tap.TAPService(TAP_SERVICE_URL)

def _run_query(tap_service, query, type_of_query, maxrec=MAXREC):
    """Dispatch the query to the appropriate synchronous or asynchronous function."""
    if type_of_query == "sync":
        return _run_query_sync(tap_service, query, maxrec)
    else:
        return _run_query_async(tap_service, query, maxrec)

def _run_query_sync(tap_service, query, maxrec=MAXREC):
    """Execute a synchronous TAP query."""
    try:
        return tap_service.search(query=query, maxrec=int(maxrec) if maxrec is not None else None).to_table()
    except (ValueError, DALQueryError, DALFormatError):
        print("Query timeout. Retrying with maxrec=100 (consider using async instead).")
        return tap_service.search(query=query, maxrec=100).to_table()

def _run_query_async(tap_service, query, maxrec=MAXREC):
    """Execute an asynchronous TAP query."""
    tap_job = tap_service.submit_job(query=query, maxrec=maxrec)
    tap_job.run()
    for status in ["EXECUTING", "COMPLETED", "ERROR", "ABORTED"]:
        tap_job.wait(phases=[status], timeout=10.0)
        print(f"Query status: {tap_job.phase}")
    tap_job.raise_if_error()
    return tap_job.fetch_result().to_table()

def _from_bytes_to_string(input_in_bytes):
    """Convert byte strings to unicode strings."""
    if isinstance(input_in_bytes, bytes):
        return input_in_bytes.decode("utf-8")
    if isinstance(input_in_bytes, np.ndarray) and input_in_bytes.dtype.type is np.bytes_:
        return np.array([x.decode("utf-8") if isinstance(x, bytes) else x for x in input_in_bytes])
    if isinstance(input_in_bytes, list):
        return [x.decode("utf-8") if isinstance(x, bytes) else x for x in input_in_bytes]
    return input_in_bytes

def _from_element_to_list(element, element_type=str):
    """
    Ensure the input is a list of elements of a given type.
    """
    if element is None:
        return None
    if isinstance(element, list):
        if all(isinstance(e, element_type) for e in element):
            return element
        else:
            raise TypeError(f"All elements must be of type {element_type}")
    if isinstance(element, np.ndarray):
        lst = element.tolist()
        if all(isinstance(e, element_type) for e in lst):
            return lst
        else:
            raise TypeError(f"All elements must be of type {element_type}")
    if isinstance(element, MaskedColumn):
        lst = element.data.data.tolist()
        if all(isinstance(e, element_type) for e in lst):
            return lst
        else:
            raise TypeError(f"All elements must be of type {element_type}")
    if isinstance(element, element_type):
        return [element]
    raise TypeError(f"Invalid type for: {element} (expected {element_type})")

def _print_query(query):
    """Print the query string."""
    if query:
        print(f"Query:\n{query}")
    else:
        print("The query is empty.")

def _create_comma_separated_list(list_of_strings):
    """Return a comma-separated string from a list (or '*' if None)."""
    return ", ".join(list_of_strings) if list_of_strings else "*"

def _condition_collections_like(collections):
    """Generate a SQL-like condition for collections."""
    cols = _from_element_to_list(collections, str) if collections is not None else ["%"]
    return " OR ".join(f"collection LIKE '{col}'" for col in cols)

def _condition_tables_like(tables):
    """Generate a SQL-like condition for tables."""
    tabs = _from_element_to_list(tables, str) if tables is not None else ["%"]
    return " OR ".join(f"table_name LIKE '{tab}'" for tab in tabs)

def _condition_order_by_like(order_by, order="ascending"):
    """Generate an ORDER BY clause if needed."""
    return f" ORDER BY {order_by} {order.upper()}" if order_by else ""

def _conditions_dict_like(conditions_dict):
    """Generate a WHERE clause from a dictionary of conditions."""
    if conditions_dict:
        return " ".join(f"WHERE {key}={value}" for key, value in conditions_dict.items())
    return ""

def _create_query_all_catalogues(all_versions, collections, tables):
    """Build the TAP query for retrieving catalogue metadata."""
    query = """
        SELECT 
            collection, title, version, table_name, filter, instrument, telescope, publication_date, 
            ref.description AS description, number_rows, number_columns, rel_descr_url, acknowledgment,
            cat_id, mjd_obs, mjd_end, skysqdeg, bibliography, document_id, kc.from_column AS from_column,
            k.target_table AS target_table, kc.target_column AS target_column, schema_name
        FROM TAP_SCHEMA.tables AS ref
        LEFT OUTER JOIN TAP_SCHEMA.keys AS k ON ref.table_name = k.from_table 
        LEFT OUTER JOIN TAP_SCHEMA.key_columns AS kc ON k.key_id = kc.key_id
        WHERE schema_name = 'safcat'
    """
    if not all_versions:
        query += """
        AND cat_id IN (
            SELECT t1.cat_id 
            FROM TAP_SCHEMA.tables t1
            LEFT JOIN TAP_SCHEMA.tables t2 ON (t1.title = t2.title AND t1.version < t2.version)
            WHERE t2.title IS NULL
        )
        """
    if collections:
        query += f" AND ({_condition_collections_like(collections)})"
    if tables:
        query += f" AND ({_condition_tables_like(tables)})"
    return query

def _create_query_all_columns(collections, tables):
    """Build the TAP query for retrieving column metadata."""
    return f"""
        SELECT table_name, column_name, ucd, datatype, description, unit
        FROM TAP_SCHEMA.columns
        WHERE table_name IN (
            SELECT table_name FROM TAP_SCHEMA.tables WHERE {_condition_collections_like(collections)}
        )
        AND ({_condition_tables_like(tables)})
    """

def _create_query_catalogues(table_name, columns, conditions_dict, order_by, order, top):
    """Build the TAP query for a specific catalogue table."""
    base = _create_query_table_base(table_name, columns, top)
    cond = _conditions_dict_like(conditions_dict)
    order_clause = _condition_order_by_like(order_by, order)
    return f"{base} {cond} {order_clause}"

def _create_query_table_base(table_name, columns, top):
    """Build the basic SELECT ... FROM ... part of a query."""
    select_clause = f"SELECT {'TOP ' + str(top) + ' ' if top else ''}{_create_comma_separated_list(columns)}"
    return f"{select_clause} FROM {table_name}"

def _is_collection_at_eso(collection):
    """Check if the collection exists in the ESO archive."""
    table = all_list_catalogues(all_versions=False, verbose=False)
    all_cols = np.unique(table["collection"].data.data).tolist()
    if collection not in all_cols:
        print(f"Warning: Collection '{collection}' not recognized. Possible values:\n{all_cols}")
        return False
    return True

def _is_table_at_eso(table_name):
    """Check if the table exists (and if it is the latest version)."""
    table = all_list_catalogues(all_versions=True, verbose=False)
    all_tables = table["table_name"].data.data.tolist()
    last_versions = table["last_version"].data.data.tolist()
    if table_name not in all_tables:
        print(f"Warning: Table '{table_name}' not recognized. Possible values:\n{all_tables}")
        return False
    if not last_versions[all_tables.index(table_name)]:
        print(f"Warning: '{table_name}' is not the most recent version of the catalogue.")
    return True

def _is_collection_list_at_eso(collections):
    """Ensure the collections input is valid and filter to known collections."""
    if collections is None:
        return None
    collections_list = _from_element_to_list(collections, str)
    return [c for c in collections_list if _is_collection_at_eso(c)]

def _is_table_list_at_eso(tables):
    """Ensure the tables input is valid and filter to known tables."""
    if tables is None:
        return None
    tables_list = _from_element_to_list(tables, str)
    return [t for t in tables_list if _is_table_at_eso(t)]

def _get_tables_from_collection(collection, all_versions=False):
    """Retrieve table names for a given collection."""
    if not _is_collection_at_eso(collection):
        return []
    table_all = all_list_catalogues(all_versions=all_versions, verbose=False)
    if table_all is None or "table_name" not in table_all.colnames:
        return []
    return table_all[table_all["collection"].data == collection]["table_name"].tolist()

def _is_collection_and_table_list_at_eso(collections, tables, all_versions=False):
    """Merge valid table names from collections and explicit table inputs."""
    clean_collections = _is_collection_list_at_eso(collections)
    clean_tables = _is_table_list_at_eso(tables) or []
    if clean_collections:
        for coll in clean_collections:
            clean_tables += _get_tables_from_collection(coll, all_versions=all_versions)
    return list(set(filter(None, clean_tables)))

def _get_catalogue_length_from_table(table_name, all_versions=False):
    """Return the number of rows for a given table."""
    if not _is_table_at_eso(table_name):
        return None
    table_all = all_list_catalogues(all_versions=all_versions, verbose=False)
    if table_all is None or "number_rows" not in table_all.colnames:
        return None
    sel = table_all[table_all["table_name"].data == table_name]
    return int(sel["number_rows"].data[0]) if len(sel) else None

def _get_catalogue_length_from_tables(tables, maxrec=None, all_versions=False):
    """Return a list of row counts (or maxrec if set) for each table."""
    if not tables:
        return []
    if maxrec is not None:
        return [maxrec] * len(tables)
    return [_get_catalogue_length_from_table(t, all_versions=all_versions) for t in tables]

def _is_column_in_catalogues(column_name, collections=None, tables=None):
    """Check if a given column exists in the catalogues."""
    table_all = list_catalogues_info(collections=collections, tables=tables, verbose=False)
    if table_all is None or "column_name" not in table_all.colnames:
        return False
    return column_name in table_all["column_name"].data.tolist()

def _is_column_list_in_catalogues(columns, collections=None, tables=None):
    """Filter a list of column names to those that exist in the catalogues."""
    if columns is None:
        return None
    columns_list = _from_element_to_list(columns, str)
    return [col for col in columns_list if _is_column_in_catalogues(col, collections, tables)]

def _get_id_ra_dec_from_columns(collections=None):
    """
    Extract the column names for Source ID, RA, and Dec based on UCD tokens.
    """
    all_columns_table = list_catalogues_info(collections, verbose=False)
    filter_tokens = ((all_columns_table["ucd"].data == "meta.id;meta.main") |
                     (all_columns_table["ucd"].data == "pos.eq.ra;meta.main") |
                     (all_columns_table["ucd"].data == "pos.eq.dec;meta.main"))
    return all_columns_table[filter_tokens]