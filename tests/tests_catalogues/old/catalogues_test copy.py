import numpy as np
from astropy.table import MaskedColumn


class ESOCatalogues:
    """
    Class for querying ESO scientific catalogues via TAP services.

    This class provides a standardized interface for querying the ESO archive,
    supporting both synchronous ('sync') and asynchronous ('async') queries.

    Attributes:
        tap_service (pyvo.dal.tap.TAPService): TAP service for querying ESO catalogues.
        query (str): The query string.
        type_of_query (str): Type of query ('sync' or 'async').
        maxrec (int, optional): Maximum number of entries returned by the query.
        result_from_query (astropy.table.Table): The result of the query.
    """

    def __init__(self, query=None, result_from_query=None, type_of_query="sync", maxrec=None):
        """Initialize an ESO catalog query with default TAP service settings."""
        self.tap_service = define_tap_service("eso_tap_cat")
        self.query = query
        self.result_from_query = result_from_query
        self.maxrec = maxrec or MAXREC

        # Validate query type
        if type_of_query not in TAP_QUERY_TYPES:
            print(f"Invalid TAP query type: {type_of_query}. Using 'sync' instead.")
            self.type_of_query = "sync"
        else:
            self.type_of_query = type_of_query

    # Query Execution Methods
    def run_query(self, to_string=True):
        """
        Execute the query and store results in result_from_query.

        Args:
            to_string (bool, optional): If True, converts columns in byte format to strings.
        """
        self.result_from_query = run_query(self.tap_service, self.query, self.type_of_query, maxrec=self.maxrec)

        # Convert bytes to strings if required
        if to_string and self.result_from_query is not None:
            for column_id in self.which_columns():
                self.result_from_query[column_id].data[:] = from_bytes_to_string(self.result_from_query[column_id].data.data)

    def print_query(self):
        """Print the current query."""
        print_query(self.query)

    def which_service(self):
        """Return and describe the TAP service in use."""
        return which_service(self.tap_service)

    def which_columns(self):
        """Return a list of column names in result_from_query."""
        if self.result_from_query is None:
            print("No results found. Run the query first.")
            return []
        return self.result_from_query.colnames

    def get_result_from_query(self):
        """Return a copy of result_from_query."""
        return self.result_from_query.copy() if self.result_from_query else None

    def clean_query(self):
        """Reset the query attribute to None."""
        self.query = None

    def clean_result_from_query(self):
        """Reset the result_from_query attribute to None."""
        self.result_from_query = None

    # ESO-Specific Functionality
    def set_last_version(self, update=True):
        """
        Adds a `last_version` column to result_from_query.

        The `last_version` column is a boolean indicating whether a catalogue version is the latest.

        Args:
            update (bool, optional): If True, updates the column if it already exists.
        """
        required_columns = ["title", "version"]

        # Ensure required columns exist
        for col in required_columns:
            if col not in self.which_columns():
                print(f"Column '{col}' missing. `last_version` will not be created.")
                return

        # Prevent duplicate column creation
        if "last_version" in self.which_columns():
            print(f"'last_version' column already exists. {'Updating' if update else 'Skipping update'}.")
            if not update:
                return

        # Determine latest versions
        unique_titles = np.unique(self.result_from_query["title"].data).tolist()
        last_version_flags = np.zeros_like(self.result_from_query["version"].data, dtype=bool)

        for title in unique_titles:
            latest_version = np.nanmax(self.result_from_query["version"].data[self.result_from_query["title"].data == title])
            last_version_flags[(self.result_from_query["title"].data == title) & 
                               (self.result_from_query["version"].data == latest_version)] = True

        # Add column to table
        self.result_from_query.add_column(
            MaskedColumn(
                data=last_version_flags,
                name="last_version",
                dtype=bool,
                description="True if this is the latest version of the catalog"
            )
        )


def query_catalogues(collections=None, tables=None, columns=None, type_of_query='sync', all_versions=False, maxrec=None, verbose=False,
                   conditions_dict=None, top=None, order_by=None, order='ascending'):
    r"""Query the ESO tap_cat service for specific catalogues

    There are two ways to select the catalogues you are interested in. Either you select directly the table_name (or the
    list of table_names) that you want to query, or you select a collection (or a list of collections). If you select
    this latter option, what happens in the background is that the code is going to search for the table(s)
    corresponding to the given collection and query them.

    If you are asking for more than one table, the result will be listed in a list of `astropy.tables` with one element
    per retrieved table

    Args:
        collections (any): list of `str` containing the names (or a single `str`) of the collections for
            which the query will be limited
        tables (any): list of `str` containing the table_name of the tables for which the query will be limited
        columns (any): list of the `column_name` that you want to download. The full list of the columns in a
            table can be found by running `list_catalogues_info()`
        all_versions (bool): if set to `True` also obsolete versions of the catalogues are searched in case
            `collections` is given
        type_of_query (str): type of query to be run
        maxrec (int, optional): define the maximum number of entries that a single query can return. If it is `None` the
            value is set by the limit of the service.
        verbose (bool): if set to `True` additional info will be displayed
        conditions_dict (dict): dictionary containing the conditions to be applied to the query
            Only currently supported is "=" condition, e.g. ID = 1234 - WORK IN PROGRESS
        top (int): number of top rows to be returned
        order_by (str): column name to be used to order the query
        order (str): order of the query (ascending or descending)

    Returns:
        any: `astropy.table` or `list` of `astropy.tables` containing the queried catalogues

    """
    # Obtain list of all tables derived from the merger of collections and tables
    clean_tables = _is_collection_and_table_list_at_eso(collections=collections, tables=tables,
                                                        all_versions=all_versions)

    # Get total number of records for each table
    totrec_list = _get_catalogue_length_from_tables(clean_tables, maxrec=None, all_versions=all_versions)

    # if maxrec is set to None, the MAXREC is used
    if maxrec is None:
        maxrec_list = [MAXREC] * len(totrec_list)
    else: 
        maxrec_list = [maxrec]

    list_of_catalogues = []
    for table_name, totrec_for_table, maxrec_for_table in zip(clean_tables, totrec_list, maxrec_list):

        # test for columns
        columns_in_table = _is_column_list_in_catalogues(columns, tables=table_name)

        # Form query
        query = create_query_catalogues(table_name, columns_in_table, conditions_dict, order_by, order, top)

        # Print query
        if verbose: 
            print_query(query)

        # instantiate ESOcatalogues
        query_table = ESOCatalogues(query=query,
                                    type_of_query=type_of_query, 
                                    maxrec=maxrec_for_table)
        
        query_table.run_query(to_string=True)
        catalogue = query_table.get_result_from_query()
        list_of_catalogues.append(catalogue)
        print('The query to {} returned {} entries out of {} (with a limit set to maxrec={})'.format(table_name,
                                                                                               len(catalogue),
                                                                                               totrec_for_table,
                                                                                               maxrec_for_table))
    if len(list_of_catalogues) == 0:
        return None
    elif len(list_of_catalogues) == 1:
        return list_of_catalogues[0]
    else:
        return list_of_catalogues


def all_list_catalogues(all_versions=False, verbose=False):
    """
    Loads an `astropy.table.Table` containing information on all catalogues present in the ESO archive.

    The returned table includes metadata about all available catalogues in the ESO archive, such as collection names, 
    instruments, telescope details, publication dates, and available data formats.

    The table contains the following columns:
        - `collection`, `title`, `version`, `table_name`, `filter`, `instrument`, `telescope`, `publication_date`
        - `description`, `number_rows`, `number_columns`, `rel_descr_url`, `acknowledgment`, `cat_id`, `mjd_obs`
        - `mjd_end`, `skysqdeg`, `bibliography`, `document_id`, `from_column`, `target_table`, `target_column`
        - `last_version`

    **Note:** This is equivalent to running:
    
        >>> list_catalogues(collections=None, tables=None)
    
    The key difference is that, since no constraints are set, this function returns the **master table** 
    containing all catalogues present in the ESO archive.

    Args:
        all_versions (bool, optional): If `True`, includes obsolete versions of the catalogues. Default is `False`.
        verbose (bool, optional): If `True`, displays additional information during execution. Default is `False`.

    Returns:
        astropy.table.Table: A table containing metadata about all catalogues in the ESO archive.
    """
    return list_catalogues(all_versions=all_versions, collections=None, tables=None, verbose=verbose)


def list_catalogues(all_versions=False, collections=None, tables=None, verbose=False):
    """
    Retrieves an `astropy.table.Table` containing metadata about selected catalogues in the ESO archive.

    Users can filter catalogues by providing either a **list of collections** or a **list of table names**. 
    If both `collections` and `tables` are set to `None`, the function returns metadata for **all ESO catalogues**.

    The output table contains the following columns:
        - `collection`, `title`, `version`, `table_name`, `filter`, `instrument`, `telescope`, `publication_date`
        - `description`, `number_rows`, `number_columns`, `rel_descr_url`, `acknowledgment`, `cat_id`, `mjd_obs`
        - `mjd_end`, `skysqdeg`, `bibliography`, `document_id`, `from_column`, `target_table`, `target_column`
        - `last_version`, `table_RA`, `table_Dec`, `table_ID`

    **Note:** If both `collections` and `tables` are provided, the conditions are **combined with an AND statement**, 
    which may result in unexpected behavior. 

    Args:
        all_versions (bool, optional): If `True`, includes obsolete versions of the catalogues. Default is `False`.
        collections (str or list, optional): A collection name or list of collection names to filter results.
        tables (str or list, optional): A table name or list of table names to filter results.
        verbose (bool, optional): If `True`, displays additional query information. Default is `False`.

    Returns:
        astropy.table.Table: A table containing metadata about the selected catalogues.
    """

    # Validate collection and table inputs
    clean_collections = _is_collection_list_at_eso(collections)
    clean_tables = _is_table_list_at_eso(tables)

    # Form query
    query = create_query_all_catalogues(all_versions=all_versions, collections=clean_collections, tables=clean_tables)

    # Print query
    if verbose: 
        print_query(query)

    # Create a query for catalogues
    query_for_catalogues = ESOCatalogues(query=query)

    # Warn if both collections and tables are provided
    if collections is not None and tables is not None:
        print("Warning: Both `collections` and `tables` are set. Ensure this is the intended behavior.")

    # Execute the query
    query_for_catalogues.run_query(to_string=True)

    # Sort results by collection, table name, and version
    query_for_catalogues.result_from_query.sort(["collection", "table_name", "version"])

    # Mark the latest versions, redundant if `all_versions=True`
    query_for_catalogues.set_last_version(update=True)
    catalogues_table = query_for_catalogues.get_result_from_query()

    # Extract Source ID, RA, and Dec columns for all collections
    id_ra_dec_table = _get_id_ra_dec_from_columns(collections=clean_collections)

    # Initialize lists for RA, Dec, and Source ID
    source_id, ra_id, dec_id = [], [], []

    # Iterate over tables to extract RA, Dec, and Source ID information
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

        # Handle cases with multiple or missing columns
        ra_id.append(ra_id_table[0] if len(ra_id_table) == 1 else None)
        dec_id.append(dec_id_table[0] if len(dec_id_table) == 1 else None)
        source_id.append(source_id_table[0] if len(source_id_table) == 1 else None)

    # Add RA, Dec, and Source ID columns to the result table
    catalogues_table.add_column(MaskedColumn(data=np.asarray(ra_id), name="table_RA", dtype=str,
                                             description="Identifier for RA in the catalog"))
    catalogues_table.add_column(MaskedColumn(data=np.asarray(dec_id), name="table_Dec", dtype=str,
                                             description="Identifier for Dec in the catalog"))
    catalogues_table.add_column(MaskedColumn(data=np.asarray(source_id), name="table_ID", dtype=str,
                                             description="Identifier for Source ID in the catalog"))

    return catalogues_table


def list_catalogues_info(collections=None, tables=None, verbose=False):
    """
    Queries the ESO archive for column metadata of selected collections or tables.

    If `collections` and `tables` are both `None`, the function retrieves metadata for all collections 
    and tables in the ESO archive.

    Args:
        collections (str or list, optional): A collection name or list of collections to filter results.
        tables (str or list, optional): A table name or list of tables to filter results.
        verbose (bool, optional): If `True`, displays additional query details.

    Returns:
        astropy.table.Table: A table containing metadata of all columns in the queried collection(s) or table(s).
        The output includes:
        - `table_name`, `column_name`, `ucd`, `datatype`, `description`, and `unit`
    """
    # Validate inputs
    clean_collections = _is_collection_list_at_eso(collections)
    clean_tables = _is_table_list_at_eso(tables)

    # Form query
    query = create_query_all_columns(collections=clean_collections, tables=clean_tables)

    # Print query
    if verbose: 
        print_query(query)

    # Create the query
    query_all_columns_info = ESOCatalogues(query=query)

    # Warn if both collections and tables are provided
    if collections is not None and tables is not None:
        print("Warning: Both `collections` and `tables` are set. Ensure this is the intended behavior.")

    # Print query if verbose mode is enabled
    if verbose:
        query_all_columns_info.print_query()

    # Execute the query
    query_all_columns_info.run_query(to_string=True)

    return query_all_columns_info.get_result_from_query()

def _is_collection_list_at_eso(collections):
    """
    Checks if a given list of collections is present in the ESO archive.

    Args:
        collections (str or list, optional): A single collection name or a list of collection names.

    Returns:
        list: A filtered list containing only valid collections present in the ESO archive.
    """
    assert collections is None or isinstance(collections, (str, list)), "`collections` must be None, str, or list"

    collections_list = from_element_to_list(collections, element_type=str)

    return [c for c in collections_list if _is_collection_at_eso(c)] if collections_list else None


def _is_table_list_at_eso(tables):
    """
    Checks if a given list of table names is present in the ESO archive.

    Args:
        tables (str or list, optional): A single table name or a list of table names.

    Returns:
        list: A filtered list containing only valid tables present in the ESO archive.
    """
    assert tables is None or isinstance(tables, (str, list)), "`tables` must be None, str, or list"

    tables_list = from_element_to_list(tables, element_type=str)

    return [t for t in tables_list if _is_table_at_eso(t)] if tables_list else None


def _is_collection_at_eso(collection):
    """
    Checks if a given collection exists in the ESO archive.

    Args:
        collection (str): The collection name to check.

    Returns:
        bool: True if the collection exists, False otherwise (with a warning message).
    """
    table_all_catalogues = all_list_catalogues(verbose=False, all_versions=False)
    all_collections_list = np.unique(table_all_catalogues["collection"].data.data).tolist()

    if collection not in all_collections_list:
        print(f"Warning: Collection '{collection}' not recognized. Possible values:\n{all_collections_list}")
        return False

    return True


def _is_table_at_eso(table_name):
    """
    Checks if a given table exists in the ESO archive.

    Args:
        table_name (str): The table name to check.

    Returns:
        bool: True if the table exists, False otherwise (with a warning message).
    """
    table_all_catalogues = all_list_catalogues(verbose=False, all_versions=True)
    all_table_list = table_all_catalogues["table_name"].data.data.tolist()
    last_version_list = table_all_catalogues["last_version"].data.data.tolist()

    if table_name not in all_table_list:
        print(f"Warning: Table '{table_name}' not recognized. Possible values:\n{all_table_list}")
        return False

    if not last_version_list[all_table_list.index(table_name)]:
        print(f"Warning: '{table_name}' is not the most recent version of the queried catalogue.")

    return True


def _get_id_ra_dec_from_columns(collections=None):
    """
    Extracts the column names corresponding to Source ID, RA, and DEC from a list of collections.

    This is based on the following UCD tokens:
        - `meta.id;meta.main` -> Source ID
        - `pos.eq.ra;meta.main` -> RA
        - `pos.eq.dec;meta.main` -> Dec

    Args:
        collections (str or list, optional): A single collection name or a list of collections.

    Returns:
        astropy.table.Table: A table containing the column names for Source ID, RA, and Dec.
    """
    all_columns_table = list_catalogues_info(collections)
    filter_tokens = (
        (all_columns_table["ucd"].data == "meta.id;meta.main") |
        (all_columns_table["ucd"].data == "pos.eq.ra;meta.main") |
        (all_columns_table["ucd"].data == "pos.eq.dec;meta.main")
    )
    return all_columns_table[filter_tokens]

def _is_collection_and_table_list_at_eso(collections=None, tables=None, all_versions=False):
    """
    Checks if lists of collections and tables exist in the ESO archive and merges them into a list of valid table names.

    Args:
        collections (str or list, optional): A collection name or list of collections to check.
        tables (str or list, optional): A table name or list of tables to check.
        all_versions (bool, optional): If `True`, includes obsolete catalogue versions when retrieving tables 
            from a given collection.

    Returns:
        list: A merged list of valid table names, with duplicates and invalid entries removed.
    """
    clean_collections = _is_collection_list_at_eso(collections)
    clean_tables = _is_table_list_at_eso(tables) or []

    if clean_collections:
        for collection in clean_collections:
            clean_tables += _get_tables_from_collection(collection, all_versions=all_versions)

    return list(filter(None, set(clean_tables)))  # Remove duplicates and None values


def _get_tables_from_collection(collection, all_versions=False):
    """
    Retrieves table names associated with a given collection in the ESO archive.

    Args:
        collection (str): The name of the collection for which tables will be retrieved.
        all_versions (bool, optional): If `True`, includes obsolete versions of catalogues.

    Returns:
        list: A list of table names corresponding to the given collection, or an empty list if no tables are found.
    """
    if not _is_collection_at_eso(collection):
        return []

    table_all_catalogues = all_list_catalogues(verbose=False, all_versions=all_versions)
    
    if table_all_catalogues is None or "table_name" not in table_all_catalogues.colnames:
        return []

    return table_all_catalogues[table_all_catalogues["collection"].data == collection]["table_name"].tolist()


def _get_catalogue_length_from_tables(tables, maxrec=None, all_versions=False):
    """
    Retrieves the number of rows for each catalogue in the provided list of tables.

    Args:
        tables (str or list): A table name or list of table names to query.
        all_versions (bool, optional): If `True`, includes obsolete versions of the catalogues.
        maxrec (int, optional): Defines the maximum number of entries a single query can return. If `None`,
            the service's default limit is used.

    Returns:
        list: A list of integers representing the row count for each table. If `maxrec` is set, all values 
        will be equal to `maxrec`.
    """
    if not tables:
        return []

    if maxrec is not None:
        return [maxrec] * len(tables)

    return [_get_catalogue_length_from_table(table, all_versions=all_versions) for table in tables]


def _get_catalogue_length_from_table(table_name, all_versions=False):
    """
    Retrieves the number of rows for a specific catalogue table.

    Args:
        table_name (str): The name of the table to query.
        all_versions (bool, optional): If `True`, includes obsolete versions of the catalogues.

    Returns:
        int: The number of rows in the given table. Returns `None` if the table does not exist.
    """
    if not _is_table_at_eso(table_name):
        return None

    table_all_catalogues = all_list_catalogues(verbose=False, all_versions=all_versions)
    
    if table_all_catalogues is None or "number_rows" not in table_all_catalogues.colnames:
        return None  # Ensuring robustness

    table_selected_catalogues = table_all_catalogues[table_all_catalogues['table_name'].data == table_name]

    if len(table_selected_catalogues) == 0:
        return None

    return int(table_selected_catalogues['number_rows'].data[0])


def _is_column_list_in_catalogues(columns, collections=None, tables=None):
    """
    Checks whether a given list of columns exists in the ESO archive.

    This function verifies if each column in the provided list exists in the ESO database.
    If `columns` is `None`, it returns `None` (keeping original behavior).

    Args:
        columns (str or list, optional): A column name or list of column names to validate.
        collections (str or list, optional): A collection name or list of collections to restrict the search.
        tables (str or list, optional): A table name or list of tables to restrict the search.

    Returns:
        list or None: A list of valid column names found in the specified collections or tables.
                     Returns `None` if `columns` is `None` (same as original function).
    """
    if columns is None:
        return None  # Preserve original behavior

    assert isinstance(columns, (str, list)), "`columns` must be `None`, a `str`, or a `list`"

    # Convert single string column to a list
    columns_list = from_element_to_list(columns, element_type=str)

    # Filter out invalid columns
    clean_columns = []
    for column in columns_list:
        if _is_column_in_catalogues(column, collections=collections, tables=tables):
            clean_columns.append(column)

    return clean_columns


def _is_column_in_catalogues(column_name, collections=None, tables=None):
    """
    Checks whether a given column exists in the specified collections or tables in the ESO archive.

    Args:
        column_name (str): The column name to check.
        collections (str or list, optional): Collection name(s) to restrict the search.
        tables (str or list, optional): Table name(s) to restrict the search.

    Returns:
        bool: `True` if the column exists, `False` otherwise.
    """
    table_all_columns = list_catalogues_info(collections=collections, tables=tables, verbose=False)

    if table_all_columns is None or "column_name" not in table_all_columns.colnames:
        return False  # Preserve original behavior

    all_column_list = table_all_columns["column_name"].data.tolist()

    return column_name in all_column_list  # Direct return without unnecessary variable assignment


from pyvo import dal
from pyvo.dal import DALQueryError, DALFormatError
import numpy as np
from astropy.table import MaskedColumn

# Supported TAP services
TAP_SERVICE = "https://archive.eso.org/tap_cat"

# Allowed query types
TAP_QUERY_TYPES = ["sync", "async"]

# Default columns for ObsCore queries
COLUMNS_FROM_OBSCORE = [
    "target_name", "dp_id", "s_ra", "s_dec", "t_exptime", "em_min", "em_max",
    "dataproduct_type", "instrument_name", "obstech", "abmaglim", "proposal_id", "obs_collection"
]

# Maximum number of records to retrieve in a single query
MAXREC = 1000

# Minimum disk space required for the query (Gb)
MIN_DISK_SPACE = 6.00

###########################
# CHECKS 
###########################

def _is_collection_list_at_eso(collections):
    """
    Checks if a given list of collections is present in the ESO archive.

    Args:
        collections (str or list, optional): A single collection name or a list of collection names.

    Returns:
        list: A filtered list containing only valid collections present in the ESO archive.
    """
    assert collections is None or isinstance(collections, (str, list)), "`collections` must be None, str, or list"

    collections_list = from_element_to_list(collections, element_type=str)

    return [c for c in collections_list if _is_collection_at_eso(c)] if collections_list else None


def _is_table_list_at_eso(tables):
    """
    Checks if a given list of table names is present in the ESO archive.

    Args:
        tables (str or list, optional): A single table name or a list of table names.

    Returns:
        list: A filtered list containing only valid tables present in the ESO archive.
    """
    assert tables is None or isinstance(tables, (str, list)), "`tables` must be None, str, or list"

    tables_list = from_element_to_list(tables, element_type=str)

    return [t for t in tables_list if _is_table_at_eso(t)] if tables_list else None


def _is_collection_at_eso(collection):
    """
    Checks if a given collection exists in the ESO archive.

    Args:
        collection (str): The collection name to check.

    Returns:
        bool: True if the collection exists, False otherwise (with a warning message).
    """
    table_all_catalogues = all_catalogues_info(verbose=False, all_versions=False)
    all_collections_list = np.unique(table_all_catalogues["collection"].data.data).tolist()

    if collection not in all_collections_list:
        print(f"Warning: Collection '{collection}' not recognized. Possible values:\n{all_collections_list}")
        return False

    return True


def _is_table_at_eso(table_name):
    """
    Checks if a given table exists in the ESO archive.

    Args:
        table_name (str): The table name to check.

    Returns:
        bool: True if the table exists, False otherwise (with a warning message).
    """
    table_all_catalogues = all_catalogues_info(verbose=False, all_versions=True)
    all_table_list = table_all_catalogues["table_name"].data.data.tolist()
    last_version_list = table_all_catalogues["last_version"].data.data.tolist()

    if table_name not in all_table_list:
        print(f"Warning: Table '{table_name}' not recognized. Possible values:\n{all_table_list}")
        return False

    if not last_version_list[all_table_list.index(table_name)]:
        print(f"Warning: '{table_name}' is not the most recent version of the queried catalogue.")

    return True



def _get_catalogue_length_from_tables(tables, maxrec=None, all_versions=False):
    """
    Retrieves the number of rows for each catalogue in the provided list of tables.

    Args:
        tables (str or list): A table name or list of table names to query.
        all_versions (bool, optional): If `True`, includes obsolete versions of the catalogues.
        maxrec (int, optional): Defines the maximum number of entries a single query can return. If `None`,
            the service's default limit is used.

    Returns:
        list: A list of integers representing the row count for each table. If `maxrec` is set, all values 
        will be equal to `maxrec`.
    """
    if not tables:
        return []

    if maxrec is not None:
        return [maxrec] * len(tables)

    return [_get_catalogue_length_from_table(table, all_versions=all_versions) for table in tables]


def _get_catalogue_length_from_table(table_name, all_versions=False):
    """
    Retrieves the number of rows for a specific catalogue table.

    Args:
        table_name (str): The name of the table to query.
        all_versions (bool, optional): If `True`, includes obsolete versions of the catalogues.

    Returns:
        int: The number of rows in the given table. Returns `None` if the table does not exist.
    """
    if not _is_table_at_eso(table_name):
        return None

    table_all_catalogues = all_catalogues_info(verbose=False, all_versions=all_versions)
    
    if table_all_catalogues is None or "number_rows" not in table_all_catalogues.colnames:
        return None  # Ensuring robustness

    table_selected_catalogues = table_all_catalogues[table_all_catalogues['table_name'].data == table_name]

    if len(table_selected_catalogues) == 0:
        return None

    return int(table_selected_catalogues['number_rows'].data[0])


def _is_column_list_in_catalogues(columns, collections=None, tables=None):
    """
    Checks whether a given list of columns exists in the ESO archive.

    This function verifies if each column in the provided list exists in the ESO database.
    If `columns` is `None`, it returns `None` (keeping original behavior).

    Args:
        columns (str or list, optional): A column name or list of column names to validate.
        collections (str or list, optional): A collection name or list of collections to restrict the search.
        tables (str or list, optional): A table name or list of tables to restrict the search.

    Returns:
        list or None: A list of valid column names found in the specified collections or tables.
                     Returns `None` if `columns` is `None` (same as original function).
    """
    if columns is None:
        return None  # Preserve original behavior

    assert isinstance(columns, (str, list)), "`columns` must be `None`, a `str`, or a `list`"

    # Convert single string column to a list
    columns_list = from_element_to_list(columns, element_type=str)

    # Filter out invalid columns
    clean_columns = []
    for column in columns_list:
        if _is_column_in_catalogues(column, collections=collections, tables=tables):
            clean_columns.append(column)

    return clean_columns


def _is_column_in_catalogues(column_name, collections=None, tables=None):
    """
    Checks whether a given column exists in the specified collections or tables in the ESO archive.

    Args:
        column_name (str): The column name to check.
        collections (str or list, optional): Collection name(s) to restrict the search.
        tables (str or list, optional): Table name(s) to restrict the search.

    Returns:
        bool: `True` if the column exists, `False` otherwise.
    """
    table_all_columns = list_catalogues_info(collections=collections, tables=tables, verbose=False)

    if table_all_columns is None or "column_name" not in table_all_columns.colnames:
        return False  # Preserve original behavior

    all_column_list = table_all_columns["column_name"].data.tolist()

    return column_name in all_column_list  # Direct return without unnecessary variable assignment

###########################
# CLEANING LISTS
###########################

def from_element_to_list(element, element_type=str):
    """
    Converts an input element into a list while ensuring all elements are of the specified type.

    This function ensures that:
    - If `element` is already a list, it is returned as is (after type validation).
    - If `element` is a NumPy array, it is converted to a list.
    - If `element` is an Astropy MaskedColumn, its data is extracted and converted to a list.
    - If `element` is a single value of the expected type, it is wrapped in a list.
    - If `element` is `None`, it returns `None`.

    Args:
        element (any): The input element that will be converted into a list.
        element_type (type, optional): The expected type of elements in the list. Defaults to `str`.

    Returns:
        list or None: A list containing the `element`, or `None` if `element` is `None`.
    """
    if element is None:
        return None

    if isinstance(element, list):
        assert all(isinstance(e, element_type) for e in element), f"All elements must be of type {element_type}"
        return element

    if isinstance(element, np.ndarray):
        element_list = element.tolist()
        assert all(isinstance(e, element_type) for e in element_list), f"All elements must be of type {element_type}"
        return element_list

    if isinstance(element, MaskedColumn):
        element_list = element.data.data.tolist()
        assert all(isinstance(e, element_type) for e in element_list), f"All elements must be of type {element_type}"
        return element_list

    if isinstance(element, element_type):
        return [element]

    raise TypeError(f"Invalid type for: {element} (expected {element_type})")


def from_bytes_to_string(input_in_bytes):
    """
    Converts an input from `bytes` to `str`. 

    This function handles:
    - A single bytes input (decodes to a string).
    - A NumPy array containing bytes (element-wise decoding).
    - A list containing bytes (element-wise decoding).
    - Any other input type remains unchanged.

    Args:
        input_in_bytes (any): The input value, which may be in bytes format.

    Returns:
        any: The input converted to a string (if applicable), otherwise returned unchanged.
    """
    if isinstance(input_in_bytes, bytes):
        return input_in_bytes.decode("utf-8")

    if isinstance(input_in_bytes, np.ndarray) and input_in_bytes.dtype.type is np.bytes_:
        return np.array([x.decode("utf-8") if isinstance(x, bytes) else x for x in input_in_bytes])

    if isinstance(input_in_bytes, list):
        return [x.decode("utf-8") if isinstance(x, bytes) else x for x in input_in_bytes]

    return input_in_bytes  # Return unchanged if not bytes-related


###########################
# TAP QUERIES
###########################


def define_tap_service(which_tap_service):
    """Load a TAP service from defaults."""
    return dal.tap.TAPService(TAP_SERVICE)


def which_service(tap_service):
    """Print a summary of the TAP service in use."""
    print("TAP Service Description:")
    tap_service.describe()


def run_query(tap_service, query, type_of_query, maxrec=MAXREC):
    """Run a query on a TAP service and return the result as an astropy Table."""
    if type_of_query not in TAP_QUERY_TYPES:
        print(f"Invalid query type: {type_of_query}. Options: {TAP_QUERY_TYPES}")

    if not query:
        print("Empty query provided.")
        return None

    return run_query_sync(tap_service, query, maxrec) if type_of_query == "sync" else run_query_async(tap_service, query, maxrec)


def run_query_sync(tap_service, query, maxrec=MAXREC):
    """Execute a synchronous TAP query."""
    try:
        return tap_service.search(query=query, maxrec=int(maxrec) if maxrec is not None else None).to_table()
    except (ValueError, DALQueryError, DALFormatError):
        print("Query timeout. Retrying with maxrec=100 (consider using async instead).")
        return tap_service.search(query=query, maxrec=100).to_table()


def run_query_async(tap_service, query, maxrec=MAXREC):
    """Execute an asynchronous TAP query."""
    tap_job = tap_service.submit_job(query=query, maxrec=maxrec)
    tap_job.run()

    for status in ["EXECUTING", "COMPLETED", "ERROR", "ABORTED"]:
        tap_job.wait(phases=[status], timeout=10.0)
        print(f"Query status: {tap_job.phase}")

    tap_job.raise_if_error()
    return tap_job.fetch_result().to_table()


def _create_comma_separated_list(list_of_strings):
    """Convert a list of strings into a comma-separated string."""
    return ", ".join(list_of_strings) if list_of_strings else "*"


def print_query(query):
    """Print the query or warn if it is empty."""
    if query:
        print(f"Query:\n{query}")
    else:
        print("The query is empty.")


def create_query_all_catalogues(all_versions=False, collections=None, tables=None):
    """Generate a TAP query to retrieve all catalogues from the ESO archive."""
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
        query += f" AND ({condition_collections_like(collections)})"
    if tables:
        query += f" AND ({condition_tables_like(tables)})"

    return query


def create_query_all_columns(collections=None, tables=None):
    """Generate a TAP query to retrieve column info from specified collections or tables."""
    return f"""
        SELECT table_name, column_name, ucd, datatype, description, unit
        FROM TAP_SCHEMA.columns
        WHERE table_name IN (
            SELECT table_name FROM TAP_SCHEMA.tables WHERE {condition_collections_like(collections)}
        )
        AND ({condition_tables_like(tables)})
    """


def create_query_catalogues(table_name, columns_in_table=None, conditions_dict=None, order_by=None, order="ascending", top=None):
    """Generate a TAP query to retrieve selected columns from a table."""
    query = "{0}{1}{2}".format(create_query_table_base(table_name, columns=columns_in_table, top=top),
                                    conditions_dict_like(conditions_dict),
                                    condition_order_by_like(order_by, order))
    return query


def create_query_table_base(table_name, columns=None, top=None):
    """Generate a TAP query to retrieve selected columns from a table."""
    select_clause = f"SELECT {f'TOP {top} ' if top else ''}{_create_comma_separated_list(columns)}"
    return f"{select_clause} FROM {table_name}"


def condition_collections_like(collections=None):
    """Generate a `LIKE` condition for a list of collections."""
    return " OR ".join(f"collection LIKE '{col}'" for col in (collections or ["%"]))


def condition_tables_like(tables=None):
    """Generate a `LIKE` condition for a list of tables."""
    return " OR ".join(f"table_name LIKE '{table}'" for table in (tables or ["%"]))


def condition_order_by_like(order_by=None, order="ascending"):
    """Generate an ORDER BY clause."""
    return f" ORDER BY {order_by} {order.upper()}" if order_by else ""


def conditions_dict_like(conditions_dict=None):
    """Generate WHERE conditions from a dictionary."""
    return " ".join(f"WHERE {key}={value}" for key, value in (conditions_dict or {}).items())