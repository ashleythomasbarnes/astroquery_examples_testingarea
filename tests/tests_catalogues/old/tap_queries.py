from pyvo import dal
from pyvo.dal import DALQueryError, DALFormatError
import old.default as default

# Initialize default values
default = default.Default()

# Supported TAP services
TAP_SERVICES = ["eso_tap_cat", "eso_tap_obs"]

# Allowed query types
TAP_QUERY_TYPES = ["sync", "async"]

# Default columns for ObsCore queries
COLUMNS_FROM_OBSCORE = [
    "target_name", "dp_id", "s_ra", "s_dec", "t_exptime", "em_min", "em_max",
    "dataproduct_type", "instrument_name", "obstech", "abmaglim", "proposal_id", "obs_collection"
]


def define_tap_service(which_tap_service):
    """Load a TAP service from defaults."""
    if which_tap_service not in TAP_SERVICES:
        print(f"Invalid TAP service: {which_tap_service}. Options: {TAP_SERVICES}")
    return dal.tap.TAPService(default.get_value(which_tap_service))


def which_service(tap_service):
    """Print a summary of the TAP service in use."""
    print("TAP Service Description:")
    tap_service.describe()


def run_query(tap_service, query, type_of_query, maxrec=default.get_value("maxrec")):
    """Run a query on a TAP service and return the result as an astropy Table."""
    if type_of_query not in TAP_QUERY_TYPES:
        print(f"Invalid query type: {type_of_query}. Options: {TAP_QUERY_TYPES}")

    if not query:
        print("Empty query provided.")
        return None

    return run_query_sync(tap_service, query, maxrec) if type_of_query == "sync" else run_query_async(tap_service, query, maxrec)


def run_query_sync(tap_service, query, maxrec=default.get_value("maxrec")):
    """Execute a synchronous TAP query."""
    try:
        return tap_service.search(query=query, maxrec=int(maxrec) if maxrec is not None else None).to_table()
    except (ValueError, DALQueryError, DALFormatError):
        print("Query timeout. Retrying with maxrec=100 (consider using async instead).")
        return tap_service.search(query=query, maxrec=100).to_table()


def run_query_async(tap_service, query, maxrec=default.get_value("maxrec")):
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