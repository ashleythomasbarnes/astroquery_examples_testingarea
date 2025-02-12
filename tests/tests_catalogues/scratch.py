def _is_collection_and_table_list_at_eso(collections=None, tables=None, all_versions=False):
    r"""Check if lists of collections and tables are present in the ESO archive and merge them in a list of tables

    Args:
        collections (any): str or list of collections to be tested
        tables (any): `list` of table_names (or a single `str`) to be tested
        all_versions (bool): if set to `True` also obsolete versions of the catalogues are searched in case
            `collections` is given
    Returns:
        list: merge of `collections` and `tables` where collections or tables not present at ESO are removed

    """
    # test on collections
    clean_collections = _is_collection_list_at_eso(collections)
    # test on tables
    clean_tables = _is_table_list_at_eso(tables)
    if clean_tables is None:
        clean_tables = []
    if clean_collections is not None:
        for clean_collection in clean_collections:
            clean_tables += _get_tables_from_collection(clean_collection, all_versions=all_versions)
    # This removes possible duplicates and removes None
    clean_tables = list(filter(None, list(set(clean_tables))))
    return clean_tables


def _get_tables_from_collection(collection, all_versions=False):
    r"""Returns the table_name corresponding to a given collection

    Args:
        collection (str): name of the collection for which the tables will be extracted
        all_versions (bool): if set to `True` also obsolete versions of the catalogues are listed

    Returns:
        list: list containing the `table_name` corresponding to the selected `collection`

    """
    if not _is_collection_at_eso(collection):
        return None
    table_all_catalogues = all_catalogues_info(verbose=False, all_versions=all_versions)
    table_selected_catalogues = table_all_catalogues[(table_all_catalogues['collection'].data == collection)]
    list_selected_tables = table_selected_catalogues['table_name'].data.data.tolist()
    return list_selected_tables


def _get_catalogue_length_from_tables(tables, maxrec=None, all_versions=False):
    r"""Returns a list with the length of catalogues given in `tables`

    Args:
        tables (any): `list` of table_names (or a single `str`) to be queried
        all_versions (bool): if set to `True` also obsolete versions of the catalogues are listed
        maxrec (int, optional): define the maximum number of entries that a single query can return. If it is `None` the
            value is set by the limit of the service.

    Returns:
        list: list of `int` containing the length of each catalogue in input. If `maxrec` is set, it will return a
            with the same length of tables, but with all entries set to `maxrec`

    """
    if maxrec is None:
        maxrec_list = []
        for table in tables:
            maxrec_list.append(_get_catalogue_length_from_table(table, all_versions=all_versions))
    else:
        maxrec_list = [maxrec] * len(tables)
    return maxrec_list


def _get_catalogue_length_from_table(table_name, all_versions=False):
    r"""Returns the length of a catalogue given a `table_name`

    Args:
        table_name (str): name of the collection for which the tables will be extracted
        all_versions (bool): if set to `True` also obsolete versions of the catalogues are listed

    Returns:
        int: `number_rows` corresponding to the selected `table_name`

    """
    if not _is_table_at_eso(table_name):
        return None
    table_all_catalogues = all_catalogues_info(verbose=False, all_versions=all_versions)
    table_selected_catalogues = table_all_catalogues[(table_all_catalogues['table_name'].data == table_name)]
    selected_number_columns = int(table_selected_catalogues['number_rows'].data.data)
    return selected_number_columns


def _is_column_list_in_catalogues(columns, collections=None, tables=None):
    r"""Check if a given list of columns is present in the ESO archive

    It is possible to test for given collection (or table) by setting the appropriate values as input

    Args:
        columns (any): list of string containing the column_name (or the single `str`) to be tested
        collections (any): list of `str` containing the names of the collections (or a single `str`) from which the
            columns will be extracted
        tables (any): list of `str`  (or a single `str`) containing the names of the tables from which the columns
            will be extracted

    Returns:
        list: same of `columns` but columns not present at in the collections/tables are removed

    """
    assert columns is None or isinstance(columns, (str, list)), r'`columns` must be `None` or a `str` or a `list`'
    columns_list = cleaning_lists.from_element_to_list(columns, element_type=str)
    if columns is not None:
        # test if it is a valid column
        clean_columns = []
        for column in columns_list:
            if _is_column_in_catalogues(column, collections=collections, tables=tables):
                clean_columns.append(column)
    else:
        clean_columns = None
    return clean_columns


def _is_column_in_catalogues(column_name, collections=None, tables=None):
    r"""Check if a given column is present in the ESO archive

    Args:
        column_name (str): column to be tested
        collections (any): list of `str` containing the names of the collections (or a single `str`) from which the
            columns will be extracted
        tables (any): list of `str`  (or a single `str`) containing the names of the tables from which the columns
            will be extracted

    Returns:
        bool: `True` if the column is present in the selected collections/tables. `False` and warning raised
            otherwise

    """
    is_at_eso = True
    table_all_columns = columns_info(collections=collections, tables=tables, verbose=False)
    all_column_list = table_all_columns['column_name'].data.data.tolist()
    if column_name not in all_column_list:
        print('Column: {} not recognized. Possible values are:\n{}'.format(column_name, all_column_list))
        is_at_eso = False
    return is_at_eso