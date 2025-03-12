
def table_to_csv(pyvo_table, filepath):
    """
    Saves the table and prints the filename to stdout
    """
    pyvo_table.write(filepath, overwrite=True)
    print(f"Table saved to {filepath}")
