import os

class Default:
    """
    A class to load and manage default values from a configuration file.

    This version is independent of ESOAsg and allows users to specify any file path.
    """

    def __init__(self, file_path="default.txt"):
        """
        Initializes the Default class and loads values from a given file.

        Parameters
        ----------
        file_path : str, optional
            The path to the default values file (default is 'default.txt').
        """
        self.default_dict = {}
        self.file_path = file_path  # Store the file path
        self._load_from_file()

    def _load_from_file(self):
        """
        Loads default values from the specified file into a dictionary.

        Only reads lines that are not empty and do not start with '#'.
        Expected format: key==value (or key::value after processing).
        """
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"Default file not found: {self.file_path}")

        with open(self.file_path, "r", encoding="utf-8") as f:
            default_list = [
                line.strip().replace("==", "::")
                for line in f
                if not line.strip().startswith("#") and line.strip() != ""
            ]

        for default in default_list:
            try:
                default_quantity, default_value = default.split("::", 1)
                self.default_dict[default_quantity.strip()] = default_value.strip()
            except ValueError:
                print(f"Warning: Skipping malformed line in {self.file_path}: {default}")

    def get_value(self, card_name):
        """
        Retrieves the default value associated with a given key.

        Parameters
        ----------
        card_name : str
            The key whose default value is needed.

        Returns
        -------
        str
            The default value of `card_name`.

        Raises
        ------
        ValueError
            If the key does not exist in the dictionary.
        """
        if not isinstance(card_name, str):
            raise TypeError(f"Expected a string key, got {type(card_name)} instead.")

        if card_name in self.default_dict:
            return self.default_dict[card_name]
        else:
            raise ValueError(f"Key '{card_name}' not found in the default dictionary.")