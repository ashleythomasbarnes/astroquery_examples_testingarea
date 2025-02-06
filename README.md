<h1 style="text-align: center;">
  <img src="http://archive.eso.org/i/esologo.png" alt="ESO Logo" width="50" style="vertical-align: middle;">
  ESO Science Archive - Jupyter Notebooks
</h1>

<!-- # ESO Science Archive - (astroquery) Notebooks -->

This repository contains Jupyter notebooks designed to introduce and interact with the [ESO Science Archive](https://archive.eso.org/cms.html). These notebooks specifically demonstrate the usage of the `astroquery.eso` module for querying, retrieving, and analyzing astronomical data from ESO's extensive archive.

## About the ESO Science Archive

The [ESO Science Archive](https://archive.eso.org/cms.html) is a vast repository of observational data collected by ESO telescopes, including VLT, ALMA, and La Silla facilities. It currently holds over 4.2 million products, including spectra, images, and raw observational data. Researchers and astronomers can access this archive to retrieve scientifically valuable datasets for analysis and discovery.

## Getting Started

### Prerequisites
Ensure you have the following installed:
- Python (>=3.9)
- Jupyter Notebook or Jupyter Lab
- Required Python libraries: `astroquery`
- Recommended Python libraries: `numpy`, `matplotlib`, `astropy`, `astroquery` 

You can install these dependencies using:
```sh
pip install numpy matplotlib astropy astroquery pandas
```

### Running the Notebook
To run the notebook locally:
```sh
jupyter notebook 01_ESO_notebook_introduction.ipynb
```
This will open the Jupyter interface in your web browser, allowing you to execute code cells interactively.

## Usage
These notebooks demonstrate how to:
- Query the ESO Science Archive using **astroquery.eso**
- Retrieve metadata and images for astronomical objects
- Analyze and visualize retrieved data using **matplotlib** and **astropy**
- Work with different types of data products, such as spectra and FITS images

## Contributing Data via ESO Phase 3

The scientific community plays a vital role in expanding and enriching the ESO Science Archive. ESO provides the **Phase 3 process**, a structured way for researchers to submit their processed data products, ensuring long-term accessibility and usability for the broader astronomical community. 

By sharing your data via Phase 3, you contribute to the collective scientific effort, allowing other researchers to build upon your work and make new discoveries. To learn more about the Phase 3 submission process, visit the [ESO Phase 3 guidelines](https://www.eso.org/sci/observing/phase3.html).

## Contributing
Contributions are welcome! If you'd like to contribute, please fork the repository and submit a pull request with your improvements.

## License
This project is licensed under the MIT License - see the LICENSE file for details.
