# Input data

The analysis scripts in this repository require survey, discrete choice experiment, coded discussion, and geospatial input data. Most input data are not stored directly in the Git repository and must be obtained from the sources listed below.

Download the required external files and place them in this `data/` directory before running the analyses. The coded discussion data required for Paper C are included directly in this directory.

## Survey and discrete choice experiment data

### Paper A

Required files:

- `DMV data.csv`
- `final design_baesyian efficient design with interaction.NGD`

Source:

https://doi.org/10.5281/zenodo.12638010

### Paper B

Required files:

- `DMV data.csv`
- `PCODE_2020_PT_SH.shp` and the associated shapefile components
- `Europe_coastline.shp` and the associated shapefile components

### Paper C

Required files:

- `DMV data_EN.csv`
- `final design_baesyian efficient design with interaction.NGD`
- `Pro_contra_statements_respondents_csv_UTF_8.csv`

The first two files are available from:

https://doi.org/10.5281/zenodo.17710200

The coded discussion data file `Pro_contra_statements_respondents_csv_UTF_8.csv` is included directly in this repository.
## Geospatial data

Some analyses in Papers A and B calculate respondents' distance from the coastline and therefore require external geospatial datasets.

### European postal-code data

Source: Eurostat GISCO

https://ec.europa.eu/eurostat/web/gisco/geodata/administrative-units/postal-codes

Download the 2020 postal-code point dataset in ESRI Shapefile format using the EPSG:3035 projection.

The analysis scripts expect the shapefile to be available as:

`PCODE_PT_2020_3035.shp`

Keep all associated shapefile components extracted from the downloaded archive together in the `data/` directory.

### European coastline data

The distance-to-coast analyses in Papers A and B use the following dataset:

**EEA coastline for analysis (line), version 2.0, September 2015**

Original source: European Environment Agency (EEA)

https://www.eea.europa.eu/en/datahub/datahubitem-view/af40333f-9e94-4926-a4f0-0a787f1d2b8f?activeAccordion=278520

The line dataset used for the analyses has since been superseded by a newer
polygon version. At the time this repository was prepared, the original EEA
download for the superseded line dataset was no longer functioning reliably.

To preserve reproducibility, the original coastline shapefile used for the
analyses is therefore included directly in this repository under the following
filenames:

- `Europe_coastline.shp`
- `Europe_coastline.shx`
- `Europe_coastline.dbf`
- `Europe_coastline.prj`

These files originate from the European Environment Agency and are redistributed
here solely to reproduce the analyses reported in the dissertation. The EEA should
be acknowledged as the original source. These files are not covered by the
GPL-3.0-or-later license applying to the analysis code in this repository.

## Expected directory contents

After downloading the required external files, the `data/` directory should contain approximately:

```text
data/
├── README.md
├── DMV data.csv
├── DMV data_EN.csv
├── final design_baesyian efficient design with interaction.NGD
├── Pro_contra_statements_respondents_csv_UTF_8.csv
├── PCODE_2020_PT_SH.shp
├── PCODE_2020_PT_SH.dbf
├── PCODE_2020_PT_SH.shx
├── PCODE_2020_PT_SH.prj
├── Europe_coastline.shp
├── Europe_coastline.dbf
├── Europe_coastline.shx
└── Europe_coastline.prj
