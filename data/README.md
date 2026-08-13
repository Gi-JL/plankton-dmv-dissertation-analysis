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

Source: European Environment Agency

https://www.eea.europa.eu/data-and-maps/data/eea-coastline-for-analysis-1

The analysis scripts use the EEA coastline for analysis (polygon), version 3.0, March 2017.

The shapefile is expected as:

`EEA_Coastline_20170228.shp`

Keep all associated shapefile components together in the `data/` directory.

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
