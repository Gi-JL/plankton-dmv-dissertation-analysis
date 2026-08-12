# Input data

The analysis scripts in this repository require survey, discrete choice experiment, and geospatial input data that are not stored directly in the Git repository.

Download the required files from the sources listed below and place them in this `data/` directory before running the analyses.

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

Source:

https://doi.org/10.5281/zenodo.17710200

## Geospatial data

Some analyses in Papers A and B calculate respondents' distance from the coastline and therefore require external geospatial datasets.

### European postal-code data

Source: Eurostat GISCO

https://ec.europa.eu/eurostat/web/gisco/geodata/administrative-units/postal-codes

The analysis scripts expect the postal-code shapefile to be available as:

`PCODE_2020_PT_SH.shp`

Remember that a shapefile consists of several files with the same base name (for example `.shp`, `.shx`, `.dbf`, and `.prj`). Keep these files together in the `data/` directory.

### European coastline data

Source: European Environment Agency

https://www.eea.europa.eu/data-and-maps/data/eea-coastline-for-analysis-1

The analysis scripts expect the coastline shapefile to be available as:

`Europe_coastline.shp`

Again, retain all associated shapefile components in the `data/` directory.

## Expected directory contents

After downloading the required files, the directory should contain approximately:

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
