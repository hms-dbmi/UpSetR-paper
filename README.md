# UpSetR Paper Figures

## Pulling Project Mutation Data From ICGC
To pull mutation data from various projects on the ICGC data portal use the `ICGC_API_pull.R` script. Since the ICGC REST API has a pagination of 100 results, this script may take several minutes to run. It is also important to note that this API is actively being updated, so data is constantly being added and removed. To avoid running this script we've provided a text file (`icgcData.txt`) containing all of the data pulled on 3/18/2016, which is the data used to generate the UpSetR plots in the paper. 

## Generating ICGC plots
The `ICGCpaperFigures.R` script can be used to generate the ICGC UpSetR plots using the data in `icgcData.txt`.

## Banana Venn
The `bananaPlot.R` script can be used to generate the UpSetR version of the "Banana Venn" plot (D’Hont,  A.,  et  al.,  2012).

## SNP Caller
The `SNPcallerPlot.R` script can be used to generate the UpSetR version of the SNP caller Euler diagram presented in Xu et al., 2012.
