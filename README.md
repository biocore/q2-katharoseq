# KatharoSeq

![](https://github.com/qiime2/q2templates/workflows/ci/badge.svg) [![DOI](https://zenodo.org/badge/459343122.svg)](https://zenodo.org/badge/latestdoi/459343122)


An implementation of the KatharoSeq protocol, originally defined in [Minich et al 2018 mSystems](https://journals.asm.org/doi/10.1128/mSystems.00218-17) and in [Minich et al 2022](https://www.biorxiv.org/content/10.1101/2022.03.07.483203v1)

## Installation

Installation assumes a working QIIME 2 environment with a minimum version of 2021.8. Details on installing QIIME 2 can be found [here](https://docs.qiime2.org/2021.11/install/).

```
git clone https://github.com/biocore/q2-katharoseq.git
cd q2-katharoseq
pip install -e .
```

## Use

Computation assumes that the user has classified their 16S features against SILVA, and that the `FeatureTable[Frequency]` has been collapsed to the genus level. Please see the [`q2-feature-classifier`](https://docs.qiime2.org/2022.2/plugins/available/feature-classifier/classify-sklearn/) for detail on how to perform taxonomy classification, and the [`q2-taxa`](https://docs.qiime2.org/2022.2/plugins/available/taxa/collapse/) plugin for information on collapsing to a taxonomic level. If you need more information on how to process your data, please refer to one of the relevant tutorials that can be found [here](https://docs.qiime2.org/2022.2/tutorials/). For these examples, data from the Fish Microbiome Project (FMP): [Fish microbiomes 101: disentangling the rules governing marine fish mucosal microbiomes across 101 species](https://www.biorxiv.org/content/10.1101/2022.03.07.483203v1) paper will be used, and can be found in the `example` folder.

For a less stringent filter, an appropriate value for `--p-threshold` would be around 50. For a more strict filter, use a value like 90, as shown below in the example.


## Read Count Threshold

In order to obtain a read count threshold, computation of a minimum read count threshold can be performed with the
`read-count-threshold` plugin action. Test data can be found under the `example` folder.

```
qiime katharoseq read-count-threshold \
    --i-table example/fmp_collapsed_table.qza \
    --m-positive-control-column-file example/fmp_metadata.tsv \
    --m-positive-control-column-column control_rct \
    --m-cell-count-column-file example/fmp_metadata.tsv \
    --m-cell-count-column-column control_cell_into_extraction \
    --p-positive-control-value control \
    --p-control classic \
    --p-threshold 90 \
    --o-visualization result_fmp_example.qzv
```

## Estimating Biomass

 Estimate the biomass of samples using KatharoSeq controls. After obtaining a read count threshold using the action above, use the same metadata and collapsed table as input. The `--p-pcr-template-vol` and `--p-dna-template-vol` values are numeric values that should come from your experimental procedures.
 
```
qiime katharoseq estimating-biomass \
    --i-table example/fmp_collapsed_table.qza \
    --m-control-cell-extraction-file example/fmp_metadata.tsv \
    --m-control-cell-extraction-column control_cell_into_extraction \
    --p-min-total-reads 1315 \
    --p-positive-control-value control \
    --m-positive-control-column-file example/fmp_metadata.tsv \
    --m-positive-control-column-column control_rct \
    --p-pcr-template-vol 5 \
    --p-dna-extract-vol 60 \
    --m-extraction-mass-g-column extraction_mass_g \
    --m-extraction-mass-g-file example/fmp_metadata.tsv \
    --o-estimated-biomass estimated_biomass_fmp_rct
```

## Biomass Plot

Finally in order to visualize the results from `estimating-biomass`, run `biomass-plot`.

```
qiime katharoseq biomass-plot \
    --i-table example/fmp_collapsed_table.qza \
    --m-control-cell-extraction-file example/fmp_metadata.tsv \
    --m-control-cell-extraction-column control_cell_into_extraction \
    --p-min-total-reads 1315 \
    --p-positive-control-value control \
    --m-positive-control-column-file example/fmp_metadata.tsv \
    --m-positive-control-column-column control_rct \
    --o-visualization biomass_plot_fmp
```


