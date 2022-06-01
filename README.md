# KatharoSeq

An implementation of the KatharoSeq protocol, originally defined in [Minich et al 2018 mSystems](https://journals.asm.org/doi/10.1128/mSystems.00218-17).

## Installation

Installation assumes a working QIIME 2 environment with a minimum version of 2021.8. Details on installing QIIME 2 can be found [here](https://docs.qiime2.org/2021.11/install/).

```
git clone https://github.com/biocore/q2-katharoseq.git
cd q2-katharoseq
pip install -e .
```

## Use

Computation of a minimum read count threshold can be performed with the
`read_count_threshold` plugin action. Computation assumes that the user has
classified their 16S features against SILVA, and that the
`FeatureTable[Frequency]` has been collapsed to the genus level. Please see the
[`q2-feature-classifier`](https://docs.qiime2.org/2021.11/plugins/available/feature-classifier/classify-sklearn) for detail on how to perform taxonomy
classification, and the [`q2-taxa`](https://docs.qiime2.org/2021.11/plugins/available/taxa/collapse) plugin for information on collapsing to a taxonomic level.

```
qiime katharoseq read_count_threshold \
    --i-table a_genus_level_table.qza \
    --p-threshold 80 \
    --p-control classic \
    --p-positive-control-value name_of_controls_in_metadata \
    --m-positive-control-column-file your_metadata.tsv \
    --m-positive-control-column-column sample_type_variable_in_metadata \
    --m-cell-count-column-file your_metadata.tsv \
    --m-cell-count-column-column cell_count_variable_in_metadata \
    --o-visualization result.qzv
```
