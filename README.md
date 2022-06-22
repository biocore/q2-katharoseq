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

Computation assumes that the user has classified their 16S features against SILVA, and that the `FeatureTable[Frequency]` has been collapsed to the genus level. Please see the [`q2-feature-classifier`](https://docs.qiime2.org/2021.11/plugins/available/feature-classifier/classify-sklearn) for detail on how to perform taxonomy
classification, and the [`q2-taxa`](https://docs.qiime2.org/2021.11/plugins/available/taxa/collapse) plugin for information on collapsing to a taxonomic level. If you need more information on how to process your data, please refer to one of the relevant tutorials that can be found [here](https://docs.qiime2.org/2022.2/tutorials/).


## Obtaining a read count threshold

Computation of a minimum read count threshold can be performed with the
`read_count_threshold` plugin action. Test data can be found under the `example ` folder.

```
qiime katharoseq read_count_threshold \
    --i-table example/example_table_genus.qza \ # a genus level table
    --p-threshold 80 \
    --p-control classic \
    --p-positive-control-value katharoseq_control \ # name of controls in metadata
    --m-positive-control-column-file example/simple_katharoseq_metadata.tsv \ # your metadata
    --m-positive-control-column-column control_or_sample \ # sample type variable in metadata
    --m-cell-count-column-file example/simple_katharoseq_metadata.tsv \
    --m-cell-count-column-column max_cell_count \ # cell count variable in metadata
    --o-visualization result.qzv
```

