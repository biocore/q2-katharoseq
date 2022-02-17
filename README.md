# katharoseq
![](https://github.com/qiime2/q2templates/workflows/ci/badge.svg)

Katharoseq is a qiime2 plugin for analyzing samples with low biomass.

## Installation
To install, download this repository and navigate to the `q2_katharoseq` folder and run:
    pip install -e .

## Usage

    qiime katharoseq read-count-threshold \
        --p-threshold 50 \
        --i-table table.qza \
        --m-positive-control-column-column p_control \
        --m-positive-control-column-file metadata.tsv \
        --p-positive-control-value katharoseq_control \
        --p-control classic \
        --m-cell-count-column-column max_cell_count \
        --m-cell-count-column-file metadata.tsv \
        --o-visualization visualization.qzv \

 
