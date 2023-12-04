Environment and modules used for this analysis:
> conda activate conda_envs/scAmpi_scRNA ; module load eth_proxy ; export PATH=path/to/cellranger/cellranger-6.0.1:$PATH

2020-12
Code and resources required for omentum sample cohort analysis

Basic single cell analysis functionality is taken from the NEXUS/single_cell_analysis git repository

2020-12-17
- novaseq preprocessing is not possible with 36 samples, the index hopping removal takes an infeasible amount of memory

- test runs for index hopping removal in:
gfb_omentum_2020/analysis/2020-12_cohort/preprocessing/

2021-01
- Rerun pilot samples with latest tupro pipeline, for updated normalization -> preparation for cell type marker classification

- Initial single sample run also with latest tupro pipeline (similar to ovca SOP 7.8)

2021-01-21
- Prepared gene set list for more comprehensive gsva analysis, rerun gsva step for cohort samples

