2021-01-21 
prepare the metadata table to merge cohort samples and visualize genes, pathways, clinical feature etc

script that will be used for cohort analysis:
git/scripts/cohort_analysis_omentum.R

directory to create metadata table:
gfb_omentum_2020/data/2021-01_metadata_table

from within the cohort analysis directory (gfb_omentum_2020/analysis/2020-12_cohort/single_samples/analysis/), get RDS files, h5 files, gsva files:

> for f in atypical_removed/*.RDS ; do filepath=$(realpath $f) ; echo $filepath >> gfb_omentum_2020/data/2021-01_metadata_table/RDS_files.txt ; done

> for f in counts_corrected/*.variable_genes.h5 ; do filepath=$(realpath $f) ; echo $filepath >> gfb_omentum_2020/data/2021-01_metadata_table/h5_variable_genes_files.txt ; done

> for f in gsva/*sce_gsva.RDS ; do filepath=$(realpath $f) ; echo $filepath >> gfb_omentum_2020/data/2021-01_metadata_table/gsva_files.txt ; done


- sort h5 files, RDS files, and gsva files in the ordering of the excel table:

> for name in $(less ../../git_gfb_omentum_2020/docs/sampleNames_2020-12-10.txt) ; do grep $name RDS_files.txt >> RDS_files_sorted.txt ; done

> for name in $(less ../../git_gfb_omentum_2020/docs/sampleNames_2020-12-10.txt) ; do grep $name h5_variable_genes_files.txt >> h5_variable_genes_files_sorted.txt ; done 

> for name in $(less ../../git_gfb_omentum_2020/docs/sampleNames_2020-12-10.txt) ; do grep $name gsva_files.txt >> gsva_files_sorted.txt ; done

- Get metadata from excel table (2020-12-10_Sample Tracking Sheet Omentum.xlsx):
Date_of_Exp>    Cohort_Type>    Sample Name>    Exp>    Patient> Type>  Diagnosis>      FIGO>   Grade>  subtype>Age>    RBC lysis>      DCR>    Viability (%)

Changed "Cohort_Type"  "Benign Omentum" to "Non-malignant Omentum"
Fixed typos, simplified descriptions, added a new column to distiguish tissue-type from region

----------
Question: B454_Tumor_BL: 
Borderline tumor (in omentum/ implant) -> what kind? Ma/Mi/Ga
-> assumed to be Majus
----------

- sanity check that alle files have been correctly combined:
> for name in $(less ../../git_gfb_omentum_2020/docs/sampleNames_2020-12-10.txt) ; do echo $name ; grep -c $name metadata_omentum_cohort.txt ; done


