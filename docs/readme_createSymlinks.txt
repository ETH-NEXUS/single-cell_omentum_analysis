#2020-12-14

>for name in $(less ../../../git_gfb_omentum_2020/sampleNames_2020-12-10.txt) ; do echo $name ; mkdir -p $name ;  done

Location of fastqs:
gfb_omentum_2020/data/2020-12_cohort/fastqs/

Create the symlinks: (project dir gfb_omentum_2020)
> sh git_gfb_omentum_2020/scripts/create_symlinks.sh gfb_omentum_2020/data/2020-12_cohort/fastqs/ git_gfb_omentum_2020/sampleNames_2020-12-10.txt data/2020-12_cohort/openbis/2020-12-15_download_info_full.txt
