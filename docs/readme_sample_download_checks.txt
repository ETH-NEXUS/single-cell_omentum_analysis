Check all data that was downloaded from openBIS:

> for f in 20201* ; do echo $f >> 2020-12-14_download_info.txt ; ll -h $f/original/* >> 2020-12-14_download_info.txt; done

Check sample name match and fastq number:

>for name in $(less ../../../git_gfb_omentum_2020/sampleNames_2020-12-10.txt) ; do echo $name >> 2020-12-14_name_sample_match.txt ; grep $name 2020-12-14_download_info.txt | grep "fastq.gz" >> 2020-12-14_name_sample_match.txt ; echo "files:" >> 2020-12-14_name_sample_match.txt ; grep $name 2020-12-14_download_info.txt | grep -c "fastq.gz" >> 2020-12-14_name_sample_match.txt ; done

> for f in 20201* ; do echo $f >> 2020-12-15_download_info_full.txt ; for files in $(ls $f/original/*/*) ; do filepath=$(realpath $files); echo $filepath >> 2020-12-15_download_info_full.txt; done ; done

> for name in $(less ../../../git_gfb_omentum_2020/sampleNames_2020-12-10.txt) ; do echo $name >> 2020-12-15_fullname_sample_match.txt ; grep $name 2020-12-15_download_info_full.txt | grep "fastq.gz" >> 2020-12-15_fullname_sample_match.txt ; done

> for name in $(less ../../../git_gfb_omentum_2020/sampleNames_2020-12-10.txt) ; do echo $name ; for i in {1..4} ; do echo $i ; patternR1=L00${i}_R1 ; myFileR1=$(grep $name 2020-12-15_download_info_full.txt | grep "fastq.gz" | grep $patternR1) ; echo $myFileR1 ; done ; done | less


2021-04-08 Also prepare fastq files of the initial sequencing run

directory: gfb_omentum_2020/data/2020-12_cohort/20201021_UM29_OMENTUM_SP

1. Create overview of the openbis files

>> for f in 2020102* ; do echo $f >> 2021-04-08_download_info_iniRun.txt ; for files in $(ls $f/original/*/*) ; do filepath=$(realpath $files); echo $filepath >> 2021-04-08_download_info_iniRun.txt ; done ; done

2. Check all fastq files

>> for name in $(less ../../../git_gfb_omentum_2020/docs/sampleNames_2020-12-10.txt) ; do echo $name >> 2021-04-08_fullname_sample_match_iniRun.txt ; grep $name 2021-04-08_download_info_iniRun.txt | grep "fastq.gz" >> 2021-04-08_fullname_sample_match_iniRun.txt ; done
