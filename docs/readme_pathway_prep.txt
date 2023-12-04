2020-01 Reformat selected genes and pathways provided by Ulli and Lilly
-> match pathways to msigDB
-> create gene sets in correct format from manual lists
-> correct naming, typos etc

2021-01-18
selected genes and pathway prep
- Download all gene sets from MsigDB (http://www.gsea-msigdb.org/gsea/downloads.jsp):
msigdb.v7.2.symbols.gmt

- from Ullis and Lilis genes and pathways, create a list of pathway names, filename: list_of_selected_pathways_from_excel.txt
https://drive.google.com/drive/u/1/folders/1H0uPJbM6STaTAuc4ij_qNWUK-4Vkqyh_

All data in directory:
gfb_omentum_2020/data/2020-12_genes_pathways/

- convert manual gene sets from Ulli into gmt format: (directory: Manual_Pathways)

> for file in *.txt ; do filename=$(basename $file .txt); echo -n $filename >> geneset_${filename}.gmt ; echo -ne '\t> ' >> geneset_${filename}.gmt ; echo -n $filename >> geneset_${filename}.gmt ;  for gene in $(cut -f 1 $file) ; do echo -ne '\t'$gene >> geneset_${filename}.gmt ; done ; done

- Macrophage subtypes were not formatted as requested -> manual distincton into the different pathways
-> no gene sets possible for therapy-resistance and tissue-residence marker (only 3 and 1 gene, respectively) -> have been included in the selected genes list (as Macro_marker)
-> CXCR4, MRC1, TIE2, TIMD4

Pathway check:
- compare list of pathwyas mentioned in the excel sheet to pathways from MSigDB, ovarian pathways, and manual pathways

Issues:
1: manual signature for Invasion_signature missing in pathway list -> removed
2: to allow a matching of gene set names, _UM has been removed from pathway excel list. Further, harmonized names between gmt file and list of pathway names (upper case, lower case + general description inconsistencies)
3:no gene set for "Mesenchymal_Stem_Cells1_UM" and "Mesenchymal_Stem_Cells2_Adipocyte_UM" provided -> removed from name list

- Combine msigb gene sets, omentum gene set, manual gene sets + remove duplicated pathways (HALLMARK gene sets and stem cell set taken out of the omentum gene set list)
NOTE: HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY is not included anymore in msigdb 7.2 -> used the omentum gene set version

Check if all pathways are found:

> for f in $(less list_of_selected_pathways_from_excel_uniq.txt) ; do echo $f >> pathway_check_tab.txt ; grep -c "$f       " combined_msigdb_omentum_manual.gmt >> pathway_check_tab.txt ; done
> grep -B 1 "^0" pathway_check_tab.txt

Not found: 9 pathways:

BIOCARTA_GLYCOLYSIS_PATHWAY
0
--
INDUCTION_OF_APOPTOSIS_BY_EXTRACELLULAR_SIGNALS
0
INDUCTION_OF_APOPTOSIS_BY_INTRACELLULAR_SIGNALS
0
--
REACTOME_EXTRINSIC_PATHWAY_FOR_APOPTOSIS
0
--
REACTOME_NEGATIVE_REGULATION_OF_MET_ACTI_ACTIVITY
0
--
REACTOME_RAF_MAP_KINASE_CASCADE
0
--
REACTOME_REGULATION_OF_APOPTOSIS
0
--
REACTOME_TNF_RECEPTOR_SUPERFAMILY_TNFSF_
0
-> changed to "REACTOME_TNF_RECEPTOR_SUPERFAMILY_TNFSF_MEMBERS_MEDIATING_NON_CANONICAL_NF_KB_PATHWAY"
--
REACTOME_APOPTOSIS
0
-> changed to "REACTOME_APOPTOSIS_INDUCED_DNA_FRAGMENTATION" 
-> is already included as a separate pathway, so REACTOME_APOPTOSIS has been excluded

Checked pathways in MSigDB, 7/9 refer to archived pathways, 2 have changed names

IN TOTAL:
8 pathways excluded, 1 with changed name

pathways with multiple hits:
CELL_MIGRATION
50 -> change to "WU_CELL_MIGRATION"
CELL_CELL_ADHESION
-> changed to  GO_CELL_CELL_ADHESION
CYTOKINE_BIOSYNTHETIC_PROCESS
2 -> include both and change name: "GO_NEGATIVE_REGULATION_OF_CYTOKINE_BIOSYNTHETIC_PROCESS" and "GO_POSITIVE_REGULATION_OF_CYTOKINE_BIOSYNTHETIC_PROCESS"
FATTY_ACID_OXIDATION
3 -> change to "GO_REGULATION_OF_FATTY_ACID_OXIDATION"
MAPK_PATHWAY
9 -> take gene set which is only called "MAPK_PATHWAY"

get relevant pathways from combined pahtway file (afterwards manually remove duplicated MAPK pathways):
> for f in $(less list_of_selected_pathways_from_excel_uniq.txt) ; do grep "$f     " combined_msigdb_omentum_manual.gmt >> cohort_omentum_genesets.gmt ; done


Dictionary for pathway to group assignment:
Some pathways had missing groups, added according to most likely group assignment (neighbors):
GO_RESPONSE_TO_GROWTH_FACTOR	Cell Response I
GO_RESPONSE_TO_GROWTH_HORMONE	Cell Response I

Some pathways were removed from dictionary because no gene sets are available (also the archived ones mentioned above):
Mesenchymal_Stem_Cells1_UM>     CellTypes_Mesenchymal
Mesenchymal_Stem_Cells2_Adipocyte_UM>   CellTypes_Mesenchymal
TherapyResistant_Macrophage_Signature_UM>       CellTypes_Macrophages


Dictionaries and grouping:
-> replace whitespace by underscore in group names

Check group sizes:
> for f in $(cut -f 2 ../../git_gfb_omentum_2020/omentum_required/dictionary_pathway_groups.txt | sort | uniq) ; do echo $f ; grep -c -w $f ../../git_gfb_omentum_2020/omentum_required/dictionary_pathway_groups.txt ; done
Adhesion
37
Antigen_Presentation
14
Cancer
16
Cancer_Ovarian
8
Cell_communication_Chemokine
28
Cell_communication_Cytokine
29
Cell_Death
5
Cell_Fusion
9
Cell_Processes
4
CellProcesses_Vesicle
26
Cell_Response_I
44
CellResponse_Inflammation
21
CellTypes_Fibroblasts
13
CellTypes_Granulocytes
1
CellTypes_Macrophages
11
CellTypes_Mesenchymal
1
CellTypes_Mesothelial
2
CellTypes_Pericytes
1
Development
5
Development_Adipocytes
12
Development_Epithelial
20
Development_Genesis
15
Development_StemCells
9
Epigenetic_regulation
5
Glycosylation
99
ligand_meso
4
Metabolism
23
Metabolism_Insulin
14
Motility
14
Proliferation
33
Protein_Folding
1
receptor_meso
16
Signaling
52
Tissue_Adaptation
23
transcription_factor_Fibr
1
transcription_factor_meso
9
VitaminD
12


- group everything below 3 group members:
CellTypes_Other = CellTypes_Granulocytes
CellTypes_Mesenchymal
CellTypes_Mesothelial
CellTypes_Pericytes
-> also add Protein_Folding

- combine transcription_factor_Fibr and transcription_factor_meso to "transcription_factor"

- split Cell_Response into 22 Cell_Response_I and 22 Cell_Response_II
- split Signaling into 26 Signaling_1 and 26 Signaling_2
- split Glycosylation into 33 Glycosylation_1 and 33 Glycosylation_2 and 33 Glycosylation_3 



Selected genes:
> for f in $(cut -f 2 ../../git_gfb_omentum_2020/omentum_required/cohort_selected_genes_omentum.txt | sort | uniq) ; do echo $f ; grep -c -w $f ../../git_gfb_omentum_2020/omentum_required/cohort_selected_genes_omentum.txt ; done
anti_inflammatory
12
Basophil_defining
5
Bcell_marker
39
cell_status
6
ColdTumor_ExclusionProgram_Short
12
cytokines
17
DendriticCell_cDC_marker
10
DendriticCell_defining
6
DendriticCell_marker
11
DendriticCell_pDC_marker
23
drug_target
38
em_transition
18
endothelial_markers
17
Eosinophil_defining
2
Epithelial_marker
9
exploratory
106
exploratory_tumorImmune
77
fibroblast_markers
27
glyco
17
Granulocyte_marker
16
HotTumor_ImmuneResponsive_Signature
32
immune_cell_type_markers
56
Immune_Regulatory
34
immunoregulatory
34
integrin_ecm
52
Macro_marker
107
Macrophage_ColdTumor_ImmuneResistant_Signature
37
MastCell_defining
2
MastCell_marker
6
MDSC_marker
15
Mesenchymal_Stem_Cell2_Adipose
6
Mesenchymal_Stem_Cells1
39
MesothelialCell_marker
53
mesothelial_FACS
4
mesothelial_Si
14
microenvironment
9
Monocyte_defining
9
Myeloid_defining
2
Monocyte_marker
18
Neutrophil_defining
4
Nkcell_marker
47
Ovarian_CSC_marker
34
Pericyte_defining
4
Pericyte_marker
8
potential_fibroblast_ligand_top100
29
potential_fibroblast_receptor_top100
20
potential_fibroblast_TF_top50
4
potential_mesothelial_ligand_top100
17
potential_mesothelial_receptor_top100
16
potential_mesothelial_TF_top50
9
pro_inflammatory
23
relevant_pathways
23
stem_cell_markers
76
stroma_cell_markers
21
---
Tcell_infiltration_Macrophage_Marker
1
Tcell_infiltration_Tcell_Marker
1
-> group into Tcell_infiltration
---
Tcell_marker
29
Th17_profile
8
Th1_profile
12
Th2_profile
8
TIME_regulation
3
Treg_profile
7
Tumor_Immunesuppression
37
Tumor_invasion
2
tumor_marker
26
