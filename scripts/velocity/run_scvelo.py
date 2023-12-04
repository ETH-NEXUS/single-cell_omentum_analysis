#! /usr/bin/env python
# -*- coding: utf-8 -*-
#a='a'
"""
Implementation of scvelo for RNA velocity estimation & analysis of scRNA-seq data w/ spliced & unspliced counts data.
https://scvelo.readthedocs.io/

Original structure: Marcus Lindberg, March 2020
Matteo Carrara November 2021

Requirements:
module load new r/3.6.0
conda activate scvelo

Important remarks:
    Dependencies:
        igraph is necessary for PAGA and for velocity clustering. The correct conda package is python-igraph
	louvain is necessary for velocity clustering. The correct conda package is louvain. WARNING: there is a package called python-louvain, that does not work. This is due to the fact that the package python-louvain is imported within the velocity clustering script as -import louvain- but the package should be imported as -community-. This cannot be changed as the import happens within the velocity clustering script from scVelo. The conda package -louvain- instead, allows to call directly -import louvain- 
    Special variables:
        The analysis includes one special variable called 'mysamplename'. Adding it at script call (--myvariable mysamplename) will trigger UMAP plots divided by sample name. This is useful when working with cohorts, in order to highlight the position of samples in the UMAPs and compare it with all other results
"""

import argparse
import sys
import scvelo as scv
from os import listdir, chdir, mkdir, getcwd, path
import anndata2ri
from rpy2.robjects import r
from rpy2.robjects import pandas2ri
import pickle
from pathlib import Path
import matplotlib
import numpy as np
from difflib import SequenceMatcher
import copy
import anndata
import loompy
import igraph
import louvain

parser = argparse.ArgumentParser(description = 'Estimates RNA velocity of spliced/unspliced counts data to generate velocity plots.')
parser.add_argument('--loom_file', dest = 'ldata_path', required = True, help = 'Path to .loom file containing spliced & unspliced counts matrices. If working with a cohort, it should be a list of paths of all loom files involved, one per line.')
parser.add_argument('--rds_file', dest = 'adata_path', required = True, help = 'Path to .RDS file containing processed SCE object w/ cell type annotations & UMAP coordinates. If working with a cohort, it should be a list of paths of all SCE object files involved, one per line, with the samples in the same order provided in the loom file list.')
parser.add_argument('--out_dir', dest = 'out_dir', required = True, help = 'Path to output directory for all plots.')
parser.add_argument('--sample_name', dest = 'sample_name', required = True, help = 'Name of sample to be used as prefix for all output files.')
parser.add_argument('--cohort_sample_names', dest = 'cohort_sample_names', required = False, help = 'Path to a list of sample names to assign to the provided cohort samplesi, one per line, with the samples in the same order provided in the loom file list.')
parser.add_argument('--colour_config_final_file', dest = 'colour_config_final_file', required = False, default = None, help = 'Path to file specifying colours to be used for specific final cell types in plots. (Default = none)')
parser.add_argument('--colour_config_major_file', dest = 'colour_config_major_file', required = False, default = None, help = 'Path to file specifying colours to be used for specific major cell types in plots. (Default = none)')
parser.add_argument('--hex_map_path', dest = 'hex_map_path', required = False, help = 'Path to file specifying mapping of colour names to hexadecimal values.')
parser.add_argument('--n_neighbours', dest = 'n_n', required = False, type = int, default = 20, help = 'Number of neighbors to be used for calculating 1st & 2nd order moments, as well as for plotting. (Default: 20)')
parser.add_argument('--n_top_gene_plots', dest = 'n_top_gene_plots', required = False, type = int, default = 10, help = 'Number of top velocity genes to generate individual velocity plots for. (Default: 10)')
parser.add_argument('--myvariable', action='append', dest = 'myvariable', required = True, help = ' Variable in the SingleCellExperiment object to be used for UMAP computation. This option can be provided multiple time and the script will loop through all of them. (Default: celltype_major) WARNING: The first element will be used for differential kinetic computation and correction')
parser.add_argument('--embedding', dest = 'myembedding', required = False, default = "umap", help = 'Embedding to use to plot the umaps. It must be an emebedding that exists in the input SingleCellEperiment object. (Default: umap)')
parser.add_argument('--min_shared_counts', dest = 'min_shared_counts', required = False, default = 30, help = 'Minimum number of counts (both unspliced and spliced) required for a gene during normalization')
parser.add_argument('--n_top_genes', dest = 'n_top_genes', required = False, default = 2000, help = 'Number of genes to consider during normalization. Not to be confused with option "n_top_gene_plot"')
parser.add_argument('--n_jobs', dest = 'n_jobs', required = False, default = 1, type = int, help = 'Number of concurrent jobs for parallelized steps')
parser.add_argument('--use_existing_merged_loom', dest = 'use_existing_merged_loom', required = False, default = 1, type = int, help = 'Wether to use a pre-existing all_looms.loom file present in the output folder to skip the loom merging procedure')
parser.add_argument('--colour_config_file', dest = 'colour_config_file', required = False, default = None, help = 'Path to file specifying colours to be used')

args = parser.parse_args()
#scv.logging.print_version()

def print_break():
    print("#"*60)

###############################################################################
################################ I N P U T ####################################
###############################################################################

################################ LOAD DATA ####################################
'''
Load necessary data files:
  1. loom file w/ spliced/unspliced counts
  2. RDS file w/ cell type annotations & UMAP coordinates.
Files will be merged into one data object (AnnData h5ad format).
'''
def load_adata(adata_path, sample_name):
    print_break()
    print("Reading anndata h5ad with processed counts & metadata: ", adata_path)
    adata = scv.read(adata_path, cache = True)
    adata.obs['mysamplename'] = sample_name
    return adata

def load_ldata_cohort(ldata_path, out_path, use_existing_merged_loom):
    with open(ldata_path) as f:
        lines = [line.rstrip() for line in f]
    if path.exists(out_path + "all_looms.loom") & use_existing_merged_loom == 1:
        print("Found exisiting merged loom file from previous run. Skipping merging. Please set --use_existing_merged_loom 0 or delete the file all_looms.loom in the output directory to avoid this")
    else:
        print("Reading and merging loom files for all cohort samples, this may take some time and a lot of resources")
        loompy.combine(lines, out_path + "all_looms.loom")
    full_loom = scv.read_loom(out_path + "all_looms.loom")
    return full_loom

def load_ldata_single_sample(ldata_path):
    print("Reading loom file: ", ldata_path)
    ldata = scv.read(ldata_path)
    return ldata

def merge_data(adata, ldata):
    barcodes_sc_last = [x[-1] for x in list(ldata.obs_names)]
    if not(any( (barcode.endswith('x') or barcode.endswith('X')) ) for barcode in barcodes_sc_last):
         return(0)
    barcodes_sc_sn = [barcode[0:barcode.rfind(':')] for barcode in list(ldata.obs_names)]
    if any(len(barcode) == 0 for barcode in barcodes_sc_sn):
         return(1)
    ldata_barcodes = set([barcode[barcode.rfind(':')+1:-1] for barcode in list(ldata.obs_names)])
    print("Number of cells in loom file (spliced & unspliced raw counts): " + str(len(ldata_barcodes)))
    adata_barcodes = set(list(adata.obs_names))
    adata_barcodes = set([barcode[barcode.rfind('-')+1:] for barcode in list(adata.obs_names)])
    print("Number of cells in processed SingleCellExperiment RDS: " + str(len(adata_barcodes)))
    matching_barcode_counts = len(ldata_barcodes.intersection(adata_barcodes))
    print(str(matching_barcode_counts) + " of " + str(len(ldata_barcodes)) + " barcodes in loom file found in processed RDS.")
    print("Corresponding to " + str("{:.2f}".format(matching_barcode_counts/len(ldata_barcodes)*100)) + "% of loom barcodes and " + str("{:.2f}".format(matching_barcode_counts/len(adata_barcodes)*100)) + "% of SingleCellExperiment RDS barcodes")
    if matching_barcode_counts == 0:
        print("No matching barcodes found. Double check input files. If data is integrated, make sure to pre-process combined loom file so that it adheres to the barcode format in the SCE RDS.")
    if matching_barcode_counts < len(adata_barcodes)/10:
        return(2)
    for barcode in ldata_barcodes.symmetric_difference(adata_barcodes):
        if barcode not in ldata_barcodes:
            print("WARNING: there are barcodes in the processed data that are not in the unprocessed data. Double check the input files.")
            break
    print_break() 
    print("Merging loom file and RDS into single file.")
    data_merged = scv.utils.merge(adata, ldata)
    print_break()
    return data_merged

########################## CALCULATE RNA VELOCITY ##############################
'''
Preprocess & process data for RNA velocity estimation (using dynamical model) w/ awareness for cell type-specific kinetics.
'''
def estimate_velocity(data, out_path, n_neighbors, min_shared_counts, n_top_genes, sample_name, myembedding, myjobs):
    print_break()
    print("Number of cells in whole sample: ", len(data.obs['barcodes']))
    #If the number of cells are limited (less than 20), bypass the user-provided number of neighbors and use n_cells-1. Send a warning if that happens
    n_n = get_n_neighbors(data, n_neighbors)
    #Print and plot the proportion of spliced and unspliced RNA
    scv.utils.show_proportions(data)
    plt = scv.pl.proportions(data)
    matplotlib_save(plt, fname=f"{out_path}/{sample_name}.proportions.pdf")
    try:
        scv.pp.filter_and_normalize(data, min_shared_counts, n_top_genes, log = True)
    except KeyError:
        print("Error in filtering step. Re-attempting without normalizing gene dispersion...")
        scv.pp.filter_genes(data)
        scv.pp.normalize_per_cell(data)
        scv.pp.log1p(data)
    #Compute first and second order moments (means and uncentered variances) for deterministic and stochastic velocity estimation respectively
    scv.pp.moments(data, n_pcs=30, n_neighbors=n_n)
    #Depending on the method used and on the sample being cohort or not, the umap may be called differently
    #We search for a perfect match between all keys and myembedding. If no match is available, find a match that contains the full string of myembedding
    #Report any other partial or full matches in any case. Fail if there are no perfect matches and more than 1 partial match
    #We use iterable unpacking with the * operator
    res = [i for i in [*data.obsm] if myembedding in i]
    res2 = [i for i in [*data.obsm] if myembedding == i]
    if len(res) == 0:
        sys.exit("FATAL: no embedding matching (completely nor partially) the provided embedding\nProvided: " + myembedding + " - Available keys: " + str([*data.obsm]))
    if (len(res) > 1) and len(res2) == 0:
        sys.exit("More than 1 embedding found partially matching the selected embedding and no perfect match. Cannot decide which one to select. Please provide a more precise embedding.\nProvided: " + myembedding + " - Available keys: " + str([*data.obsm]))
    if (len(res) > 1) and len(res2) > 1:
        sys.exit("Multiple perfect matches of the selected embedding. This should never happen, please check original objects.\nProvided: " + myembedding + " - Available keys: " + str([*data.obsm]))
    if (len(res) > 1) and len(res2) == 1:
        print("One or more partial matches of the selected embedding detected alongside one perfect match. Continuing with the perfect match.\nProvided: " + myembedding + " - Available keys: " + str([*data.obsm]))
        myembedding = res2[0]
    if len(res2) == 1:
        print("One perfect match found of the selected embedding\nProvided: " + myembedding + " - Available keys: " + str([*data.obsm]))
        myembedding = res2[0]
    if len(res) == 1 and len(res2) ==0:
        print("One partial match found of the selected embedding\nProvided: " + myembedding + " - Available keys: " + str([*data.obsm]))
        myembedding = res[0]
    #create a dedicated X_umap for the analysis where the embedding is stored
    data.obsm["X_umap"] = data.obsm[myembedding]
    myembedding="X_umap"
    myempty = np.empty((len(data.obsm['X_umap']), 2), np.float32)
    for i in range(len(myempty)):
        myempty[i][0] = np.array(data.obsm['X_umap']["umap1"][i])
        myempty[i][1] = np.array(data.obsm['X_umap']["umap2"][i])
    data.obsm["X_umap_corrected"] = myempty
    myembedding = "X_umap_corrected"
    scv.tl.recover_dynamics(data, n_jobs=myjobs)  # this step required for dynamical mode; time-intensive
    #Start by computing the velocity without differential kinetics. This can be done afterwards
    scv.tl.velocity(data, mode='dynamical', diff_kinetics=False)
    scv.tl.velocity_graph(data, n_jobs=myjobs)
    scv.tl.recover_latent_time(data)
    print_break()
    return [data, myembedding]


###############################################################################
############################### O U T P U T ###################################
###############################################################################

################################## PLOTS ######################################
'''
Generate & save velocity plots: UMAP embedding stream, heatmap of top genes,
latent-time scatter plots, & individual gene plots.
Uses color config & labels from config passed by plot_scvelo fxn.
Only to be called by plot_scvelo fxn.
'''
def save_plots(sample_name, data, config, n_top_gene_plots, n_neighbors, myembedding, out_path, myvariable):
    top_genes = data.var['fit_likelihood'].sort_values(ascending=False).index
    if(n_top_gene_plots > len(top_genes)):
        print(f'WARNING: Asked to plot more genes ({n_top_genes_plots}) than the maximum number of top variable genes detected ' + len(top_genes) + ". Dropped the number of genes to plot to the total amount of top variable genes available. Please consider reducing the value of the variable n_top_genes_plots")
        max_genes = len(top_genes)
    else:
        max_genes = n_top_gene_plots
    matplotlib.use('Agg')
    n_n = get_n_neighbors(data, n_neighbors)
    if n_n < 20:
        print("Sample size too small. Skipping plots...")
        pass
    if len(top_genes) > len(set(top_genes)):
        print("WARNING: possible duplicated genes:")
        print(set([gene for gene in top_genes if top_genes.count(gene) > 1]))
    plt = scv.pl.scatter(data, basis=myembedding, color='latent_time', #size=100,
        color_map='gnuplot', perc=[2, 98], colorbar=True, rescale_color=[0, 1],
        show=False, legend_loc=config['legend_loc'],
        title=config['title'] + " latent time", figsize=(14,10))
    matplotlib_save(plt, fname=f"{out_path}/{sample_name}.latent_time.pdf")
    gene_path = out_path + "top_genes/"
    Path(gene_path).mkdir(parents=True, exist_ok=True)
    for myvar in myvariable:
        print(f'Plotting for variable {myvar}')
        data.obs[myvar]=data.obs[myvar].astype("string")
        plt = scv.pl.velocity_embedding_stream(data, basis=myembedding, color=config[myvar+'_colours'],fontsize=20,
            figsize=(14, 10), show=False, legend_loc=config['legend_loc'],
            title=config['title']+" velocity", palette=config[myvar+'_palette'], n_neighbors=n_n)
        print("Warning: due to matplotlib limitations, the velocity plot cannot be output as pdf. Using svg instead")
        matplotlib_save(plt, fname=f"{out_path}/{sample_name}.velocity.{myvar}.svg")
        plot_speed_coherence(data, myvar, out_path, sample_name, myembedding)
        plt = scv.pl.heatmap(data, var_names=top_genes[0:300], tkey='latent_time', n_convolve=100, col_color=myvar, show=False, font_scale=0.5, xticklabels=False)
        print("Warning: due to matplotlib limitations, the top driver genes heatmap can only be output as png")
        matplotlib_save(plt, fname=f"{out_path}/{sample_name}.top_driver_genes_heatmap.{myvar}.png")
        plot_velocity_clustering(data, myvar, out_path, sample_name, myembedding, config)
        print('Generating plots for top variable genes.')
        gene_counter = 1
        for gene in top_genes[0:max_genes]:
            print(f'Working on gene {gene}')
            Path(gene_path + str(gene_counter) + "-" + gene).mkdir(parents=True, exist_ok=True)
            plt = scv.pl.scatter(data, basis=myembedding, var_names=gene,color=config[myvar+'_colours'], colorbar=True, palette=config[myvar+'_palette'], size = 30, linewidth=2, legend_loc="right", legend_fontsize = 12, figsize=(13,10), show=False)
            matplotlib_save(plt, fname=f"{gene_path}/{str(gene_counter)}-{gene}/{sample_name}.{str(gene_counter)}-{gene}.phase_portrait.{myvar}.pdf")
            plt = scv.pl.scatter(data, basis=myembedding, x='latent_time', y=gene, color=config[myvar+'_colours'], colorbar=True, palette=config[myvar+'_palette'], size = 30, linewidth=2, legend_loc="right", legend_fontsize = 12, figsize=(13,10), show=False)
            matplotlib_save(plt, fname=f"{gene_path}/{str(gene_counter)}-{gene}/{sample_name}.{str(gene_counter)}-{gene}.phase_latent_time.{myvar}.pdf")
            plt = scv.pl.velocity(data, basis=myembedding, var_names=gene,color=config[myvar+'_colours'], color_map = "RdBu_r", colorbar=True, legend_loc='none', size = 15, linewidth=2, figsize=(20,15), show=False)
            matplotlib_save(plt, fname=f"{gene_path}/{str(gene_counter)}-{gene}/{sample_name}.{str(gene_counter)}-{gene}.multiplot.{myvar}.pdf")
            gene_counter = gene_counter + 1

'''
Configure plot parameters (colors, margin, labels, size, title).
Passes config to save_plots.
'''
def plot_configuration(sample_name, data, out_path, myvariable, hex_map_path = None, colour_config_final_path = None, colour_config_major_path = None):
    print("Setting up the plot configuration")
    plot_config = {  # "legend_loc":"upper right",
        "legend_loc": "right",
        "min_mass": "0",  # 0 - all trajectories; 100 - large only
        "title": sample_name.replace('_', '-'),
        "layer": ""
        }
    phenograph_clusters_colours = ["red", "green", "blue", "cyan", "yellow", "purple", "brown","chocolate", "chartreuse", "darkgoldenrod", "steelblue", "slateblue", "olivedrab", "gold", "violetred", "darkcyan","orchid", "darksalmon", "darkslategrey", "khaki", "indianred","magenta", "slategray", "darkolivegreen", "mediumaquamarine", "hotpink", "darkorange", "bisque", "darkseagreen", "dodgerblue", "deeppink", "tan", "mediumorchid"]
    for myvar in myvariable:
        plot_config[myvar+'_colours'] = myvar
        plot_config[myvar+'_palette'] = phenograph_clusters_colours
    if(colour_config_final_path): 
        cell_types = list(data.obs["celltype_final"].cat.categories)
        colour_config = read_colour_config(colour_config_path = colour_config_final_path, hex_map_path = hex_map_path)
        data.uns['celltype_final_colours'] = []    
        colour_config_map = {}
        plot_config["celltype_final_palette"] = map_type_colours(cell_types, colour_config)
    else:
        #Create and use a standard palette
        data.uns['celltype_final'] = phenograph_clusters_colours
    if(colour_config_major_path): 
        cell_types = list(data.obs["celltype_major"].cat.categories)
        colour_config = read_colour_config(colour_config_path = colour_config_major_path, hex_map_path = hex_map_path)
        data.uns['celltype_major_colours'] = []    
        colour_config_map = {}
        plot_config["celltype_major_palette"] = map_type_colours(cell_types, colour_config)
    else:
        plot_config["celltype_major_palette"] = phenograph_clusters_colours
    plot_config["phenograph_clusters_palette"] = read_colour_config(hex_map_path = hex_map_path, manual_colour_list = phenograph_clusters_colours)
    plot_config['generic_palette'] = read_colour_config(hex_map_path = hex_map_path, manual_colour_list = phenograph_clusters_colours)
    print(plot_config)
    return(plot_config)

########################## Plotting helper fxns ###############################
'''
Helper fxn to set appropriate number of neighbours if sample size is small.
Default is 20. (Can be passed as command-line arg.)
'''
def get_n_neighbors(data, n_neighbors):
    n_cells = data.n_obs
    if n_cells < 20:
        print(f"Warning: The sample size is below 20. The value of variable n_n passed at command line is bypassed and set to {n_cells -1}")
        n_n = n_cells - 1
    else:
        n_n = n_neighbors
    return n_n

'''
Helper function for calculating similarities between cell type labels.
Matches are determined once maximal similarity scores are found.
'''
def match_cell_type(type_to_match, reference_list):
    from difflib import SequenceMatcher
    def similar(a, b):
        return SequenceMatcher(None, a, b).ratio()
    cell_type = ''
    previous_sim_score = 0
    for each in reference_list:
        sim_score = similar(each, type_to_match)
        if (sim_score < previous_sim_score):
            continue
        else:
            cell_type = each
            previous_sim_score = sim_score
    return cell_type

'''
Maps color names from color config file to color hexadecimal values (required by plotting functions).
'''
def read_colour_config(hex_map_path, colour_config_path = None, manual_colour_list = None):
    colour_map = {}
    print("Colour config used: ", colour_config_path)
    print("Hex colour mapping used: ", hex_map_path)
    if colour_config_path:
        with open(colour_config_path) as f:
            for line in f:
                t = line.split()
                colour_map[t[0]] = t[1]
    elif manual_colour_list:
        for i in range(len(manual_colour_list)):
            colour_map[i+1] = manual_colour_list[i]
            colour_map[str(i+1)] = manual_colour_list[i]
    else:
        sys.exit("ERROR: No colour config file provided for function read_colour_config call")
    if('cell_type' in colour_map.keys()):
        colour_map.pop('cell_type')
    hex_map = open(hex_map_path, "r").read()
    for line in hex_map.split('\n')[:-1]:
        name = line.split(',')[0]
        colour = line.split(',')[1]
        for k, v in colour_map.items():
            if v.lower() == name.lower() or v.lower() == name.replace(" ", "").lower():
                colour_map[k] = colour
            elif v.lower() == name.lower()[:-1] or v.lower() == name.replace(" ", "").lower()[:-1]:
                colour_map[k] = colour
    colour_map['unknown'] = "#000000"
    return colour_map

'''
Assigns colors to cell type labels based on hexadecimal value mapping.
Similar cell types (but w/ arbitrarily different annotation) will be assigned the same color.
'''
def map_type_colours(types_to_map, hex_colour_config):
    colour_map = {}
    print(types_to_map)
    for cell_type in types_to_map:
        matched = match_cell_type(cell_type, list(hex_colour_config.keys()))
        colour_map[cell_type] = hex_colour_config[matched]
    return colour_map


############################### MATPLOTLIB SAVE ##############################
def matplotlib_save(plt, fname):
    matplotlib.pyplot.savefig(fname, bbox_inches='tight')
    print(f"Saving figure to {fname}")
    matplotlib.pyplot.close()


############################### GENE LISTS ###################################
'''
Save top 300 most-likely velocity driver genes to text file.
Each gene is separated by newline.
'''
def save_top_genes_list(out_path, sample_name, data, save=False):
    top_likelihood_genes = data.var['fit_likelihood'].sort_values(ascending=False)[:300].to_frame()
    top_likelihood_genes = top_likelihood_genes.rename_axis("gene").reset_index()
    top_genes = top_likelihood_genes['gene'].values
    print(f"{sample_name} - top velocity genes: ", list(top_genes)[:11])
    if save:
        top_likelihood_genes.to_csv(f'{out_path}/{sample_name}.top_likelihood-ranked_genes.tsv', sep = '\t')
        print(f"{sample_name} velocity gene list saved.")


############################### GENE LISTS BY VARIABLE ######################
def save_top_genes_list_byvar(out_path, sample_name, data, myvar, save=False):
    scv.tl.rank_dynamical_genes(data, groupby=myvar)
    df = scv.get_df(data, 'rank_dynamical_genes/names')
    if save:
        df.to_csv(f"{out_path}/{sample_name}.top_likelyhood-ranked_genes_by_{myvar}.tsv", index = False, sep = '\t')
        print(f"{sample_name} velocity gene list by {myvar} saved.")


############################### SPEED AND COHERENCE #########################

def plot_speed_coherence(data, myvar, out_path, sample_name, myembedding):
    print('Plotting speed of rate of differentiation and coherence of the vector field')
    scv.tl.velocity_confidence(data)
    keys = 'velocity_length', 'velocity_confidence'
    plt = scv.pl.scatter(data, basis=myembedding, c=keys,cmap='coolwarm', perc=[2,98], size = 30, linewidth=1, legend_loc="right", legend_fontsize = 10, figsize=(13,10), show=False)
    matplotlib_save(plt, fname=f"{out_path}/{sample_name}.velocity-confidence.pdf")
    print (f'Saving the speed and coherence table divided by {myvar}')
    df = data.obs.groupby(myvar)[keys].mean().T
    df.style.background_gradient(cmap='coolwarm',axis=1)
    print(f'Saving the velocity length confidence table divided by {myvar}')
    df.to_csv(f'{out_path}/{sample_name}_velocity_length_confidence_{myvar}.tsv', sep = '\t')
    return(data)


############################# DIFFERENTIAL KINETICS #################
def diff_kinetics(data, sample_name, myvar, n_top_gene_plots, out_path, config, myjobs, myembedding, n_neighbors):
    print(f'Starting computation of differential kinetics. Cells in UMAPS will be colored by {myvar}')
    n_n = get_n_neighbors(data, n_neighbors)
    if n_top_gene_plots > 30:
        print(f'WARNING: n_top_gene_plots variable set too high, the resulting plot may be too big or difficult to read. Please consider lowering the value to 30 or lower')
    kwargs = dict(linewidth=2, add_linfit=True, frameon=False)
    #Extract the top 300 most varible genes and compute the differential kinetic test on those. 300 is chosen to reduce computational time and avoid to introduce noise by including genes with lower likelihoods
    top_genes = data.var['fit_likelihood'].sort_values(ascending=False).index[:300]
    scv.tl.differential_kinetic_test(data, var_names=top_genes, groupby=myvar)
    plt = scv.pl.scatter(data, basis=top_genes[:n_top_gene_plots], ncols=5, add_outline='fit_diff_kinetics', color=config[myvar+'_colours'], colorbar=True, palette=config[myvar+'_palette'], size = 30, legend_loc="right", legend_fontsize = 12, figsize=(13,10), show=False, fontsize=20, **kwargs)
    matplotlib_save(plt, fname=f"{out_path}/{sample_name}.top_genes_diff_kinetics.{myvar}.pdf")
    print("Recomputing velocity based on the differential kinetics")
    scv.tl.velocity(data, diff_kinetics=True)
    scv.tl.velocity_graph(data, n_jobs=myjobs)
    plt = scv.pl.velocity_embedding_stream(data, basis=myembedding, color=config[myvar+'_colours'],fontsize=20,
            figsize=(14, 10), show=False, legend_loc=config['legend_loc'],
            title=config['title']+" velocity", palette=config[myvar+'_palette'], n_neighbors=n_n)
    matplotlib_save(plt, fname=f"{out_path}/{sample_name}.velocity_kinetic_corrected.{myvar}.svg")
    return data


############################### VELOCITY CLUSTERING #########################
def plot_velocity_clustering(data, myvar, out_path, sample_name, myembedding, config):
    print(f'Clustering velocities and plotting the clusters against the values of {myvar}')
    ##Transform the variable in categorical to avoid clashes in case of categories with numbers as names
    data.obs[myvar] = data.obs[myvar].astype('category')
    scv.tl.velocity_clusters(data, match_with=myvar)
    mycolours = copy.deepcopy(config['generic_palette'])
    j=1
    for i in set(data.obs['velocity_clusters']):
        mycolours[i] = mycolours[j]
        del mycolours[j]
        j=j+1
    plt = scv.pl.scatter(data, basis=myembedding, color='velocity_clusters', size=30, legend_loc="right", legend_fontsize=10, figsize=(13,10), show=False, palette=mycolours)
    matplotlib_save(plt, fname=f"{out_path}/{sample_name}.velocity_clusters.{myvar}.pdf")
    return(data)

############################### RUN ALL ######################################

'''
Run main scvelo preprocessing, processing, & plotting workflow.
Velocities are estimated for whole sample & by cell type.
'''
def run_scvelo(ldata_path, adata_path, out_path, sample_name, colour_config_final_path, colour_config_major_path, hex_map_path, n_top_gene_plots, myvariable, n_neighbors, min_shared_counts, n_top_genes, myembedding, myjobs, cohort_sample_names, use_existing_merged_loom):
    print("Running the velocity analysis module!")
    out_path = out_path + "/"
    print("Output dir: ", out_path)
    print("Testing if the provided data are for single samples or a cohort")
    try:
        with open(ldata_path) as f:
            ldata = f.readlines()
        ldata = [ name.rstrip() for name in ldata ]
        print("Detected cohort data. Proceeding accordingly...")
        ldata = load_ldata_cohort(ldata_path, out_path, use_existing_merged_loom)
    except (UnicodeDecodeError) as error:
        print("Detected single sample. Proceeding accordingly...")
        ldata = load_data_single_sample(ldata_path)
    adata = load_adata(adata_path, sample_name)
    data = merge_data(adata,ldata)
    print("Created the object necessary for the velocity analysis")
    del adata
    del ldata
    
    print(data)
    output = estimate_velocity(data, out_path, n_neighbors, min_shared_counts, n_top_genes, sample_name, myembedding, myjobs)
    myembedding = output[1]
    output = output[0]
    print("Number of cells in whole sample: ", len(output.obs['barcodes']))
   
    print_break()
    print("Current sample: ", sample_name)

    plot_path = out_path + "plots/"
    Path(plot_path).mkdir(parents=True, exist_ok=True)
    print(f"Saving plots (whole sample) to {out_path}/plots/...")
    plot_config = plot_configuration(sample_name, output, plot_path, myvariable, hex_map_path, colour_config_final_path, colour_config_major_path)
    save_plots(sample_name, output, plot_config, n_top_gene_plots, n_neighbors, myembedding, plot_path, myvariable)
    print("Plots saved.")

    print_break()
    top_genes_path = plot_path + "top_genes/"
    Path(top_genes_path).mkdir(parents=True, exist_ok=True)
    save_top_genes_list(top_genes_path, sample_name, output, save=True)

    diff_kin_path = plot_path + "differential_kinectics/"
    Path(diff_kin_path).mkdir(parents=True, exist_ok=True)
    for myvar in myvariable:
        print(f"Saving differential kinetics for variable {myvar}...")
        save_top_genes_list_byvar(top_genes_path, sample_name, output, myvar, save=True)
        diff_kinetics(output, sample_name, myvar[0:], n_top_gene_plots, diff_kin_path, plot_config, myjobs, myembedding, n_neighbors)

    writeLines("Run completed successfully", out_path + "/" + sample_name + ".velocity_success.txt")
    print_break()
    print(f"Completed velocity analysis. Output is available in {out_path}")
    
##############################################################################################

run_scvelo(args.ldata_path, args.adata_path, args.out_dir, args.sample_name, args.colour_config_final_file, args.colour_config_major_file, args.hex_map_path, args.n_top_gene_plots, args.myvariable, args.n_n, args.min_shared_counts, args.n_top_genes, args.myembedding, args.n_jobs, args.cohort_sample_names, args.use_existing_merged_loom)
