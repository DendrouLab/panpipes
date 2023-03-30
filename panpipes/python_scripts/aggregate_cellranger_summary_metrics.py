#!/bin/py
# aggregate_cellranger.py

# takes cellranger multi outputs and cellranger count outputs, and makes them consistent for plotting.


import pandas as pd
import glob
import os
import re
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--pipe_df',
                    default='',
                    help='')
parser.add_argument('--output_file',
                    default='',
                    help='')
parser.add_argument('--figdir',
                    default='./figures',
                    help='')
parser.add_argument('--cellranger_column_conversion_df',
                    default='',
                    help='')
parser.set_defaults(verbose=True)
args = parser.parse_args()

L.info(args)

def get_metrics_summary_path(path,sample_id=None):
    """ infers the path to metrics_summary.csv based on what cellranger algorithm is used 
    and checks the file exists

    Args:
        path (str):path to the folder called 'outs'
        sample_id (str): Required if path is to cellranger multi outputs, defaults to None
    """    
    # subset path to only go up to 'outs'
    if 'outs' not in path:
        print('the outs folder must be included in the path')
    else:
        # split to make sure that outs is at the end of the path
        path = path.split("outs")[0] + "outs"
    outpath=None
    # use the path to cellranger count or vdj outputs as default
    if os.path.exists(os.path.join(path, 'metrics_summary.csv') ):
        outpath = os.path.join(path, 'metrics_summary.csv') 
    elif sample_id is not None and os.path.exists(os.path.join(path, 'per_sample_outs', sample_id, 'metrics_summary.csv') ):
        outpath = os.path.join(path, 'per_sample_outs', sample_id, 'metrics_summary.csv')
    elif sample_id is None and os.path.exists(os.path.join(path, 'per_sample_outs')):
        print('input folder appears to be from cellranger multi but no sample_id is given')
    else:
        # use the alternative path from cellranger_multi outputs
        print('path not found')
    return outpath


def detect_cellranger_algorithm(pth):
    if 'per_sample_outs' in pth:
        filetype='multi'
    else:
        filetype='count'
    return filetype


def get_all_unique_paths(pipe_df):
    """subsets all columns of pd.Dataframe which are suffixed 'path' and stacks the result

    Args:
        pipe_df (pd.DataFrame): dataframne that contains at least one columns suffixed 'path'

    Returns:
        pd.DataFrame: stacked dataframe of paths
    """    
    pipe_df = pipe_df.set_index('sample_id')
    pipe_df = pipe_df.loc[:,pipe_df.columns.str.endswith('path')]
    all_paths = pipe_df.stack().reset_index()
    all_paths.columns = ['sample_id', 'path_type', 'path']
    all_paths = all_paths.drop_duplicates()
    # rename path_type to match cellranger terminolgy
    recode_dict = {'gex_path': "Gene Expression", 
                   'adt_path': "Antibody Capture", 
                   'tcr_path': 'VDJ T', 
                   'bcr_path': 'VDJ B'}
    all_paths['path_type'] = all_paths['path_type'].replace(recode_dict)
    return all_paths


def parse_10x_cellranger_multi(path_df,path_col='metrics_summary_path'):
    """takes the metrics_summary.csv from cellranger multi for each sample 
    and concantenates into on pd.DataFrame (long format)

    Args:
        path_df (pd.DataFrame): pandas dataframe with two columns sample_id and and 
                            path which contains the path to the outs folder.
        path_col (str): colname for the columns containing the paths

    Returns:
        pd.DataFrame: 7 columns ['sample_id', 'category', 'library_type', 'grouped_by', 
        'group_name', metric_name', 'metric_value']
    """   
    # subset to unique rows to save loading in duplicate data
    path_df = path_df[['sample_id', path_col]].drop_duplicates()
    # read and concat
    msums = pd.concat([pd.read_csv(f) for f in path_df[path_col]], 
                    keys=path_df['sample_id'], 
                    names=['sample_id']).reset_index().drop(columns='level_1')
    # remove spaces and capitils from columns names
    msums.columns = [re.sub(' ', '_', x.lower()) for x in msums.columns]
    # convert percentages from string to numeric
    msums['metric_value'] = [re.sub(",|%", "", x) for x in msums['metric_value']]
    msums['metric_value'] = msums['metric_value'].astype(float)
    # sort outputs
    msums = msums.sort_values(['metric_name', 'sample_id'])
    return msums


def parse_10x_cellranger_count(path_df, convert_df,  path_col='metrics_summary_path'):
    """takes the metrics_summary.csv from cellranger count output for each sample 
    and concantenates into on pd.DataFrame (long format)

    Args:
        path_df (pd.DataFrame): pandas dataframe with three columns sample_id, path 
                                path which contains the path to the outs folder.and 
                                path_type which contains the type of library (one of ['Gene Expression', 'Antibody Capture', 'VDJ T', 'VDJ B'])
        path_col (str): colname for the columns containing the paths

    Returns:
        pd.DataFrame: 7 columns ['sample_id', 'category', 'library_type', 'grouped_by', 
        'group_name', metric_name', 'metric_value']
    """   
    # read and concat
    msums = pd.concat([pd.read_csv(f) for f in path_df[path_col]], 
                    keys=[(x, y) for x, y in zip(path_df['sample_id'], path_df['path_type'])], 
                    names=['sample_id','path_type'])
    msums = msums.stack().reset_index().drop(columns='level_2')
    # make the columns match the conversion df
    msums.columns = ['sample_id', 'library_type', 'count_metric_name', 'metric_value']
    # update the columns names
    msums = msums.merge(convert_df).drop(columns=['count_metric_name'])
    msums['group_name'] = msums['sample_id']
    msums = msums[['sample_id','category', 'library_type', 'grouped_by', 'group_name', 'metric_name', 'metric_value']]
    # convert percentages from string to numeric
    msums['metric_value'] = [re.sub(",|%", "", str(x)) for x in msums['metric_value']]
    msums['metric_value'] = msums['metric_value'].astype(float)
    msums = msums.sort_values(['library_type', 'metric_name', 'sample_id'])
    return msums
    

convert_df = pd.read_csv(args.cellranger_column_conversion_df, sep='\t')   

# get the paths
pipe_df = pd.read_csv(args.pipe_df, sep='\t')
L.info('finding all the cellranger paths')
all_paths_df = get_all_unique_paths(pipe_df)
all_paths_df['metrics_summary_path'] = all_paths_df.apply(lambda x: get_metrics_summary_path(path=x.path, sample_id=x.sample_id), axis=1)
# all_paths_df_uniq = all_paths_df.drop(columns=['path', 'path_type']).drop_duplicates().reset_index(drop=True)
all_paths_df['cellranger_type'] = [detect_cellranger_algorithm(x) for x in all_paths_df['metrics_summary_path']]
L.info('cellranger metrics_summary files found')
L.info(all_paths_df)

# parse and concatenate the tenx x metrics summary files into long format
tenx_metrics = []
if any(all_paths_df['cellranger_type']=='multi'):
    tenx_metrics.append(parse_10x_cellranger_multi(all_paths_df[all_paths_df['cellranger_type']=='multi']))

if any(all_paths_df['cellranger_type']=='count'):
    tenx_metrics.append(parse_10x_cellranger_count(all_paths_df[all_paths_df['cellranger_type']=='count'], convert_df))



tenx_metrics_full = pd.concat(tenx_metrics)
tenx_metrics_full = tenx_metrics_full.sort_values(['library_type', 'metric_name'])
tenx_metrics_full.to_csv(args.output_file, index=False)
L.info('done')


# tenx_metrics['metric_value'] = [re.sub(",|%", "", x) for x in tenx_metrics['metric_value']]
# tenx_metrics['metric_value'] = tenx_metrics['metric_value'].astype(float)
# split by library_type and metric_name 
for idx, row in tenx_metrics_full[['library_type','metric_name']].drop_duplicates().iterrows():
    mn = row['metric_name']
    lt = row['library_type']
    plt_df = tenx_metrics_full[(tenx_metrics_full['library_type'] == lt) & (tenx_metrics_full['metric_name'] == mn)]  
    if len(plt_df.category.unique()) > 1:
        # if multiple categories i.e. Cells and Libraries, these are likely duplicate rows
        plt_df = plt_df[['sample_id', 'metric_name', 'library_type' , 'metric_value']].drop_duplicates()
    if len(plt_df['sample_id']) == len(plt_df['sample_id'].unique()):
        fig = sns.barplot(plt_df, x='sample_id', y='metric_value', color='grey')
        fig.set_xticklabels(fig.get_xticklabels(), rotation=90)
        fig.set_title(lt + ':' + mn)
        plt.savefig(os.path.join(args.figdir, lt + '-' + mn + '.png'), bbox_inches='tight')
    else:
        # do a boxplot instead
        fig = sns.boxplot(plt_df, x='sample_id', y='metric_value', color='grey')
        fig.set_xticklabels(fig.get_xticklabels(), rotation=90)
        fig.set_title(lt + ':' + mn)
        plt.savefig(os.path.join(args.figdir, lt + '-' + mn + '.png'), bbox_inches='tight')
    plt.clf()


## gene expression  extra plots -----------------------------

gex_tenx_metrics = tenx_metrics_full[tenx_metrics_full['library_type'] == 'Gene Expression']

# sequencing_saturaion plot scatter plot
plot_metrics = ['Sequencing saturation', 'Mean reads per cell', 'Estimated number of cells', 'Number of reads', 'Median UMI counts per cell']
plt_df = gex_tenx_metrics[gex_tenx_metrics.metric_name.isin(plot_metrics)]
plt_df = plt_df[['sample_id', 'metric_name', 'metric_value']]
plt_tab = plt_df.pivot_table(index='sample_id', columns='metric_name', values='metric_value', aggfunc=sum)
# sns.scatterplot(data=plt_tab, 
#                 x='Sequencing saturation', 
#                 y='Number of reads', size='Mean reads per cell',sizes=(20, 200))


f, ax = plt.subplots()
points = ax.scatter(x=plt_tab['Sequencing saturation'], 
            y=plt_tab['Number of reads'], 
            c=plt_tab['Estimated number of cells'])
# define labels
plt.xlabel('Sequencing saturation')
plt.ylabel('Number of reads')
# define color bar
cbar = f.colorbar(points)
cbar.ax.get_yaxis().labelpad = 15
cbar.ax.set_ylabel('Estimated # of cells', rotation=270)
# put it on a grid
plt.grid(axis='both', color='0.95')
# save
plt.savefig(os.path.join(args.figdir,'10x_sequencing_saturation_summary.png'))
plt.clf()

# plot 2 -  cells vs UMI counts per cell

fig, ax = plt.subplots()
x=list(np.log10(plt_tab['Estimated number of cells']))
y=list(np.log10(plt_tab['Median UMI counts per cell']))
ax.scatter(x=x, 
            y=y)
# annotate each point
for i, txt in enumerate(plt_tab.index):
    ax.annotate(txt, (x[i], y[i]), size=8, ha='center')
plt.xlabel('log10(Estimated number of cells)')
plt.ylabel('log10(Median UMI counts per cell)')
# put a boundary box around the good cells
plt.vlines(x=(np.log10(500), np.log10(20000)), 
           ymin=min(1.5, min(y)),
           ymax=max(5, max(y)),
           colors='black', linestyles='dotted'
           )
plt.hlines(y=(3, 4), 
            xmin=min(1.5, min(x)),
           xmax=max(5, max(x)),
           colors='black', linestyles='dotted'
)

# put it on a grid
plt.grid(axis='both', color='0.95')
# save
plt.savefig(os.path.join(args.figdir,'10x_cells_by_UMIs.png'))
plt.clf()