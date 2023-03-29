


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt



tenx_metrics = pd.read_csv("/well/cartography/users/zsj686/src/sandbox/10x_metrics_concat.csv")
tenx_metrics = tenx_metrics.sort_values(['library_type', 'metric_name'])

# tenx_metrics['metric_value'] = [re.sub(",|%", "", x) for x in tenx_metrics['metric_value']]
# tenx_metrics['metric_value'] = tenx_metrics['metric_value'].astype(float)
# split by library_type and metric_name 
for idx, row in tenx_metrics[['library_type','metric_name']].drop_duplicates().iterrows():
    mn = row['metric_name']
    lt = row['library_type']
    plt_df = tenx_metrics[(tenx_metrics['library_type'] == lt) & (tenx_metrics['metric_name'] == mn)]
    if len(plt_df['sample_id']) == len(plt_df['sample_id'].unique()):
        fig = sns.barplot(plt_df, x='sample_id', y='metric_value', color='grey')
        fig.set_xticklabels(fig.get_xticklabels(), rotation=90)
        fig.set_title(lt + ':' + mn)
        plt.savefig('figures/' + lt + '-' + mn + '.png', bbox_inches='tight')
    else:
        # do a boxplot instead
        fig = sns.boxplot(plt_df, x='sample_id', y='metric_value', color='grey')
        fig.set_xticklabels(fig.get_xticklabels(), rotation=90)
        fig.set_title(lt + ':' + mn)
        plt.savefig('figures/' + lt + '-' + mn + '.png', bbox_inches='tight')
    plt.clf()
