import argparse
import pandas as pd
import re


from panpipes.funcs.processing import  extract_parameter_from_fname

import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)
L.debug("test logging message")

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--input_files_str',
                    default='',
                    help='')
parser.add_argument('--output_file',
                    default='',
                    help='')
parser.add_argument('--sample_prefix',
                    default='',
                    help='')
parser.add_argument('--clusters_or_markers',
                    default='',
                    help='')
parser.set_defaults(verbose=True)
args = parser.parse_args()

L.info(args)

infiles = re.split(',', args.input_files_str)
if args.clusters_or_markers == "clusters":
    combined_csv = pd.concat([pd.read_csv(f, sep='\t', index_col=0) for f in infiles], axis=1)
    # get colnames
    cnames = []
    for f in infiles:
        alg = extract_parameter_from_fname(f, 'alg', prefix=args.sample_prefix)
        res = extract_parameter_from_fname(f, 'res', prefix=args.sample_prefix)
        cnames.append(alg + '_res_' + str(res))
    combined_csv.to_csv(args.output_file, sep='\t', header=cnames, index=True)


if args.clusters_or_markers == "markers":
    li = []
    all_markers_file = re.sub("_top", "_all", args.output_file)
    excel_file = re.sub("_top.txt.gz", "_all.xlsx", args.output_file)
    excel_file_top = re.sub("_top.txt.gz", "_top.xlsx", args.output_file)
    with pd.ExcelWriter(excel_file) as writer:
        with pd.ExcelWriter(excel_file_top) as writer2:
            for ff in infiles:
                # print(os.path.join(in_path, ff))
                df = pd.read_csv(ff, sep='\t')
                # add a cluster column
                clust_val = extract_parameter_from_fname(ff, "cluster", prefix=args.sample_prefix)
                sname = "cluster" + str(clust_val)
                df['cluster'] = clust_val
                li.append(df)
                df.to_excel(writer, sheet_name=sname, index=False)
                df_sub = df[df['p.adj.bonferroni'] < 0.05]
                df_sub = df_sub[df_sub['avg_logFC'] > 0]
                df_sub.to_excel(writer2, sheet_name=sname, index=False)
    frame = pd.concat(li, axis=0, ignore_index=True)
    frame.to_csv(all_markers_file, sep='\t', header=True, index=False)
    frame_sub = frame[frame['p.adj.bonferroni'] < 0.05]
    frame_sub = frame_sub[frame_sub['avg_logFC'] > 0]
    frame_sub.to_csv(args.output_file, sep='\t', header=True, index=False)


L.info("Done")

