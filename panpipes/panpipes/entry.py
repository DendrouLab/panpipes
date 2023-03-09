# entry .py script adapted from COMBAT clusthub repo, originally by Adam Cribbs

'''
Panpipes - Single Cell analysis workflows 
============================================================

Please remember to acknowledge Charlotte Rich-Griffin and Fabiola Curion 
if you use this pipeline in your research

To use a specific workflow, type::
    panpipes <workflow> [workflow options] [workflow arguments]
For this message and a list of available keywords type::
    panpipes --help
To get help for a specify workflow, type::
    panpipes <workflow> --help
'''

import os
import sys
import re
import glob
import imp
import panpipes


def printListInColumns(l, ncolumns):
    '''output list *l* in *ncolumns*.'''
    ll = len(l)

    if ll == 0:
        return

    max_width = max([len(x) for x in l]) + 3
    n = ll // ncolumns
    if ll % ncolumns != 0:
        n += 1

    # build columns
    columns = [l[x * n:x * n + n] for x in range(ncolumns)]

    # add empty fields for missing columns in last row
    for x in range(ncolumns - (len(l) % ncolumns)):
        columns[-(x + 1)].append('')

    # convert to rows
    rows = list(zip(*columns))

    # build pattern for a row
    p = '%-' + str(max_width) + 's'
    pattern = ' '.join([p for x in range(ncolumns)])

    # put it all together
    return '\n'.join([pattern % row for row in rows])

# l = ['a', 'b', 'c', 'd']
# ncolumns = 2
# printListInColumns(l, ncolumns)

def main(argv=None):

    argv = sys.argv

    # paths to look for pipelines:
    #print(pipelines.__file__)
    path = os.path.abspath(os.path.dirname(panpipes.__file__))
    relpath = os.path.abspath("../")

    paths = [path, relpath]

    if len(argv) == 1 or argv[1] == "--help" or argv[1] == "-h":
        pipelines = []
        for path in paths:
            pipelines.extend(glob.glob(os.path.join(path, "pipeline_*.py")))
        print((globals()["__doc__"]))
        print("The list of available workflows are:\n")
        # pipelines_list = sorted([os.path.basename(x)[len("pipeline_"):-len(".py")] for x in pipelines])
        pipelines_list = ['1. "qc_mm" : for the ingestion of data and computation of QC metrics' , 
                          '2. "preprocess" : for filtering and normalising of each modality',
                          '3. "integration" : integrate and batch correction using  single and multimodal methods', 
                          '4. "clustering" : cell clustering on single modalities', 
                          '5. "refmap" : transfer scvi-tools models from published data to your data', 
                          '6. "vis" : visualise metrics from other pipelines in context of experiment metadata']
        print("{}\n".format(
            printListInColumns(pipelines_list,1)))
        return

    command = argv[1]
    command = re.sub("-", "_", command)
    pipeline = "pipeline_{}".format(command)

    # remove 'cellhub' from sys.argv
    del sys.argv[0]

    (file, pathname, description) = imp.find_module(pipeline, paths)

    module = imp.load_module(pipeline, file, pathname, description)

    module.main(sys.argv)


if __name__ == "__main__":
    sys.exit(main())
