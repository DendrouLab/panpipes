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
import imp
import panpipes


def main(argv=None):
    argv = sys.argv
    # paths to look for pipelines:
    path = os.path.abspath(os.path.join(os.path.dirname(panpipes.__file__), "panpipes"))
    relpath = os.path.abspath("./")
    paths = [path, relpath]
    if len(argv) == 1 or argv[1] == "--help" or argv[1] == "-h":
        print("The list of available workflows are:\n")
        pipelines_list = ['1. "ingest" : for the ingestion of data and computation of QC metrics' , 
                          '2. "preprocess" : for filtering and normalising of each modality',
                          '3. "integration" : integrate and batch correction using  single and multimodal methods', 
                          '4. "clustering" : cell clustering on single modalities', 
                          '5. "refmap" : transfer scvi-tools models from published data to your data', 
                          '6. "vis" : visualise metrics from other pipelines in context of experiment metadata']
        print(*pipelines_list, sep="\n")
        return
    command = argv[1]
    command = re.sub("-", "_", command)
    pipeline = "pipeline_{}".format(command)
    del sys.argv[0]
    (file, pathname, description) = imp.find_module(pipeline, paths)
    module = imp.load_module(pipeline, file, pathname, description)
    module.main(sys.argv)


if __name__ == "__main__":
    sys.exit(main())
