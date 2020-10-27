# GEMPE
The Java Code have been build with MAVEN using IntellJ IDE. Please see the pom.xml file for all necessary libraries. Input files in the form of Matlab matrices and edge lists are supported. See meshes.mat (Matlab matix) and openFlightJ.txt (edge List; first line: number of nodes; subsequent lines of the form "node_i node_j") for example input files. The main classes to run the code are nature.ParallelMain and visu.DisplayEmbedding for embedding without or with visualization. For more details please refer to the appendix of the GEMPE paper.

Please find C++ code in the corresponding folder. This zip folder also includes a readme and an example input data set.

You find the experimental data as Matlab matrices in data.zip.

The folder additionalExperiments contains more experiments and the folder pythonScripts the scripts to reproduce the representation learning experiments in the paper (https://dl.acm.org/doi/10.1145/3394486.3403174, open access).

## Citation 

    @inproceedings{DBLP:conf/kdd/PlantBB20,
        author    = {Claudia Plant and
                     Sonja Biedermann and
                     Christian B{\"{o}}hm},
        title     = {Data Compression as a Comprehensive Framework for Graph Drawing and
                     Representation Learning},
        booktitle = {{KDD} '20: The 26th {ACM} {SIGKDD} Conference on Knowledge Discovery
                     and Data Mining, Virtual Event, CA, USA, August 23-27, 2020},
        pages     = {1212--1222},
        year      = {2020},
        crossref  = {DBLP:conf/kdd/2020},
        url       = {https://dl.acm.org/doi/10.1145/3394486.3403174},
        timestamp = {Mon, 24 Aug 2020 14:13:33 +0200},
        biburl    = {https://dblp.org/rec/conf/kdd/PlantBB20.bib},
        bibsource = {dblp computer science bibliography, https://dblp.org}
    }
