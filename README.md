# GEMPE
This repository contains Java and C++ implementations of the algorithm GEMPE introduced in the KDD 2020 paper "Data Compression as a Comprehensive Framework for Graph Drawing and Representation Learning" by Claudia Plant, Sonja Biedermann and Christian Böhm. 

The Java Code have been build with MAVEN using IntellJ IDE. Please see the pom.xml file for all necessary libraries. Input files in the form of Matlab matrices and edge lists are supported. See meshes.mat (Matlab matix) and openFlightJ.txt (edge List; first line: number of nodes; subsequent lines of the form "node_i node_j") for example input files. The main classes to run the code are nature.ParallelMain and visu.DisplayEmbedding for embedding without or with visualization. For more details please refer to the appendix of the GEMPE paper.

Please find C++ code in the corresponding folder. This zip folder also includes a readme and an example input data set.

You find the experimental data as Matlab matrices in data.zip.

The folder additionalExperiments contains more experiments and the folder pythonScripts the scripts to reproduce the representation learning experiments in the paper.
