/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package utils;
import utils.IO;
import utils.DataUtils;

import edu.uci.ics.jung.graph.Graph;

/**
 *
 * @author plantc59cs
 */
public class IOMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws InterruptedException {

        // TODO code application logic here
        IO ea = new IO();
        DataUtils du = new DataUtils();
        Graph g = ea.readEdgeList("openFlightJ.txt");

        //Graph g = ea.matlabToGraph("meshes.mat", "smallmesh"); //random 70 unverfaltet; ist es durch Grid schlechter?


//Graph g = ea.matlabToGraph("meshes.mat", "tapir");
        //Graph g = ea.matlabToGraph("actFull.mat", "graphFull"); //random 70 unverfaltet; ist es durch Grid schlechter?
        //Graph g = ea.matlabToGraph("dti.mat", "graph"); //random 70 unverfaltet; ist es durch Grid schlechter?
        //Graph g = ea.matlabToGraph("g10frisiert.mat", "graph"); //random 70 unverfaltet; ist es durch Grid schlechter?
        //g = du.getLargestComponent(g);
        //Graph g = ea.matlabToGraph("football.mat", "graph"); //random 70 unverfaltet; ist es durch Grid schlechter?
        // Graph g = ea.matlabToGraph("airflights.mat", "graph"); //random 70 unverfaltet; ist es durch Grid schlechter?
        //Graph g = ea.matlabToGraph("netScienceGraph.mat", "graph");
        //Graph g = ea.matlabToGraph("meshes.mat", "eppstein");
        //Graph g = ea.matlabToGraph("dblp.mat", "graph");
        //double[][] init = ea.readMatlabMatrix("init.mat", "init");

        // double[][] init = ea.readMatlabMatrix("coordi_4490.mat", "coord");
        // double[][] init = ea.readMatlabMatrix("coordi_100.mat", "coord");


        // Graph g = ea.matlabToGraph("reducedGraph.mat", "graph");


        System.out.println(" finished");
    }
}
