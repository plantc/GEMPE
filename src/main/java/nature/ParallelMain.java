/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;
import utils.IO;
import utils.DataUtils;

import edu.uci.ics.jung.graph.Graph;

/**
 *
 * @author plantc59cs
 */
public class ParallelMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws InterruptedException {

        // TODO code application logic here
        IO ea = new IO();
        //DataUtils du = new DataUtils();
        //Graph g = ea.matlabToGraph("minnesotaL.mat", "graph");
        //Graph g = ea.matlabToGraph("meshes.mat", "smallmesh"); //random 70 unverfaltet; ist es durch Grid schlechter?
        //Graph g = ea.readEdgeList("californiaJ.txt");
        Graph g = ea.readEdgeList("oldenburgJ.txt");




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

        //   double[][] labels = ea.readMatlabMatrix("labels_airflights.mat", "labels");
        System.out.println("numNodes " + g.getVertexCount() + " numEdges: " + g.getEdgeCount());


        int d = 2;
        int numT = 2;
        Parallel gp = new Parallel(g, numT, d);
        // gp.twoDGridStart();
        gp.twoD(); //--this is to run Gempe
       // gp.threeD();
        //  gp.initForMovie();
        //  gp.twoDFinishOnly(init);

//        Visualization v = new Visualization(g);
//        int numT = 2;
//        int iterSigmoid = 100;
//        int localMinEscape = 3;
//        int numInits = 10;
//        boolean visualize = true;
//        boolean save = false;
//        int d = 2;
//        Parallel gp = new Parallel(g, numT, d);
//        double[][] init = gp.getBestInitialization(numInits, iterSigmoid, localMinEscape);
//        if (save) {
//            du.saveAsMatlab(init, "init", "init.mat");
//        }
//        if (visualize) {
//            v.displayCoordNew(init, " ");
//        }
//        gp.setCoord(init);
//        if(d == 2)
//        gp.startWithGrid();
//        double[][] stab = gp.stabilizationPhase(iterSigmoid, localMinEscape);
//
//         if(save)
//         du.saveAsMatlab(stab, "stab", "stab.mat");
//        if(visualize){
//             v.displayCoordNew(stab, " ");
//        }
//
//        gp.setCoord(init);
//        if (d == 2) {
//            gp.startWithGrid();
//        }
//        gp.useGrid = false;
//        double[][] res = gp.finishPhase(iterSigmoid, localMinEscape, false);
//        if (save) {
//            du.saveAsMatlab(res, "res", "res.mat");
//        }
//        if (visualize) {
//            v.displayCoordNew(res, " ");
//        }
//    }
        System.out.println("finished");
    }
}
