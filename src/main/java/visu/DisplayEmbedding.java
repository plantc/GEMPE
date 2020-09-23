/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package visu;

import Jama.Matrix;
import edu.uci.ics.jung.graph.Graph;
import utils.DataUtils;
import utils.IO;

/**
 *
 * @author claudia.plant
 */
public class DisplayEmbedding {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        IO ea = new IO();
        DataUtils du = new DataUtils();
        Graph g = ea.matlabToGraph("meshes.mat", "smallmesh");
        //  Graph g = ea.matlabToGraph("locationGraph.mat", "graph");
        // g = du.getLargestComponent(g);
        //System.out.println("m");
        // Graph g = ea.matlabToGraph("meshes.mat", "tapir");
        // Graph g = ea.matlabToGraph("meshes.mat", "eppstein");
        //Graph g = ea.matlabToGraph("football.mat", "graph");
        // Graph g = ea.matlabToGraph("reducedGraph.mat", "graph");
        // Graph g = ea.matlabToGraph("netScienceGraph.mat", "graph");


        // Graph g = ea.matlabToGraph("airflights.mat", "graph");
//        Graph g = ea.matlabToGraph("graphLevel23.mat", "graph");
//        double[][] nodeweights = du.readMatlabMatrix("nw23.mat", "nodeweights");
//        double[] nw = new double[nodeweights.length];
//        for(int i = 0; i < nw.length; i++)
//            nw[i] = nodeweights[i][0];
//
//        double[][] edgemultiplicities = du.readMatlabMatrix("nw23.mat", "edgemultiplicities");
//         double[] em = new double[edgemultiplicities.length];
//         for(int i = 0; i < em.length; i++)
//            em[i] = edgemultiplicities[i][0];
//        double[][] initCoord = du.readMatlabMatrix("nw23.mat", "coord");


//Graph g = ea.matlabToGraph("minnesotaL.mat", "graph");
        //Graph g = ea.matlabToGraph("netScienceGraph.mat", "graph");
        //  Graph g = ea.matlabToGraph("dti.mat", "graph");
        //Graph g = ea.matlabToGraph("twoMoons.mat", "g6");

        //Graph g = ea.matlabToGraph("polbooks.mat", "graph");
        //  double[][] coord = ea.readMatlabMatrix("init.mat", "init");
//


        // Graph g = ea.matlabToGraph("smallTest.mat", "graph");
        // Graph g = ea.matlabToGraph("sphere3d.mat", "graph");
        //Graph g = ea.matlabToGraph("meshes.mat", "tapir");
        //Graph g = ea.matlabToGraph("graph.mat", "graph");
//         DataUtils du = new DataUtils();
//         g = du.getLargestComponent(g);
//        double[][] labels = ea.readMatlabMatrix("labels_airflights.mat", "labels");
//        int[] l = new int[labels.length];
//        for(int i  = 0; i < l.length; i++)
//            l[i] = (int) labels[i][0];
//int[] l = new int[g.getVertexCount()];
//for(int i = 0; i < l.length; i++){
//    l[i] = i;
//}
//
//        GraphClustering gc = new GraphClustering(g, 0.9);
//        gc.run();
//        gc.getClusterIds();
//        int[][] ids = gc.clusterIds;
//        int[] numCl = gc.numClusters;
//        int[][] idsd = new int[ids.length+1][ids[0].length];
//        int[] numCld = new int[numCl.length +1];
//        numCld[0] = g.getVertexCount();
//        for(int i = 1; i < numCld.length; i++)
//            numCld[i] = numCl[i-1];
//
//

        //  System.out.println(g.isNeighbor(21, 23));
        Visualization v = new Visualization(g);
//      int[] clId = new int[g.getVertexCount()];
//
//      int numIds = 136;
//      int[] dummyId = new int[g.getVertexCount()];
//      for(int i = 0; i < dummyId.length; i++){
//          dummyId[i] = i;
//      }
//
//      idsd[0] = dummyId;
//      for(int i = 1; i < idsd.length; i++)
//          idsd[i] = ids[i-1];

        // double[][] coord_mj = v.getCoordinatesItMaj(ids[1], numCl[1]);

        //double[][] coord_mj = v.getCoordinatesItMaj(idsd, numCld, l);

        //double[][] coord_mj = v.getCoordinatesItMaj();


        double[][] coord_mj = v.getCoordinatesStochasticGD();
        //du.saveAsMatlab(coord_mj, "coord", "init.mat");


        //     double[][] coord_mj = v.getCoordinatesItMaj(dummyId, numIds);

        //double[][] coord_mj = v.getCoordinatesItMajWithInit(coord);

        //double[][] coord_mj = v.getCoordinatesItMaj(l);
        //double[][] coord_mj = v.getCoordinatesItMaj();

        // double[][] coord_mj = v.getCoordinatesWeighted(initCoord, nw, em);
        //  ea.writeDoubleToMatlab(coord_mj, "finalResult");

//        GraphCompression gk_iso = new GraphCompression(g, coord_mj);
//        double cost_iso = gk_iso.mdlFunction();
//        double[] costs = gk_iso.edgeCosts;
//        System.out.println(cost_iso);

        //v.displayCoordNew(coord_mj, "mj");
    }
}
