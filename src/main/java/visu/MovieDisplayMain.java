/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package visu;

import edu.uci.ics.jung.graph.Graph;
import utils.IO;
import utils.DataUtils;
/**
 *
 * @author plantc59cs
 */
public class MovieDisplayMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        IO ea = new IO();
        DataUtils du = new DataUtils();
        Graph g = ea.matlabToGraph("minnesotaL.mat", "graph");
        System.out.println(g.getVertexCount());
        double[][] gtCoord = ea.readMatlabMatrix("groundTruthScaled.mat", "gtScaled");
        Visualization v = new Visualization(g);
        v.displayMovie(gtCoord);

    }

}
