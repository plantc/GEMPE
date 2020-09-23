/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import edu.uci.ics.jung.graph.util.Pair;
import java.util.Random;
import java.util.concurrent.Callable;

/**
 *
 * @author plantc59cs
 */
public class CoordThread implements Callable<CoordUpdateResult> {

    int[] nodes; //the subset of nodes this thread updates, nodenames
    double[][] positions;
    Parallel p;
    CoordUpdateResult res;
    int n;
    int d;
    boolean verbose = false;

    public CoordThread(int[] nodes, Parallel p) {
        n = p.n;
        d = p.d;
        this.nodes = nodes;
        positions = new double[nodes.length][d];
        this.p = p;

    }

    public CoordUpdateResult call() {
        if (verbose) {
            System.out.println("thread  " + this.hashCode() + " started. Processing " + nodes.length + " nodes.");
        }
        for (int i = 0; i < nodes.length; i++) {
            iterate(nodes[i], i);
        }

        if (verbose) {
            System.out.println(this.hashCode() + " finished.");
        }
//        p.computeSigmoid();
//        Visualization v = new Visualization(p.g);
//        v.displayCoordNew(p.coord, Integer.toString(seed));
        CoordUpdateResult cu = new CoordUpdateResult(nodes, positions);
        return cu;
    }

    //get new coordinates for a single point; implement this in the embedding display class
    public void iterate(int i, int index) {
        double sumWeights = 0.0;
        for (int j = 0; j < n; j++) {
            if (i != j) {
//                    if(j == 23)
//                        System.out.println("m");
                sumWeights += p.weights[i][j];
                double sij = 0;
                for (int k = 0; k < d; k++) {
                    sij += Math.pow((p.coord[i][k] - p.coord[j][k]), 2);
                }
                //   sij = Math.max(sij, 1e-200);
                if (sij != 0) {
                    sij = p.distances[i][j] / Math.sqrt(sij);
                }

                for (int k = 0; k < d; k++) {
                    positions[index][k] += p.weights[i][j] * (p.coord[j][k] + sij * (p.coord[i][k] - p.coord[j][k]));
                    if (Double.isNaN(positions[index][0])) {
                        System.out.println("Coord Thread: index " + index + " k: " + k + " j: " + j);
                    }
                }
            }
        }//j
        for (int k = 0; k < d; k++) {
            positions[index][k] /= sumWeights;
        }

    }

}
