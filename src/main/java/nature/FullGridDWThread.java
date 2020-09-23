/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import edu.uci.ics.jung.graph.util.Pair;
import java.util.ArrayList;
import java.util.concurrent.Callable;

/**
 *
 * @author plantc59cs
 */
public class FullGridDWThread implements Callable<DWUpdateResult> {

    Pair<Integer>[] edges; //the subset of edges this thread updates
    int[] notEdges;
    Parallel p;
    int n;
    int d;
    int counter; //place to write the result
    DWUpdateResult res;
    boolean verbose = false;

    public FullGridDWThread(Pair<Integer>[] edges, int[] notEdges, Parallel p) {
        this.edges = edges;
        this.p = p;
        n = p.n;
        d = p.d;
        this.notEdges = notEdges;
        int internalNotEdges = notEdges.length * (notEdges.length - 1) / 2;
        int maxOutgoingNotEdges = notEdges.length * p.g.getVertexCount() - notEdges.length;
        res = new DWUpdateResult(edges.length + internalNotEdges + maxOutgoingNotEdges);
        counter = 0;
    }

    public DWUpdateResult call() {
        if (verbose) {
            System.out.println("thread  " + this.hashCode() + " started. Processing " + edges.length + " edges.");
        }

        for (int j = 0; j < edges.length; j++) {
//                if (j == 3 && i == 7) {
//                    System.out.println(j);
//                }
            update(edges[j], true, counter);
            counter++;
        }
        includeNotEdges();

        if (verbose) {
            System.out.println(this.hashCode() + " finished.");
        }
        return res;

    }

    private void includeNotEdges() {
        //int counter = edges.length;
        for (int i = 0; i < notEdges.length; i++) {
            includeNotEdges(notEdges[i]);
        }
    }

    private void includeNotEdges(int endpoint) {
////        //DEBUG
//        if(endpoint == 34)
//       System.out.println("ep "+ endpoint);
////        //DEBUG
        int iIndex = Math.min((int) Math.floor((p.coord[endpoint][0] - p.minMax[0][0]) / p.cutOff), p.numBinsX - 1);
        int jIndex = Math.min((int) Math.floor((p.coord[endpoint][1] - p.minMax[1][0]) / p.cutOff), p.numBinsY - 1);
        if (!p.useGrid) {
            iIndex = 0;
            jIndex = 0;
        }

        //points in the same gridcell not connected to endpoint
        for (int i = 0; i < p.pId[iIndex][jIndex].length; i++) {
            // System.out.println(i);
            if ((p.pId[iIndex][jIndex][i] != endpoint) && (!p.g.isNeighbor(p.pId[iIndex][jIndex][i], endpoint))) {
                Pair<Integer> e = new Pair<Integer>(endpoint, p.pId[iIndex][jIndex][i]);
                update(e, false, counter);
                counter++;
            }
        }
        if (p.useGrid == true) {

            //i, j+1
            if (jIndex + 1 <= p.numBinsY - 1) {
                for (int i = 0; i < p.pId[iIndex][jIndex + 1].length; i++) {
                    if (!p.g.isNeighbor(p.pId[iIndex][jIndex + 1][i], endpoint)) {
                        Pair<Integer> e = new Pair<Integer>(endpoint, p.pId[iIndex][jIndex + 1][i]);
                        update(e, false, counter);
                        counter++;
                    }
                }
            }
            //i, j-1
            if (jIndex - 1 >= 0) {
                for (int i = 0; i < p.pId[iIndex][jIndex - 1].length; i++) {
                    if (!p.g.isNeighbor(p.pId[iIndex][jIndex - 1][i], endpoint)) {
                        Pair<Integer> e = new Pair<Integer>(endpoint, p.pId[iIndex][jIndex - 1][i]);
                        update(e, false, counter);
                        counter++;
                    }
                }
            }
            //i-1, j-1
            if (jIndex - 1 >= 0 && iIndex - 1 >= 0) {
                for (int i = 0; i < p.pId[iIndex - 1][jIndex - 1].length; i++) {
                    if (!p.g.isNeighbor(p.pId[iIndex - 1][jIndex - 1][i], endpoint)) {
                        Pair<Integer> e = new Pair<Integer>(endpoint, p.pId[iIndex - 1][jIndex - 1][i]);
                        update(e, false, counter);
                        counter++;
                    }
                }
            }
            //i-1, j
            if (iIndex - 1 >= 0) {
                for (int i = 0; i < p.pId[iIndex - 1][jIndex].length; i++) {
                    if (!p.g.isNeighbor(p.pId[iIndex - 1][jIndex][i], endpoint)) {
                        Pair<Integer> e = new Pair<Integer>(endpoint, p.pId[iIndex - 1][jIndex][i]);
                        update(e, false, counter);
                        counter++;
                    }
                }
            }
            //i-1, j+1
            if (iIndex - 1 >= 0 && jIndex + 1 <= p.numBinsY - 1) {
                for (int i = 0; i < p.pId[iIndex - 1][jIndex + 1].length; i++) {
                    if (!p.g.isNeighbor(p.pId[iIndex - 1][jIndex + 1][i], endpoint)) {
                        Pair<Integer> e = new Pair<Integer>(endpoint, p.pId[iIndex - 1][jIndex + 1][i]);
                        update(e, false, counter);
                        counter++;
                    }
                }
            }
            //i+1, j-1
            if (iIndex + 1 <= p.numBinsX - 1 && jIndex - 1 >= 0) {
                for (int i = 0; i < p.pId[iIndex + 1][jIndex - 1].length; i++) {
                    if (!p.g.isNeighbor(p.pId[iIndex + 1][jIndex - 1][i], endpoint)) {
                        Pair<Integer> e = new Pair<Integer>(endpoint, p.pId[iIndex + 1][jIndex - 1][i]);
                        update(e, false, counter);
                        counter++;
                    }
                }
            }
            //i+1, j
            if (iIndex + 1 <= p.numBinsX - 1) {
                for (int i = 0; i < p.pId[iIndex + 1][jIndex].length; i++) {
                    if (!p.g.isNeighbor(p.pId[iIndex + 1][jIndex][i], endpoint)) {
                        Pair<Integer> e = new Pair<Integer>(endpoint, p.pId[iIndex + 1][jIndex][i]);
                        update(e, false, counter);
                        counter++;
                    }
                }
            }
            //i+1, j+1
            if (iIndex + 1 <= p.numBinsX - 1 && jIndex + 1 <= p.numBinsY - 1) {
                for (int i = 0; i < p.pId[iIndex + 1][jIndex + 1].length; i++) {
                    if (!p.g.isNeighbor(p.pId[iIndex + 1][jIndex + 1][i], endpoint)) {
                        Pair<Integer> e = new Pair<Integer>(endpoint, p.pId[iIndex + 1][jIndex + 1][i]);
                        update(e, false, counter);
                        counter++;
                    }
                }
            }
        }

    }

    public void update(Pair<Integer> edge, boolean connected, int counter) {
        double x = p.dist(p.coord[edge.getFirst()], p.coord[edge.getSecond()]);
        double[] dw = p.s.parabola(x, connected);
//        if(counter == 2488)
//            System.out.println("m");
        res.d[counter] = dw[0];
        res.w[counter] = dw[1];
        res.ij[counter][0] = edge.getFirst();
        res.ij[counter][1] = edge.getSecond();
    }

}
