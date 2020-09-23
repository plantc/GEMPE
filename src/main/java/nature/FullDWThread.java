/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import edu.uci.ics.jung.graph.util.Pair;
import java.util.ArrayList;
import java.util.Collection;
import java.util.concurrent.Callable;

/**
 *
 * @author plantc59cs
 */
public class FullDWThread implements Callable<DWUpdateResult> {

    Pair<Integer>[] edges; //the subset of edges this thread updates
    int[] notEdges;
    Parallel p;
    int n;
    int d;
    DWUpdateResult res;
    boolean verbose = false;

    public FullDWThread(Pair<Integer>[] edges, int[] notEdges, Parallel p) {
        this.edges = edges;
        this.notEdges = notEdges;
        this.p = p;
        n = p.n;
        d = p.d;
        int internalNotEdges = notEdges.length * (notEdges.length - 1) / 2;
        int maxOutgoingNotEdges = notEdges.length * p.g.getVertexCount() - notEdges.length;
        res = new DWUpdateResult(edges.length + internalNotEdges + maxOutgoingNotEdges);
    }

    public DWUpdateResult call() {
        if (verbose) {
            System.out.println("thread  " + this.hashCode() + " started. Processing " + edges.length + " edges.");
        }
        int counter = 0;

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
        int counter = edges.length;
        //connections within pointset of the thread
        for (int i = 0; i < notEdges.length; i++) {
            for (int j = 0; j < i; j++) {
                if (!p.g.isNeighbor(notEdges[i], notEdges[j])) {
                    Pair e = new Pair<Integer>(notEdges[i], notEdges[j]);
                    update(e, false, counter);
                    counter++;
                }
            }
            //connections of i to other points processed by different threads, consider only points with smaller id than the current ones
            for (int j = 0; j < notEdges[i]; j++) {
                if (!p.g.isNeighbor(notEdges[i], j)) {
                    Pair e = new Pair<Integer>(notEdges[i], j);
                    update(e, false, counter);
                    counter++;
                }
            }

        }

    }

    public void update(Pair<Integer> edge, boolean connected, int counter) {
        double x = p.dist(p.coord[edge.getFirst()], p.coord[edge.getSecond()]);
        double[] dw = p.s.parabola(x, connected);
        res.d[counter] = dw[0];
        res.w[counter] = dw[1];
        res.ij[counter][0] = edge.getFirst();
        res.ij[counter][1] = edge.getSecond();
    }

}
