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
public class SamplingDWThread implements Callable<DWUpdateResult> {

    Pair<Integer>[] edges; //the subset of edges this thread updates
    Parallel p;
    int n;
    int d;
    Random r;
    DWUpdateResult res;
    boolean verbose = false;

    public SamplingDWThread(Pair<Integer>[] edges, Parallel p, Random r) {
        this.edges = edges;
        this.p = p;
        this.r = r;
        n = p.n;
        d = p.d;
        res = new DWUpdateResult(3 * edges.length);

    }

    public DWUpdateResult call() {
        if (verbose) {
            System.out.println("thread  " + this.hashCode() + " started. Processing " + edges.length + " edges.");
        }
        int counter = 0;
        // r = new Random(seed + i);
        for (int j = 0; j < edges.length; j++) {
//                if (j == 3 && i == 7) {
//                    System.out.println(j);
//                }
            update(edges[j], true, counter);
            counter++;
            int notEdge1 = r.nextInt(n);

            while (p.g.isNeighbor(edges[j].getFirst(), notEdge1) || edges[j].getFirst() == notEdge1) {
                notEdge1 = r.nextInt(n);
            }
            Pair<Integer> e = new Pair<Integer>(edges[j].getFirst(), notEdge1);
            update(e, false, counter);
            counter++;
            int notEdge2 = r.nextInt(n);

            while (p.g.isNeighbor(edges[j].getSecond(), notEdge2) || edges[j].getSecond() == notEdge2) {
                notEdge2 = r.nextInt(n);
            }
            e = new Pair<Integer>(edges[j].getSecond(), notEdge2);
            update(e, false, counter);
            counter++;
        }

        if (verbose) {
            System.out.println(this.hashCode() + " finished.");
        }
//        p.computeSigmoid();
//        Visualization v = new Visualization(p.g);
//        v.displayCoordNew(p.coord, Integer.toString(seed));
        return res;

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
