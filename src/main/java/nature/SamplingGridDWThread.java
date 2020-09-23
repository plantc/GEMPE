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
public class SamplingGridDWThread implements Callable<DWUpdateResult> {

    Pair<Integer>[] edges; //the subset of edges this thread updates
    Parallel p;
    int n;
    int d;
    Random r;
    DWUpdateResult res;
    boolean verbose = false;

    public SamplingGridDWThread(Pair<Integer>[] edges, Parallel p, Random r) {
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
        int nonEdgeFailed = 0;
        // r = new Random(seed + i);
        for (int j = 0; j < edges.length; j++) {
//                if (j == 3 && i == 7) {
//                    System.out.println(j);
//                }
            update(edges[j], true, counter);
            counter++;
            //sample not edge from the same and adjacent grid cells

            int notEdge1 = sampleNotEdge(edges[j].getFirst());
            if (notEdge1 > 0) {
                Pair<Integer> e = new Pair<Integer>(edges[j].getFirst(), notEdge1);
                update(e, false, counter);
                counter++;
            } else {
                nonEdgeFailed++;
                notEdge1 = sampleRandomNotEdge(edges[j].getFirst());
                Pair<Integer> e = new Pair<Integer>(edges[j].getFirst(), notEdge1);
                update(e, false, counter);

            }
            int notEdge2 = sampleNotEdge(edges[j].getSecond());
            if (notEdge2 > 0) {
                Pair<Integer> e = new Pair<Integer>(edges[j].getSecond(), notEdge2);
                update(e, false, counter);
                counter++;
            } else {
                nonEdgeFailed++;
                notEdge2 = sampleRandomNotEdge(edges[j].getSecond());
                Pair<Integer> e = new Pair<Integer>(edges[j].getSecond(), notEdge2);
                update(e, false, counter);

            }
        }

        if (verbose) {
            System.out.println(this.hashCode() + " finished. For " + nonEdgeFailed + " endpoints no non-edge found.");
        }
//        p.computeSigmoid();
//        Visualization v = new Visualization(p.g);
//        v.displayCoordNew(p.coord, Integer.toString(seed));
        return res;

    }

    private int sampleRandomNotEdge(int endpoint) {
        int notEdge = r.nextInt(n);
        while (p.g.isNeighbor(endpoint, notEdge) || endpoint == notEdge) {
            notEdge = r.nextInt(n);
        }
        return notEdge;
    }

    private int sampleNotEdge(int endpoint) {
        int res = -1;
        int iIndex = Math.min((int) Math.floor((p.coord[endpoint][0] - p.minMax[0][0]) / p.cutOff), p.numBinsX - 1);
        int jIndex = Math.min((int) Math.floor((p.coord[endpoint][1] - p.minMax[1][0]) / p.cutOff), p.numBinsY - 1);
        if (iIndex < 0)
            iIndex = 0;
        if(jIndex < 0)
            jIndex = 0;
        int numTry = 100;
        int numLocal = 5;
        int tt = 0;
        boolean iplusExcluded = false;
        boolean iminusExcluded = false;
        boolean jplusExcluded = false;
        boolean jminusExcluded = false;
        if (iIndex == 0) {
            iminusExcluded = true;
        }
        if (jIndex == 0) {
            jminusExcluded = true;
        }
        if (iIndex == p.numBinsX - 1) {
            iplusExcluded = true;
        }
        if (jIndex == p.numBinsY - 1) {
            jplusExcluded = true;
        }
        boolean success = false;
        while (!success && tt < numTry) {
            int i = r.nextInt(3);
            int j = r.nextInt(3);
            if (!(i == 0 && j == 0)) {
                if (iplusExcluded && !iminusExcluded) {
                    double d = r.nextDouble();
                    if (d < 0.5) {
                        i = 0;
                    } else {
                        i = 2;
                    }
                }
                if (iminusExcluded && !iplusExcluded) {
                    double d = r.nextDouble();
                    if (d < 0.5) {
                        i = 0;
                    } else {
                        i = 1;
                    }
                }
                if (iminusExcluded && iplusExcluded) {
                    i = 0;
                }
                if (jplusExcluded && !jminusExcluded) {
                    double d = r.nextDouble();
                    if (d < 0.5) {
                        j = 0;
                    } else {
                        j = 2;
                    }
                }
                if (jminusExcluded && !jplusExcluded) {
                    double d = r.nextDouble();
                    if (d < 0.5) {
                        j = 0;
                    } else {
                        j = 1;
                    }
                }
                if (jminusExcluded && jplusExcluded) {
                    j = 0;
                }

            }
            if (i == 0 && j == 0) {  //0: same cell, 1: i+1, 2: i-1
                int cc = 0;
                int index = r.nextInt(p.pId[iIndex][jIndex].length);
                while (p.g.isNeighbor(p.pId[iIndex][jIndex][index], endpoint) && cc < numLocal) {
                    index = r.nextInt(p.pId[iIndex][jIndex].length);
                    cc++;
                }
                if (!p.g.isNeighbor(p.pId[iIndex][jIndex][index], endpoint)) {
                    res = p.pId[iIndex][jIndex][index];
                    success = true;
                }
            }
            //i, j+1
            if (i == 0 && j == 1) {
                if (p.pId[iIndex][jIndex + 1].length > 0) {
                    int cc = 0;
                    int index = r.nextInt(p.pId[iIndex][jIndex + 1].length);
                    while (p.g.isNeighbor(p.pId[iIndex][jIndex + 1][index], endpoint) && cc < numLocal) {
                        index = r.nextInt(p.pId[iIndex][jIndex + 1].length);
                        cc++;
                    }
                    if (!p.g.isNeighbor(p.pId[iIndex][jIndex + 1][index], endpoint)) {
                        res = p.pId[iIndex][jIndex + 1][index];
                        success = true;
                    }
                } else {
                    success = false;
                }
            }
            //i, j-1
            if (i == 0 && j == 2) {
                if (p.pId[iIndex][jIndex - 1].length > 0) {
                    int cc = 0;
                    int index = r.nextInt(p.pId[iIndex][jIndex - 1].length);
                    while (p.g.isNeighbor(p.pId[iIndex][jIndex - 1][index], endpoint) && cc < numLocal) {
                        index = r.nextInt(p.pId[iIndex][jIndex - 1].length);
                        cc++;
                    }
                    if (!p.g.isNeighbor(p.pId[iIndex][jIndex - 1][index], endpoint)) {
                        res = p.pId[iIndex][jIndex - 1][index];
                        success = true;
                    }
                } else {
                    success = false;
                }
            }
            //i-1, j-1
            if (i == 2 && j == 2) {
                if (p.pId[iIndex - 1][jIndex - 1].length > 0) {
                    int cc = 0;
                    int index = r.nextInt(p.pId[iIndex - 1][jIndex - 1].length);
//                //Debug
//                System.out.println(index);
                    while (p.g.isNeighbor(p.pId[iIndex - 1][jIndex - 1][index], endpoint) && cc < numLocal) {
                        index = r.nextInt(p.pId[iIndex - 1][jIndex - 1].length);
                        cc++;
                    }
                    if (!p.g.isNeighbor(p.pId[iIndex - 1][jIndex - 1][index], endpoint)) {
                        res = p.pId[iIndex - 1][jIndex - 1][index];
                        success = true;
                    }
                } else {
                    success = false;
                }
            }
            //i-1, j
            if (i == 2 && j == 0) {
                // System.out.println(endpoint);
                if (p.pId[iIndex - 1][jIndex].length > 0) {
                    int cc = 0;
                    int index = r.nextInt(p.pId[iIndex - 1][jIndex].length);
                    while (p.g.isNeighbor(p.pId[iIndex - 1][jIndex][index], endpoint) && cc < numLocal) {
                        index = r.nextInt(p.pId[iIndex - 1][jIndex].length);
                        cc++;
                    }
                    if (!p.g.isNeighbor(p.pId[iIndex - 1][jIndex][index], endpoint)) {
                        res = p.pId[iIndex - 1][jIndex][index];
                        success = true;
                    }
                } else {
                    success = false;
                }
            }
            //i-1, j+1
            if (i == 2 && j == 1) {
                if (p.pId[iIndex - 1][jIndex + 1].length > 0) {
                    int cc = 0;
                    int index = r.nextInt(p.pId[iIndex - 1][jIndex + 1].length);
                    while (p.g.isNeighbor(p.pId[iIndex - 1][jIndex + 1][index], endpoint) && cc < numLocal) {
                        index = r.nextInt(p.pId[iIndex - 1][jIndex + 1].length);
                        cc++;
                    }
                    if (!p.g.isNeighbor(p.pId[iIndex - 1][jIndex + 1][index], endpoint)) {
                        res = p.pId[iIndex - 1][jIndex + 1][index];
                        success = true;
                    }
                } else {
                    success = false;
                }
            }
            //i+1, j-1
            if (i == 1 && j == 2) {
                if (p.pId[iIndex + 1][jIndex - 1].length > 0) {
                    int cc = 0;
                    int index = r.nextInt(p.pId[iIndex + 1][jIndex - 1].length);
                    while (p.g.isNeighbor(p.pId[iIndex + 1][jIndex - 1][index], endpoint) && cc < numLocal) {
                        index = r.nextInt(p.pId[iIndex + 1][jIndex - 1].length);
                        cc++;
                    }
                    if (!p.g.isNeighbor(p.pId[iIndex + 1][jIndex - 1][index], endpoint)) {
                        res = p.pId[iIndex + 1][jIndex - 1][index];
                        success = true;
                    }
                } else {
                    success = false;
                }
            }
            //i+1, j
            if (i == 1 && j == 0) {
                if (p.pId[iIndex + 1][jIndex].length > 0) {
                    int cc = 0;
                    int index = r.nextInt(p.pId[iIndex + 1][jIndex].length);
                    while (p.g.isNeighbor(p.pId[iIndex + 1][jIndex][index], endpoint) && cc < numLocal) {
                        index = r.nextInt(p.pId[iIndex + 1][jIndex].length);
                        cc++;
                    }
                    if (!p.g.isNeighbor(p.pId[iIndex + 1][jIndex][index], endpoint)) {
                        res = p.pId[iIndex + 1][jIndex][index];
                        success = true;
                    }
                } else {
                    success = false;
                }
            }
            //i+1, j+1
            if (i == 1 && j == 1) {
                if (p.pId[iIndex + 1][jIndex + 1].length > 0) {
                    int cc = 0;
                    int index = r.nextInt(p.pId[iIndex + 1][jIndex + 1].length);
//                //Debug
//                System.out.println(index);
                    while (p.g.isNeighbor(p.pId[iIndex + 1][jIndex + 1][index], endpoint) && cc < numLocal) {
                        index = r.nextInt(p.pId[iIndex + 1][jIndex + 1].length);
                        cc++;
                    }
                    if (!p.g.isNeighbor(p.pId[iIndex + 1][jIndex + 1][index], endpoint)) {
                        res = p.pId[iIndex + 1][jIndex + 1][index];
                        success = true;
                    }
                } else {
                    success = false;
                }
            }

        }
        return res;

    }
    //        //same cell
    //        int iIndex = Math.min((int) Math.floor((p.coord[endpoint][0] - p.minMax[0][0]) / p.cutOff), p.numBinsX - 1);
    //        int jIndex = Math.min((int) Math.floor((p.coord[endpoint][1] - p.minMax[1][0]) / p.cutOff), p.numBinsY - 1);
    //        int ij = p.pId[iIndex][jIndex].length;
    //
    //        // i+1, j
    //        if (iIndex + 1 < p.numBinsX) {
    //            int iplusj = p.pId[iIndex + 1][jIndex].length;
    //            pRelevant += iplusj;
    //            if (jIndex - 1 >= 0) {
    //                // i+1, j-1
    //                int iplusjminus = p.pId[iIndex + 1][jIndex - 1].length;
    //                pRelevant += iplusjminus;
    //            }
    //        }
    //
    //        // i-1, j
    //        int iminusj = p.pId[iIndex - 1][jIndex].length;
    //        pRelevant += iminusj;
    //        // i-1, j-1
    //        int iminusjminus = p.pId[iIndex - 1][jIndex - 1].length;
    //        pRelevant += iminusjminus;
    //
    //        // i, j+1
    //        int ijplus = p.pId[iIndex][jIndex + 1].length;
    //        pRelevant += ijplus;
    //        // i, j-1
    //        int ijminus = p.pId[iIndex][jIndex - 1].length;
    //        pRelevant += ijminus;
    //
    //        //i+1, j+1
    //        int iplusjplus = p.pId[iIndex + 1][jIndex + 1].length;
    //        pRelevant += iplusjplus;

    public void update(Pair<Integer> edge, boolean connected, int counter) {
        double x = p.dist(p.coord[edge.getFirst()], p.coord[edge.getSecond()]);
        double[] dw = p.s.parabola(x, connected);
        res.d[counter] = dw[0];
        res.w[counter] = dw[1];
        res.ij[counter][0] = edge.getFirst();
        res.ij[counter][1] = edge.getSecond();
    }

}
