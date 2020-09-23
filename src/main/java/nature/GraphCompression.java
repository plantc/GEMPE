/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import edu.uci.ics.jung.graph.Graph;
import java.util.Collection;
import utils.DataUtils;

/**
 *
 * @author claudia assume that vertices in g are labeled corresponding to the
 * order of coordinates.
 */
public class GraphCompression {

    Graph<Integer, Integer> g;
    double[][] coord;
    int numVertices;
    int numEdges;
    int numPossEdges; //possible number of edges
    int dim;
    int numBins;
    double[] edgeCosts;
    double[] nodeWeights;
    double[] edgeMultiplicities;
    SimpleSigmoid s;
    WeightedSigmoid w;
    boolean verbose = true;

    public GraphCompression(Graph g, double[][] coord, int numBins) {
        this.g = g;
        this.coord = coord;
        numVertices = g.getVertexCount();
        numEdges = g.getEdgeCount();
        numPossEdges = (numVertices * (numVertices - 1)) / 2;
        dim = coord[0].length;
        this.numBins = numBins;
    }

    public GraphCompression(Graph g, double[][] coord) {
        this.g = g;
        this.coord = coord;
        numVertices = g.getVertexCount();
        numEdges = g.getEdgeCount();
        numPossEdges = (numVertices * (numVertices - 1)) / 2;
        dim = coord[0].length;
    }

    public GraphCompression(Graph g) {
        this.g = g;
        numVertices = g.getVertexCount();
        numEdges = g.getEdgeCount();
        numPossEdges = (numVertices * (numVertices - 1)) / 2;
    }

    public void setEdgeMultiplicities(double[] edgeMultiplicities) {
        this.edgeMultiplicities = edgeMultiplicities;
    }

    public void setNodeWeights(double[] nodeWeights) {
        this.nodeWeights = nodeWeights;
    }

    public void setCoord(double[][] coord) {
        this.coord = coord;
    }

    public double paramCostBic() {
        return (numVertices / 2 * lg2(numPossEdges)) * dim;
        //return numVertices  * lg2(numPossEdges);
    }

    public double codingCostNoEmbedding() {
        double pEdge = (double) numEdges / (double) numPossEdges;
        double pNoEdge = 1.0 - pEdge;
        double pEdger = Math.round(1000.0 * pEdge) / 1000.0; // Ergebnis: 1.23
        double pNoEdger = Math.round(1000.0 * pNoEdge) / 1000.0; // Ergebnis: 1.23
        System.out.println(pEdger + " " + pNoEdger);
        double entropy = -(pEdge * lg2(pEdge) + pNoEdge * lg2(pNoEdge));
        return entropy * numPossEdges;
    }

    private double probFunction(double dist) {
        return 1.0 / (1 + Math.pow(dist, 2));
    }

    //cost computed with function defined by the paper of Tang et al.
    public double probCost(int[][] relevantNotLinks) {
        double sumLog = 0.0;
        for (int i = 0; i < numVertices; i++) {
            double linkCost = 0.0;
            Collection<Integer> neighbors = g.getNeighbors(i);
            for (int n : neighbors) {
                double dist = euclideanDistance(i, n);
                linkCost += lg2(probFunction(dist));
            }
            sumLog += linkCost;
//            if(Double.isInfinite(linkCost))
//                System.out.println("m");
            double notLinkCost = 0.0;
            for (int j = 0; j < relevantNotLinks[i].length; j++) {
                double dist = euclideanDistance(i, relevantNotLinks[i][j]);
                notLinkCost += lg2(1 - probFunction(dist));
            }
            if(Double.isInfinite(notLinkCost))
                System.out.println("n");
            sumLog += notLinkCost;
        }
        return sumLog;
    }

    public void setEdgeCosts() {
        this.edgeCosts = new double[numPossEdges];
        DataUtils du = new DataUtils();
        int n = g.getVertexCount();
        int m = n * (n - 1) / 2;

        PairOwn[] p = new PairOwn[m];
        double sigmoidCost = 0.0;
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < i; j++) {
                p[i * (i - 1) / 2 + j] = new PairOwn(distance(i, j), g.isNeighbor(i, j));

            }
        }
        SimpleSigmoid s = new SimpleSigmoid(p);
        for (int i = 0; i < g.getVertexCount(); i++) {
            for (int j = 0; j < i; j++) {
                // double aktCost = 0.0;
                if (g.isNeighbor(i, j)) {
                    //   aktCost = s.costEdge(distance(i, j));
                    //sigmoidCost += s.costEdge(distance(i, j));
                    double dist = distance(i, j);
                    if (dist > s.mu) {
                        edgeCosts[du.getIndex(i, j, numVertices)] = 100.0;
                    }

                } else {
                    // aktCost = s.costNoEdge(distance(i, j));
                    //sigmoidCost += s.costNoEdge(distance(i, j));
                }
                //sigmoidCost += aktCost;
//                cost[i][j] = aktCost;
//                cost[j][i] = aktCost;
                // edgeCosts[du.getIndex(i, j, numVertices)] = aktCost;
            }

        }
//        if (verbose) {
//            DataUtils du = new DataUtils();
//            du.saveAsMatlab(cost, "cost", "sigmoidCost.mat");
//        }

    }

    public double mdlFunction() {
        //double[][] cost = new double[numVertices][numVertices];
        this.edgeCosts = new double[numPossEdges];
        DataUtils du = new DataUtils();
        int n = g.getVertexCount();
        int m = n * (n - 1) / 2;

        PairOwn[] p = new PairOwn[m];
        double sigmoidCost = 0.0;
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < i; j++) {
                p[i * (i - 1) / 2 + j] = new PairOwn(distance(i, j), g.isNeighbor(i, j));

            }
        }
        Sigmoid s = new Sigmoid(p);
        for (int i = 0; i < g.getVertexCount(); i++) {
            for (int j = 0; j < i; j++) {
                double aktCost = 0.0;
                if (g.isNeighbor(i, j)) {
                    aktCost = s.costEdge(distance(i, j));
                    //sigmoidCost += s.costEdge(distance(i, j));
                    double dist = distance(i, j);
                    if (dist > s.mu) {
                        edgeCosts[du.getIndex(i, j, numVertices)] = 100.0;
                    }

                } else {
                    aktCost = s.costNoEdge(distance(i, j));
                    //sigmoidCost += s.costNoEdge(distance(i, j));
                }
                sigmoidCost += aktCost;
//                cost[i][j] = aktCost;
//                cost[j][i] = aktCost;
                // edgeCosts[du.getIndex(i, j, numVertices)] = aktCost;
            }

        }
//        if (verbose) {
//            DataUtils du = new DataUtils();
//            du.saveAsMatlab(cost, "cost", "sigmoidCost.mat");
//        }

        return sigmoidCost;

    }

    public void mdlFunctionSimpleSigmoidComparisonMethods() {
        //double[][] cost = new double[numVertices][numVertices];
        this.edgeCosts = new double[numPossEdges];
        DataUtils du = new DataUtils();
        int n = g.getVertexCount();
        int m = n * (n - 1) / 2;

        PairOwn[] p = new PairOwn[m];

        for (int i = 1; i < n; i++) {
            for (int j = 0; j < i; j++) {
                p[i * (i - 1) / 2 + j] = new PairOwn(distance(i, j), g.isNeighbor(i, j));

            }
        }

        SimpleSigmoidComparisonMethods sc = new SimpleSigmoidComparisonMethods(p);
        System.out.println("Sigmoid computed");
        // s.sigma = Math.max(s.sigma, 10.0);
        //DEBUG
        //System.out.println("variance after update " + s.sigma);
        //DEBUG
        double sigmoidCost = sc.costAllPairs(p);
        System.out.println(sc.sigma + " " + sigmoidCost);

        //return sigmoidCost;
    }

    public double mdlFunctionSimpleSigmoid() {
        //double[][] cost = new double[numVertices][numVertices];
        this.edgeCosts = new double[numPossEdges];
        DataUtils du = new DataUtils();
        int n = g.getVertexCount();
        int m = n * (n - 1) / 2;

        PairOwn[] p = new PairOwn[m];

        for (int i = 1; i < n; i++) {
            for (int j = 0; j < i; j++) {
                p[i * (i - 1) / 2 + j] = new PairOwn(distance(i, j), g.isNeighbor(i, j));

            }
        }

        s = new SimpleSigmoid(p);
        // s.sigma = Math.max(s.sigma, 10.0);
        //DEBUG
        //System.out.println("variance after update " + s.sigma);
        //DEBUG
        double sigmoidCost = s.costAllPairs(p);

        return sigmoidCost;

    }

    public double mdlFunctionWeightedSigmoid(double maxSigma) {
        //double[][] cost = new double[numVertices][numVertices];
        this.edgeCosts = new double[numPossEdges];
        DataUtils du = new DataUtils();
        int n = g.getVertexCount();
        int m = n * (n - 1) / 2;

        WeightedPair[] p = new WeightedPair[m];

        for (int i = 1; i < n; i++) {
            for (int j = 0; j < i; j++) {
                double alpha = 0.0;
                double beta = 0.0;
                if (g.isNeighbor(i, j)) {
                    int index = g.findEdge(i, j);
                    alpha = edgeMultiplicities[index];
                    beta = nodeWeights[i] * nodeWeights[j] - alpha;
                } else {
                    beta = nodeWeights[i] * nodeWeights[j];

                }

                p[i * (i - 1) / 2 + j] = new WeightedPair(distance(i, j), alpha, beta);

            }
        }

        w = new WeightedSigmoid(p);
        if (w.sigma > maxSigma) {
            w.sigma = maxSigma;
        }
        // s.sigma = Math.max(s.sigma, 10.0);
        //DEBUG
        //System.out.println("variance after update " + s.sigma);
        //DEBUG
        double sigmoidCost = w.costAllPairs(p);

        return sigmoidCost;

    }

    public double mdlFunctionSimpleSigmoid(boolean cheat) {
        //double[][] cost = new double[numVertices][numVertices];
        this.edgeCosts = new double[numPossEdges];
        DataUtils du = new DataUtils();
        int n = g.getVertexCount();
        int m = n * (n - 1) / 2;

        PairOwn[] p = new PairOwn[m];

        for (int i = 1; i < n; i++) {
            for (int j = 0; j < i; j++) {
                p[i * (i - 1) / 2 + j] = new PairOwn(distance(i, j), g.isNeighbor(i, j));

            }
        }

        SimpleSigmoid s = new SimpleSigmoid(p);
        if (cheat) {
            s.sigma = Math.max(s.sigma, 5.0);
        }
        //DEBUG
        //System.out.println("variance after update " + s.sigma);
        //DEBUG
        double sigmoidCost = s.costAllPairs(p);

        return sigmoidCost;

    }

    //3.7.14: write matrix C for code of Dhillon
    public void writeContraintsMatrix(double[] dist) {
        double[][] c = new double[numPossEdges][4];
        DataUtils du = new DataUtils();
        //determine cutoff
        int n = g.getVertexCount();
        int m = n * (n - 1) / 2;

        PairOwn[] p = new PairOwn[m];
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < i; j++) {
                p[i * (i - 1) / 2 + j] = new PairOwn(dist[du.getIndex(i, j, numVertices)], g.isNeighbor(i, j));

            }
        }
        Sigmoid s = new Sigmoid(p);
        double cutOff = s.mu;
        double stdv = Math.sqrt(s.sigma);
        double smallDist = cutOff - stdv;
        double largeDist = cutOff + stdv;
        int counter = 0;
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < i; j++) {
                if (g.isNeighbor(i, j) && dist[du.getIndex(i, j, numVertices)] > cutOff) {
                    c[counter][0] = i + 1;
                    c[counter][1] = j + 1;
                    c[counter][2] = 1.0;
                    c[counter][3] = smallDist;
                    counter++;
                }
                if (!g.isNeighbor(i, j) && dist[du.getIndex(i, j, numVertices)] < cutOff) {
                    c[counter][0] = i + 1;
                    c[counter][1] = j + 1;
                    c[counter][2] = -1.0;
                    c[counter][3] = largeDist;
                    counter++;
                }

            }
        }
        System.out.println("number of constraints: " + counter);
        du.saveAsMatlab(c, "c", "c.mat");

    }

    //1.7.14: evaluate the cost based on given metric or non-metric distance function
    public double mdlFunction(double[] dist) {
        DataUtils du = new DataUtils();
        double[][] cost = new double[numVertices][numVertices];
        int n = g.getVertexCount();
        int m = n * (n - 1) / 2;
        this.edgeCosts = new double[numPossEdges];

        PairOwn[] p = new PairOwn[m];
        double sigmoidCost = 0.0;
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < i; j++) {
                p[i * (i - 1) / 2 + j] = new PairOwn(dist[du.getIndex(i, j, numVertices)], g.isNeighbor(i, j));

            }
        }
        Sigmoid s = new Sigmoid(p);
        for (int i = 0; i < g.getVertexCount(); i++) {
            for (int j = 0; j < i; j++) {
                double aktCost = 0.0;
                if (g.isNeighbor(i, j)) {
//                    aktCost = s.costEdge(dist[du.getIndex(i, j, numVertices)]);
//                    if(dist[du.getIndex(i, j, numVertices)] > s.mu])
//                        costs[du.getIndex(i, j, numVertices)] = 1.0;
                    //sigmoidCost += s.costEdge(distance(i, j));

                } else {
                    aktCost = s.costNoEdge(dist[du.getIndex(i, j, numVertices)]);
                    //sigmoidCost += s.costNoEdge(distance(i, j));
                }
                sigmoidCost += aktCost;
                cost[i][j] = aktCost;
                cost[j][i] = aktCost;
            }

        }
        if (verbose) {

            du.saveAsMatlab(cost, "cost", "sigmoidCost.mat");
        }

        return sigmoidCost;

    }

    private double distance(int ii, int jj) {
        double dist = 0.0;
        for (int i = 0; i < dim; i++) {
            dist += (coord[ii][i] - coord[jj][i]) * (coord[ii][i] - coord[jj][i]);
        }
        return Math.sqrt(dist);
    }

    public static double sigmoid(double x, double mu, double stddev, double a, double b) {
        return (a - b) * (1 - normp((x - mu) / stddev)) + b;
    }

    public double mdlHisto() {
        double codingCost = codingCostEmbedding();
        double paramCost = paramCostBic();
        double overallCost = codingCost + paramCost;
        double costWithoutEmbedding = codingCostNoEmbedding() * numPossEdges;
        System.out.println("mdl: " + overallCost + "(codingCost: " + codingCost + " paramCost: " + paramCost + "), without embedding: " + costWithoutEmbedding + " savedCost: " + (costWithoutEmbedding - overallCost));
        return codingCost;
    }

    public static double normp(double z) {

        double zabs;
        double p;
        double expntl, pdf;

        final double p0 = 220.2068679123761;
        final double p1 = 221.2135961699311;
        final double p2 = 112.0792914978709;
        final double p3 = 33.91286607838300;
        final double p4 = 6.373962203531650;
        final double p5 = .7003830644436881;
        final double p6 = .3526249659989109E-01;

        final double q0 = 440.4137358247522;
        final double q1 = 793.8265125199484;
        final double q2 = 637.3336333788311;
        final double q3 = 296.5642487796737;
        final double q4 = 86.78073220294608;
        final double q5 = 16.06417757920695;
        final double q6 = 1.755667163182642;
        final double q7 = .8838834764831844E-1;

        final double cutoff = 7.071;
        final double root2pi = 2.506628274631001;

        zabs = Math.abs(z);
        if (z > 37.0) {
            p = 1.0;
            return p;
        }
        if (z < -37.0) {
            p = 0.0;
            return p;
        }
        expntl = Math.exp(-.5 * zabs * zabs);
        pdf = expntl / root2pi;
        if (zabs < cutoff) {
            p = expntl * ((((((p6 * zabs + p5) * zabs + p4) * zabs + p3)
                    * zabs
                    + p2) * zabs + p1) * zabs + p0) / (((((((q7 * zabs
                    + q6) * zabs
                    + q5) * zabs + q4) * zabs + q3) * zabs + q2) * zabs
                    + q1) * zabs
                    + q0);
        } else {
            p = pdf / (zabs + 1.0 / (zabs + 2.0 / (zabs + 3.0 / (zabs + 4.0
                    / (zabs + 0.65)))));
        }
        if (z < 0.0) {
            return p;
        } else {
            p = 1.0 - p;
            return p;
        }
    }

    public double codingCostEmbedding() {
        double maxDist = -Double.MAX_VALUE;
        double minDist = Double.MAX_VALUE;

        double[][] distm = new double[numVertices][numVertices];
        for (int i = 0; i < numVertices; i++) {
            for (int j = 0; j < numVertices; j++) {
                if (i < j) {
                    distm[i][j] = euclideanDistance(i, j);
                    if (distm[i][j] < minDist) {
                        minDist = distm[i][j];
                    }
                    if (distm[i][j] > maxDist) {
                        maxDist = distm[i][j];
                    }
                    distm[j][i] = distm[i][j];
                }

            }
        }

        maxDist *= 1.00000001;

        int[] histEdge = new int[numBins];
        int[] histAll = new int[numBins];
        for (int i = 1; i < numVertices; i++) {
            for (int j = 0; j < i; j++) {
                int bin = (int) ((distm[i][j] - minDist) / (maxDist - minDist) * numBins);
                histAll[bin]++;
                if (g.isNeighbor(i, j)) {
                    histEdge[bin]++;
                }
            }
        }
        double[] pEdge = new double[numBins];
        for (int i = 0; i < numBins; i++) {
            if (histAll[i] > 0) {
                pEdge[i] = (double) histEdge[i] / (double) histAll[i];
            }
        }
        double sumEntropy = 0.0;
        for (int i = 0; i < numBins; i++) {
            if ((pEdge[i] > 0) && (pEdge[i] < 1.0)) {
                double pNotEdge = 1.0 - pEdge[i];
                sumEntropy -= (pEdge[i] * lg2(pEdge[i]) + pNotEdge * lg2(pNotEdge)) * histAll[i];
            }

        }
        return sumEntropy;

    }

    private double lg2(double d) {
        return Math.log(d) / Math.log(2.0);
    }

    private double euclideanDistance(int o, int p) {
        double dist = 0.0;

        for (int i = 0; i
                < dim; i++) {
            dist += (coord[o][i] - coord[p][i]) * (coord[o][i] - coord[p][i]);

        }
        return Math.sqrt(dist);

    }
}
