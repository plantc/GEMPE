package visu;

//package edu.uci.ics.jung.algorithms.layout;
import Jama.Matrix;
import edu.uci.ics.jung.algorithms.layout.AbstractLayout;
import edu.uci.ics.jung.algorithms.layout.GraphElementAccessor;
import edu.uci.ics.jung.algorithms.layout.util.RandomLocationTransformer;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraShortestPath;
import edu.uci.ics.jung.algorithms.shortestpath.Distance;
import edu.uci.ics.jung.algorithms.util.IterativeContext;
import org.apache.commons.math3.special.Erf;

import java.awt.Dimension;
import java.awt.geom.Point2D;
import java.util.ConcurrentModificationException;

//import edu.uci.ics.jung.algorithms.*; //import edu.uci.ics.jung.algorithms.GraphStatistics;
//import edu.uci.ics.jung.algorithms.IterativeContext;
//import edu.uci.ics.jung.algorithms.shortestpath.Distance;
//import edu.uci.ics.jung.algorithms.shortestpath.UnweightedShortestPath;
//import edu.uci.ics.jung.algorithms.util.RandomLocationTransformer;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.util.Pair;
import edu.uci.ics.jung.visualization.VisualizationViewer;
import edu.uci.ics.jung.visualization.control.AbstractPopupGraphMousePlugin;
import edu.uci.ics.jung.visualization.control.DefaultModalGraphMouse;
import java.awt.Color;
import java.awt.Paint;
import java.awt.event.ActionEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.IOException;
import java.util.List;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.AbstractAction;
import javax.swing.JPopupMenu;
import nature.SimpleSigmoid;
import nature.GraphCompression;
import nature.PairOwn;
import nature.Sigmoid;
import utils.DataUtils;
import utils.IO;

import org.apache.commons.collections15.Transformer;

/**
 * Implements the Kamada-Kawai algorithm for node layout. Does not respect
 * filter calls, and sometimes crashes when the view changes to it.
 *
 * @see "Tomihisa Kamada and Satoru Kawai: An algorithm for drawing general
 * indirect graphs. Information Processing Letters 31(1):7-15, 1989"
 * @see "Tomihisa Kamada: On visualization of abstract objects and relations.
 * Ph.D. dissertation, Dept. of Information Science, Univ. of Tokyo, Dec. 1988."
 *
 * @author Masanori Harada
 */
public class MovieDisplay extends AbstractLayout<Integer, Integer> implements IterativeContext {

    private int currentIteration;
    private int maxTry = 100;
    private double maxWeightChange = 5.0;
    private int loop;
    private int bestloop;
    private double bestCost; //overall best cost
    private double lastCost = Double.MAX_VALUE; //cost in previous interation
    private double lastSigma; // sigma in previous iteration.
    private boolean cheatVar; // enlarge variance in next iteration
    private int cheatCounter; //how often cheating can happen until converged
    private int d = 2;
    private String status = "OwnLayout";
    private boolean adjustForGravity = true;
    private int[] vertices;
    private Point2D[] xydata;
    double costWithoutEmbedding;
    double paramCost;
    double savedBits;
    double[][] bestDb;
    double[][] distances;
    double[][] coord;
    double[][] weights;
    double[][] groundTruth; //ground truth coordinates
    boolean isoInit = false;
    boolean mdsInit = false;
    boolean groundTruthInit = true;
    boolean randomInit = false;
    int convergeCounter = 0;
    boolean varianceFree = false;
    Graph g;
    GraphCompression gk;
    SimpleSigmoid s;
    Random r;
    // int d;
    int n;
    int m;
    //double factor; //how much to keep of old coordinate
    Dimension dd; //size of display
    private static double PRECISION = 1E-12;
    static int maxIteration = 6990;

    /**
     * Retrieves graph distances between vertices of the visible graph
     */
    public void setXydata(Point2D[] xydata) {
        this.xydata = xydata;
        scaleToDisplaySize();
    }

    public MovieDisplay(Graph<Integer, Integer> g) {
        super(g);
        this.g = g;
        r = new Random(20);
        //this.factor = factor;
        n = g.getVertexCount();
        m = n * (n - 1) / 2;

    }

    public double getBestllh() {
        return bestCost;
    }

    public double[][] getCoordinates() {

        return bestDb;
    }

    public void setGroundTruth(double[][] groundTruth) {
        this.groundTruth = groundTruth;
    }

    public MovieDisplay(Graph<Integer, Integer> g, double factor) {
        super(g);
        r = new Random(20);
        n = g.getVertexCount();
        m = n * (n - 1) / 2;
    }

    /**
     *
     *
     *
     *
     *
     * public String getStatus() { return status + this.getSize(); }
     *
     * public void setMaxIterations(int maxIterations) { this.maxIterations =
     * maxIterations; }
     *
     * /**
     * This one is an incremental visualization.
     */
    public boolean isIncremental() {
        return true;
    }

    /**
     * Returns true once the current iteration has passed the maximum count.
     */
//    public boolean done() {
//        if ((loop - bestloop) < maxTry) {
//            return false;
//        } else {
//            return true;
//        }
//    }
//
    public boolean done() {
        if (currentIteration > maxIteration) {
            return true;
        }
        return false;
    }

    public static double dist(double[] x, double[] y) {
        int d = x.length;
        double result = 0;
        for (int i = 0; i < d; i++) {
            result += (x[i] - y[i]) * (x[i] - y[i]);
        }
        return Math.sqrt(result);
    }

    public void initialize() {
        // r = new Random(1);
        currentIteration = 0;
        bestloop = -1;
        loop = 0;
        dd = getSize();
        currentIteration = 0;
        cheatCounter = 0;
        cheatVar = false;
        vertices = new int[n];
        for (int i = 0; i < n; i++) {
            vertices[i] = i;
        }

        coord = new double[n][d];
        //init coord with isomap
        if (isoInit) {

            distances = pathdist();
            //coord = new Matrix(MDSJ.classicalScaling(distances)).transpose().getArrayCopy();
        }
        if (groundTruthInit) {
            coord = groundTruth;

//             try {
//
//                    System.in.read();
//                      wait(300);
//                        scaleToDisplaySize();
//                        repaint();
//                } catch (Exception ex) {
//                    Logger.getLogger(MovieDisplay.class.getName()).log(Level.SEVERE, null, ex);
//                }

            distances = new double[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    distances[i][j] = distances[j][i] = dist(coord[i], coord[j]);
                }
            }

        }

        if (randomInit) {
            Random r = new Random(1);
            for (int i = 0; i < coord.length; i++) {
                for (int j = 0; j < d; j++) {
                    coord[i][j] = r.nextDouble();
                }
            }
            DataUtils du = new DataUtils();
            //du.saveAsMatlab(coord, "randomInit", "randomInit.mat");
            du.saveResult(coord, "randomInit.txt");
            distances = new double[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    distances[i][j] = distances[j][i] = dist(coord[i], coord[j]);
                }
            }
        }

        if (mdsInit) {
            distances = new double[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    if (g.isNeighbor(i, j)) {
                        distances[i][j] = distances[j][i] = 1.0;
                    } else {
                        distances[i][j] = distances[j][i] = 2.0;
                    }
                }
              //  coord = new Matrix(MDSJ.stressMinimization(distances)).transpose().getArrayCopy();
            }

        }

        weights = new double[n][n];
//        for (int i = 0; i < n; i++) {
//            for (int j = 0; j < i; j++) {
//                weights[i][j] = weights[j][i] = 1.0;
//            }
//        }

//
//        gk = new GraphCompression(g, coord);
//        costWithoutEmbedding = gk.codingCostNoEmbedding();
//        System.out.println("costWithoutEmbedding: " + costWithoutEmbedding);

        xydata = new Point2D[n];
        bestDb = new double[n][d];
        for (int i = 0; i < n; i++) {
            xydata[i] = transform(vertices[i]);
            xydata[i].setLocation(coord[i][0], coord[i][1]);
            bestDb[i][0] = coord[i][0];
            bestDb[i][1] = coord[i][1];
        }
//        computeSigmoid();
//        gk.setCoord(coord);
//        loop = 0;
//        bestCost = gk.mdlFunctionSimpleSigmoid();
//        paramCost = gk.paramCostBic();
//
//        System.out.println("pCostsBic: " + paramCost);
//        savedBits = costWithoutEmbedding - (bestCost + paramCost);
//        System.out.println("init: coding Costs " + bestCost + " saved bits: " + savedBits);
        scaleToDisplaySize();
    }



    private double[][] pathdist() {
        DijkstraShortestPath<Integer, Integer> alg = new DijkstraShortestPath(g);
        double[][] dist = new double[n][n];
        double maxDist_mds = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i > j) {
                    List<Integer> l = alg.getPath(i, j);
                    if (l.size() > 0) {
                        dist[i][j] = l.size();
                        dist[j][i] = dist[i][j];
                        if (dist[i][j] > maxDist_mds) {
                            maxDist_mds = dist[i][j];
                        }
//                    } else {
//                        dist[i][j] = maxDist;
//                        dist[j][i] = maxDist;
//                    }
                    }
                }
            }
        }
        for (int i = 0; i < dist.length; i++) {
            for (int j = 0; j < dist.length; j++) {
                if (dist[i][j] == 0.0 && i != j) {
                    dist[i][j] = maxDist_mds + 1;
                }
            }
        }
        return dist;
    }

    //scales symmetric distance matrix and weight matrix between 0...1. Assumes positive values only
    private double[][] scale(double[][] d) {
        double[][] ds = new double[n][n];
        double max = -Double.MAX_VALUE;
        double min = Double.MAX_VALUE;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                if (d[i][j] > max) {
                    max = d[i][j];
                }
                if (d[i][j] < min) {
                    min = d[i][j];
                }
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                double s = (d[i][j] - min) / (max - min);
                s *= 100;
                ds[i][j] = ds[j][i] = Math.max(s, 0.000000001);
            }
        }

        return ds;
    }

    private void computeSigmoid() {
        //compute sigmoid
        PairOwn[] p = new PairOwn[m];
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < i; j++) {
                p[i * (i - 1) / 2 + j] = new PairOwn(dist(coord[i], coord[j]), graph.isNeighbor(i, j));
            }
        }
        s = new SimpleSigmoid(p);
    }

    private double getAvgDist() {
        double avgDist = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i < j) {
                    avgDist += dist(coord[i], coord[j]) / n;
                }
            }
        }
        return avgDist;
    }


    private synchronized void update() {
        try {

            DataUtils du = new DataUtils();
            int wait = 0;
            if (currentIteration == 2) {

                try {

                    System.in.read();
                    wait(300);
                } catch (Exception ex) {
                    Logger.getLogger(MovieDisplay.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
            if (currentIteration < 5) {
                wait = 300;
            }
            if (currentIteration  >= 5 && currentIteration < 10000) {
                // wait = 100;
                wait=50;
            }


            if (currentIteration < 10) {

                currentIteration++;
            } else {

                currentIteration = currentIteration + 20;
            }
            if(currentIteration > 100)
                currentIteration = currentIteration + 80;

            if(wait > 0){
                try {
                    wait(wait);
                } catch (InterruptedException ex) {
                    Logger.getLogger(Visualization.class.getName()).log(Level.SEVERE, null, ex);
                }
            }


            if (currentIteration <= maxIteration) {
                String s = "coordi_" + Integer.toString(currentIteration) + ".mat";
                if(currentIteration == 0 || currentIteration == 1 || currentIteration == 2)
                    s = "coordi_0.mat";
                System.out.println(s);
                double[][] cc = du.readMatlabMatrix(s, "coord");
                for (int i = 0; i < coord.length; i++) {
                    coord[i][0] = cc[i][0];
                    coord[i][1] = cc[i][1];
                }
                //coord = cc;

                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < d; j++) {
                        // bestDb[i][j] = coord[i][j];
                    }
                    xydata[i].setLocation(coord[i][0], coord[i][1]);
                }

                scaleToDisplaySize();
                //Test

                //String name22 = "coord_" + loop;
                //ea.writeDoubleToMatlab(coord, name22);
                //Test
//            if (loop == 1) {
//                String fName = name + ".mat";
//
//                du.saveAsMatlab4(distances, weights, du.colScaleData(coord), du.colScaleData(groundTruth), "dist", "weights", "coord", "gt", fName);
//            }
                // currentIteration++;
            } else {

            }

            loop++;
        } catch (ConcurrentModificationException cme) {
        }

    }

    public double getDistThreshold(double sigma) {
        double w = 1e-6;
        double arg = -Math.PI * Math.pow(Math.log(2), 2) * Math.pow(sigma, 4) * Math.pow(w, 2);
        double lower = Math.exp(0.5 * LambertW(arg));
        double upper = Math.sqrt(2) * w * Math.pow(sigma, 3) * Math.sqrt(Math.PI) * Math.log(2) - 1.5 * lower;
        return -upper / lower;

    }

    public static double LambertW(double z) {
        double S = 0.0;
        for (int n = 1; n <= 100; n++) {
            double Se = S * Math.pow(Math.E, S);
            double S1e = (S + 1)
                    * Math.pow(Math.E, S);
            if (PRECISION > Math.abs((z - Se) / S1e)) {
                return S;
            }
            S
                    -= (Se - z) / (S1e - (S + 2) * (Se - z) / (2 * S + 2));
        }
        return S;
    }

    private synchronized void updateOld() {
        try {
            double[][] weights_new = new double[n][n];
            double sumWeights_new = 0.0;
            double sumWeights_old = 0.0;
            Sigmoid s = new Sigmoid(0.5);

            PairOwn[] p = new PairOwn[m];
            for (int i = 1; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    p[i * (i - 1) / 2 + j] = new PairOwn(dist(coord[i], coord[j]), graph.isNeighbor(i, j));
                }
            }
            //if (loop%10==0)
            s = new Sigmoid(p);
            double sumllh = 0;
            for (int i = 0; i < g.getVertexCount(); i++) {
                for (int j = 0; j < i; j++) {
                    double h = 0;
                    if (graph.isNeighbor(i, j)) {
                        h = -log2(s.f(dist(coord[i], coord[j])));
                    } else {
                        h = -log2(1.0 - s.f(dist(coord[i], coord[j])));
                    }
                    //weights[i][j] = weights[j][i] = weights[j][i] + h;
                    weights_new[i][j] = weights_new[j][i] = h;
                    sumWeights_new += h;
                    sumWeights_old += weights[i][j];
                    sumllh += h;
                }
            }
            double weightChange = sumWeights_new / sumWeights_old;
            double factor = Math.max(weightChange, maxWeightChange);
            double meanWeightsNew = 0.0;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    weights[i][j] = weights[j][i] = weights[j][i] + factor * weights_new[i][j];
                    meanWeightsNew += (weights[i][j] / (double) m);

                }
            }
            double stdWeightsNew = 0.0;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    stdWeightsNew += ((weights[i][j] - meanWeightsNew) * (weights[i][j] - meanWeightsNew)) / (double) m;
                }
            }
            stdWeightsNew = Math.sqrt(stdWeightsNew);

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    double threshold = meanWeightsNew + stdWeightsNew;
                    if (weights[i][j] > threshold) {
                        weights[i][j] = weights[j][i] = 1.0;
                    }
                }
            }

            if (sumllh < bestCost) {
                bestCost = sumllh;
                savedBits = costWithoutEmbedding - (sumllh + paramCost);
                bestloop = loop;
                System.out.println(loop + ", " + sumllh + " saved bits: " + savedBits);
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < d; j++) {
                        bestDb[i][j] = coord[i][j];
                    }
                    xydata[i].setLocation(coord[i][0], coord[i][1]);
                }
                scaleToDisplaySize();
                //Test
                IO ea = new IO();
                String name = "coord_" + loop;
                ea.writeDoubleToMatlab(coord, name);
                //Test

            }

            for (int innerloop = 0; innerloop < 10; innerloop++) {
                double[][] dbRestoredNew = new double[n][d];
                double[][] sij = new double[n][n];
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        sij[i][j] = 0;
                        double h = dist(coord[i], coord[j]);
                        if (h != 0) {
                            sij[i][j] = distances[i][j] / h;
                        }
                    }
                }
                for (int i = 0; i < n; i++) {
                    double sumWeight = 0;
                    for (int j = 0; j < n; j++) {
                        if (i != j) {
                            double h = weights[i][j]; //0;
                            //if (adjacency[i][j]) {
                            //    h = -log2(s.f(dist(dbRestored[i], dbRestored[j])));
                            //} else {
                            //    h = -log2(1.0 - s.f(dist(dbRestored[i], dbRestored[j])));
                            //}
                            sumWeight += h;
                            for (int jj = 0; jj < d; jj++) {
                                dbRestoredNew[i][jj] += h * (coord[j][jj] + sij[i][j] * (coord[i][jj] - coord[j][jj]));
                            }
                        }
                    }
                    for (int jj = 0; jj < d; jj++) {
                        dbRestoredNew[i][jj] /= sumWeight;
                    }
                }
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < d; j++) {
                        coord[i][j] = dbRestoredNew[i][j];
                    }
                }
            }

            loop++;
        } catch (ConcurrentModificationException cme) {
        }

    }

    public void step() {
        //randomShuffle();

        update();
        //TEST
//        String name = new Integer(currentIteration).toString();
//        IO ea = new IO();
//        ea.writeDoubleToMatlab(weights, "w_" + name);
//        currentIteration++;
        //TEST

    }
//finally {
//		fireStateChanged();
//		}

    public static double log2(double x) {
        final double l2 = 1.0 / Math.log(2.0);
        return Math.log(x) * l2;
    }

    private synchronized void randomShuffle() {
        for (int i = 0; i < n; i++) {
            xydata[i].setLocation(r.nextDouble(), r.nextDouble());
        }
        scaleToDisplaySize();
    }

    public void scaleToDisplaySizeLargest() {
        dd = getSize();
        double height = dd.getHeight();
        double width = dd.getWidth();
        //scale coordinates between 0 and 1
        double max_x = -java.lang.Double.MAX_VALUE;
        double max_y = -java.lang.Double.MAX_VALUE;
        double min_x = -max_x;
        double min_y = -max_y;
        for (int i = 0; i < n; i++) {

            if (xydata[i].getX() > max_x) {
                max_x = xydata[i].getX();
            }
            if (xydata[i].getX() < min_x) {
                min_x = xydata[i].getX();
            }
            if (xydata[i].getY() > max_y) {
                max_y = xydata[i].getY();
            }
            if (xydata[i].getY() < min_y) {
                min_y = xydata[i].getY();
            }
        }
        double max = Math.max(max_x, max_y);
        double min = Math.min(min_x, min_y);
        for (int i = 0; i < xydata.length; i++) {
            double x_new = ((xydata[i].getX() - min) / (max - min)) * 0.99 * width;
            double y_new = ((xydata[i].getY() - min) / (max_y - min)) * 0.99 * height;
            xydata[i].setLocation(x_new, y_new);
        }
    }

    /**
     * Shift all vertices so that the center of gravity is located at the center
     * of the screen.
     */
    public void scaleToDisplaySize() {
        dd = getSize();
        double height = dd.getHeight();
        double width = dd.getWidth();
        //scale coordinates between 0 and 1
        double max_x = -java.lang.Double.MAX_VALUE;
        double max_y = -java.lang.Double.MAX_VALUE;
        double min_x = -max_x;
        double min_y = -max_y;
        for (int i = 0; i < n; i++) {

            if (xydata[i].getX() > max_x) {
                max_x = xydata[i].getX();
            }
            if (xydata[i].getX() < min_x) {
                min_x = xydata[i].getX();
            }
            if (xydata[i].getY() > max_y) {
                max_y = xydata[i].getY();
            }
            if (xydata[i].getY() < min_y) {
                min_y = xydata[i].getY();
            }
        }

        for (int i = 0; i < xydata.length; i++) {
            double x_new = ((xydata[i].getX() - min_x) / (max_x - min_x)) * 0.99 * width + 20;
            double y_new = ((xydata[i].getY() - min_y) / (max_y - min_y)) * 0.99 * height + 20;
            xydata[i].setLocation(x_new, y_new);
        }
    }

    /* (non-Javadoc)
     * @see edu.uci.ics.jung.visualization.layout.AbstractLayout#setSize(java.awt.Dimension)
     */
    @Override
    public void setSize(Dimension size) {
        setInitializer(new RandomLocationTransformer<Integer>(size));
        super.setSize(size);
    }

    /**
     * Enable or disable gravity point adjusting.
     */
    public void setAdjustForGravity(boolean on) {
        adjustForGravity = on;
    }

    /**
     * Returns true if gravity point adjusting is enabled.
     */
    public boolean getAdjustForGravity() {
        return adjustForGravity;
    }

    public void reset() {
        currentIteration = 0;
    }
}
