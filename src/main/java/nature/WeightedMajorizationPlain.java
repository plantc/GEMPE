/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import Jama.Matrix;
import edu.uci.ics.jung.graph.Graph;
import java.util.Random;
import static nature.WeightedMajorizationIncremental.dist;
import utils.IO;
/**
 *
 * @author claudia.plant
 */
public class WeightedMajorizationPlain {

    double[][] distances;
    double[][] weights;
    double[][] coord;
    double[][] positions; //updated coords
    int n;
    int d;

    public WeightedMajorizationPlain() {
    }

    public WeightedMajorizationPlain(double[][] distances, double[][] coord, double[][] weights) {
        this.distances = distances;
        this.weights = weights;
        this.coord = new Matrix(coord).transpose().getArrayCopy();
        n = this.coord.length;
        d = this.coord[0].length;
        positions = new double[n][d];
    }

    public double[][] getPositions() {
        return new Matrix(positions).transpose().getArrayCopy();
    }

    public double[][] initAndTwoInterations(Graph g) {
        this.n = g.getVertexCount();
        d = 2;
        Random r = new Random(1);
        double[][] initc = new double[n][d];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < d; j++) {
                initc[i][j] = r.nextDouble();
            }
        }
        SimpleSigmoid s = new SimpleSigmoid(1.0);
        double[][] dd = new double[n][n];
        double[][] ww = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                double[] dw = s.parabola(dist(initc[i], initc[j]), g.isNeighbor(i, j));
                dd[i][j] = dd[j][i] = dw[0];
                ww[i][j] = ww[j][i] = dw[1];

            }
        }
        WeightedMajorizationPlain wp = new WeightedMajorizationPlain(dd, new Matrix(initc).transpose().getArrayCopy(), ww);
        wp.iterate();
        double[][] bla = wp.getPositions();
        double[][] coordFirst = new Matrix(bla).transpose().getArrayCopy();
        IO ea = new IO();
        ea.writeDoubleToMatlab(coordFirst, "init", "initPlain");
        dd = new double[n][n];
        ww = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                double[] dw = s.parabola(dist(coordFirst[i], coordFirst[j]), g.isNeighbor(i, j));
                dd[i][j] = dd[j][i] = dw[0];
                ww[i][j] = ww[j][i] = dw[1];

            }
        }
        WeightedMajorizationPlain wp1 = new WeightedMajorizationPlain(dd, new Matrix(coordFirst).transpose().getArrayCopy(), ww);
        wp1.iterate();
        return new Matrix(wp1.getPositions()).transpose().getArrayCopy();
    }

    //compute coords for all objects
    public void iterate() {
        for (int i = 0; i < n; i++) {
            //DEBUG
//            if (i == 21) {
//                System.out.println("m");
//            }
            //DEBUG
            double sumWeights = 0.0;
            for (int j = 0; j < n; j++) {
                if (i != j) {
//                    if(j == 23)
//                        System.out.println("m");
                    sumWeights += weights[i][j];
                    double sij = 0;
                    for (int k = 0; k < d; k++) {
                        sij += Math.pow((coord[i][k] - coord[j][k]), 2);
                    }
                    //   sij = Math.max(sij, 1e-200);
                    if (sij != 0) {
                        sij = distances[i][j] / Math.sqrt(sij);
                    }

                    for (int k = 0; k < d; k++) {
                        positions[i][k] += weights[i][j] * (coord[j][k] + sij * (coord[i][k] - coord[j][k]));
//                        if (Double.isNaN(positions[i][k])) {
//                           System.out.println("i " + i + " k: " + k + " j: " + j);
//                        }
                    }
                }
            }//j
            for (int k = 0; k < d; k++) {
                positions[i][k] /= sumWeights;
            }
        }//i

    }

    //get new coordinates for a single point; implement this in the embedding display class
    public void iterate(int i) {
        double sumWeights = 0.0;
        for (int j = 0; j < n; j++) {
            if (i != j) {
//                    if(j == 23)
//                        System.out.println("m");
                sumWeights += weights[i][j];
                double sij = 0;
                for (int k = 0; k < d; k++) {
                    sij += Math.pow((coord[i][k] - coord[j][k]), 2);
                }
                //   sij = Math.max(sij, 1e-200);
                if (sij != 0) {
                    sij = distances[i][j] / Math.sqrt(sij);
                }

                for (int k = 0; k < d; k++) {
                    positions[i][k] += weights[i][j] * (coord[j][k] + sij * (coord[i][k] - coord[j][k]));
//                        if (Double.isNaN(positions[i][k])) {
//                           System.out.println("i " + i + " k: " + k + " j: " + j);
//                        }
                }
            }
        }//j
        for (int k = 0; k < d; k++) {
            positions[i][k] /= sumWeights;
        }

    }

}
