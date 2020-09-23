/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import Jama.Matrix;
import edu.uci.ics.jung.graph.Graph;
import java.util.Random;
import utils.IO;

/**
 *
 * @author claudia.plant
 */
public class WeightedMajorizationIncremental {

    double[][] coordOld;
    double[][] coord;
    double[][] coordNew;
    double[] sumWeights;
    //double[] sumWeightsNew;
    SimpleSigmoid s;
    Graph g;
    int n;
    int d;

    public WeightedMajorizationIncremental() {
    }

    public WeightedMajorizationIncremental(Graph g) {
        n = g.getVertexCount();
        d = 2;
        this.g = g;
        coord = new double[n][d];
        coordOld = new double[n][d];
        coordNew = new double[n][d];
        sumWeights = new double[n]; //weights used to compute coord from coordOld
        // sumWeightsNew = new double[n]; //weights used to compute coordNew from coord
        s = new SimpleSigmoid(1.0);
        //initialize();
    }

    // initialize coordOld randomly
    //one full round of weighted majorization in order to initialize coord and sumWeights
    public void initialize() {
        Random r = new Random(1);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < d; j++) {
                coordOld[i][j] = r.nextDouble();
            }
        }
        double[][] distances = new double[n][n];
        double[][] weights = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                double[] dw = s.parabola(dist(coordOld[i], coordOld[j]), g.isNeighbor(i, j));
                distances[i][j] = distances[j][i] = dw[0];
                weights[i][j] = weights[j][i] = dw[1];
                sumWeights[i] += dw[1];
                // sumWeightsNew[i] += dw[1];
                sumWeights[j] += dw[1];
                //sumWeightsNew[j] += dw[1];
            }
        }
        WeightedMajorizationPlain wp = new WeightedMajorizationPlain(distances, new Matrix(coordOld).transpose().getArrayCopy(), weights);
        wp.iterate();
        double[][] bla = wp.getPositions();
        coord = new Matrix(bla).transpose().getArrayCopy();
        coordNew = new Matrix(bla).transpose().getArrayCopy();
        IO ea = new IO();
        ea.writeDoubleToMatlab(coord, "initIncremental", "initIncremental");

    }

    //update the coordinates of point i considering the influence of point j on it. writes coordNew[i][], weightsNew[i]
    public void update(int i, int j, boolean connected) {
        //compute influence of j on coord[i][] based on coordOld
        double[] infOld = new double[d];
        double[] dwOld = s.parabola(dist(coordOld[i], coordOld[j]), connected);
        double sijOld = 0;
        for (int k = 0; k < d; k++) {
            sijOld += Math.pow((coordOld[i][k] - coordOld[j][k]), 2);
        }
        //   sij = Math.max(sij, 1e-200);
        if (sijOld != 0) {
            sijOld = dwOld[0] / Math.sqrt(sijOld);
        }
        for (int k = 0; k < d; k++) {
            infOld[k] = dwOld[1] * (coordOld[j][k] + sijOld * (coord[i][k] - coord[j][k]));
            infOld[k] /= sumWeights[i];
//                        if (Double.isNaN(positions[i][k])) {
//                           System.out.println("i " + i + " k: " + k + " j: " + j);
//                        }
        }
        //compute influence of j on coordNew[i][] based on coord, update sumWeightsNew
        double[] infNew = new double[d];
        double[] dw = s.parabola(dist(coord[i], coord[j]), connected);
        double sij = 0;
        for (int k = 0; k < d; k++) {
            sijOld += Math.pow((coord[i][k] - coord[j][k]), 2);
        }
        //   sij = Math.max(sij, 1e-200);
        if (sij != 0) {
            sij = dw[0] / Math.sqrt(sij);
        }
        double sumWeightsUpdated = sumWeights[i] - dwOld[1] + dw[1];
        //sumWeightsNew[i] = sumWeights[i] - 2* dwOld[1] + 2* dw[1];
        for (int k = 0; k < d; k++) {
            infNew[k] = dw[1] * (coord[j][k] + sij * (coord[i][k] - coord[j][k]));
            infNew[k] /= sumWeightsUpdated;
//                        if (Double.isNaN(positions[i][k])) {
//                           System.out.println("i " + i + " k: " + k + " j: " + j);
//                        }
        }
        for (int k = 0; k < d; k++) {
            coordNew[i][k] -= infOld[k];
            coordNew[i][k] += infNew[k];
        }
        sumWeights[i] = sumWeightsUpdated;

    }

    //one full round of updates
    public void updateAll() {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    update(i, j, g.isNeighbor(i, j));
                }
            }
        }
        for(int i = 0; i < n; i++){
            for(int j = 0; j < d; j++)
                coord[i][j] = coordNew[i][j];
        }

    }

    public static double dist(double[] x, double[] y) {
        int d = x.length;
        double result = 0;
        for (int i = 0; i < d; i++) {
            result += (x[i] - y[i]) * (x[i] - y[i]);
        }
        return Math.sqrt(result);
    }
}
