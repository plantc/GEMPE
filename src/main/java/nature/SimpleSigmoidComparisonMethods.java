/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import static java.lang.Math.*;
import static nature.SimpleSigmoid.erf;
import org.apache.commons.math3.special.Erf;

/**
 *
 * @author boehm
 */
public class SimpleSigmoidComparisonMethods {

    public double sigma;
    private static final double isq2 = 0.7071067814, sq2ipi = 0.7978845605, is2spil2 = 0.5755520493, ln2 = 0.6931471806;
    private static final double minWeight = 1e-9, maxWeight = 2., minSigma = 1E-6;
    //isq2pi2 = 0.1994711401, isqpi = 0.5641895835,
    public double mu;

    public SimpleSigmoidComparisonMethods(double s) {
        sigma = s;
    }

    //15.5.15: determine also mu for the comparison methods
    public SimpleSigmoidComparisonMethods(PairOwn[] pair) {
        //get the largest distance of edge as initial guess
        double maxDistEdge = -Double.MAX_VALUE;
        double meanDistEdge = 0.0;
        int edgeCount = 0;
        for (int i = 0; i < pair.length; i++) {
            if (pair[i].isEdge) {
                meanDistEdge += pair[i].dist;
                edgeCount++;
                if (pair[i].dist > maxDistEdge) {
                    maxDistEdge = pair[i].dist;
                }
            }
        }
        meanDistEdge /= edgeCount;
        double muMin = meanDistEdge / 3;
        double muMax = 1.0 * maxDistEdge;
        double stepSize = (muMax - muMin) / 3;
        int numSteps = (int) Math.floor((muMax - muMin) / stepSize);

        double[] costs = new double[numSteps];
        double bestMu = 0.0;
        double bestSigma = 0.0;
        double bestCost = Double.MAX_VALUE;
        System.out.println("numSteps " + numSteps);

        for (int s = 0; s < numSteps; s++) {
            mu = muMin + s * stepSize;
            //mu = 1.5;
            PairOwn[] p = pair.clone();
            for (int i = 0; i < p.length; i++) {
                p[i] = new PairOwn((p[i].dist - mu) * isq2, p[i].isEdge);
            }

            double upper = 20.;
            double lower = 0.;
            for (int i = 0; i < 30; i++) {
                sigma = (upper + lower) / 2.;
                double sum1 = 0;
                for (PairOwn x : p) {
                    double h = x.dist / sigma;
                    if (x.isEdge) {
                        sum1 += x.dist * exp(-h * h) / (1. - erf(h));
                    } else {
                        sum1 -= x.dist * exp(-h * h) / (1. + erf(h));
                    }
                }
                // System.out.println(sigma+", "+log(abs(sum1))) ;
                if (sum1 < 0) {
                    upper = sigma;
                } else {
                    lower = sigma;
                }
            }
            sigma = (upper + lower) / 2.;
            // System.out.println(sigma) ;
            for (int i = 0; i < 0; i++) {
                double s2invc = 1. / sigma / sigma;
                double s3inv = s2invc / sigma;
                s2invc *= sq2ipi;
                double sum1 = 0;
                double sum2 = 0;
                for (PairOwn x : p) {
                    double h = x.dist / sigma;
                    if (x.isEdge) {
                        h = x.dist * exp(-h * h) / (1. - erf(h));
                        sum1 += h;
                        sum2 += h * (x.dist * x.dist * s3inv - x.dist * s2invc);
                    } else {
                        h = x.dist * exp(-h * h) / (1. + erf(h));
                        sum1 -= h;
                        sum2 -= h * (x.dist * x.dist * s3inv + x.dist * s2invc);
                    }
                }
                sigma -= sum1 / sum2;
                if (sigma < minSigma || Double.isNaN(sigma)) {
                    sigma = minSigma;
                    break;
                }
                // System.out.println(sigma+", "+log(abs(sum1))) ;
            }
            costs[s] = costAllPairs(pair);
            //System.out.println(costs[s] + " " + sigma + " " + mu);
            if (costs[s] < bestCost) {
                bestCost = costs[s];
                bestMu = mu;
                bestSigma = sigma;
            }
        }
        mu = bestMu;
        sigma = bestSigma;
        System.out.println("bestCost: " + bestCost + " mu: " + bestMu + " sigma: " + bestSigma);
    }

    public double[] parabola(double x, boolean isEdge) {
        double h = (x - mu) * isq2 / sigma;
        double EXP = exp(-h * h);
        double ERF = erf(h);
        if (isEdge) {
            ERF -= 1.;
            ERF = min(ERF, -1e-20);
        } else {
            ERF += 1.;
            ERF = max(ERF, 1e-20);
        }

        double w = ((is2spil2 * x - 0.8633280737) / sigma
                + 0.4592240941 * EXP / ERF)
                * EXP / sigma / sigma / ERF;
//        w = min(w, maxWeight) ;
//        w = max(w, minWeight) ;
        double d = x;
        if (w != 0) {
            d += is2spil2 * EXP / w / sigma / ERF;
        }
//        if (d < 0 || w < 0) {
//            d = 0;
//            w = -is2spil2 * EXP / sigma / x / ERF;
//            w = min(w, maxWeight) ;
//        }
        return new double[]{d, w};
    }

    public double costEdge(double dist) {
//        double h = (dist - mu) * isq2 / sigma;
//        return -log(0.5 - 0.5 * erf(h)) / ln2;
        double h = 0.5 - 0.5 * erf((dist - mu) * isq2 / sigma);
        if (h < 1e-20) {
            h = 1e-20;
        }
        return -log(h) / ln2;
    }

    public double costNoEdge(double dist) {
//        double h = (dist - mu) * isq2 / sigma;
//        return -log(0.5 + 0.5 * erf(h)) / ln2;
        double h = 0.5 + 0.5 * erf((dist - mu) * isq2 / sigma);
        if (h < 1e-20) {
            h = 1e-20;
        }
        return -log(h) / ln2;
    }

    public double costAllPairs(PairOwn[] p) {
        double sum = 0.;
        for (int i = 0; i < p.length; i++) //DEBUG
        {
            if (p[i].isEdge) {
                sum += costEdge(p[i].dist);
            } else {
                sum += costNoEdge(p[i].dist);
            }
            if (Double.isNaN(sum)) {
                // System.out.println("m");
            }
        }

        return sum;

    }

    public static double erf(double x) {
        return Erf.erf(x);
//        if (x >= 0) {
//            double t = 1. / (1. + .5 * x);
//            double t2 = t * t, t3 = t2 * t, t4 = t3 * t, t5 = t4 * t,
//                    t6 = t5 * t, t7 = t6 * t, t8 = t7 * t, t9 = t8 * t;
//            return 1. - t * exp(-x * x - 1.26551223 + 1.00002368 * t
//                    + 0.37409196 * t2
//                    + 0.09678418 * t3 - 0.18628806 * t4 + 0.27886807
//                    * t5 - 1.13520398 * t6
//                    + 1.48851587 * t7 - 0.82215223 * t8 + 0.17087277 * t9);
//        } else {
//            double t = 1. / (1. - .5 * x);
//            double t2 = t * t, t3 = t2 * t, t4 = t3 * t, t5 = t4 * t,
//                    t6 = t5 * t, t7 = t6 * t, t8 = t7 * t, t9 = t8 * t;
//            return t * exp(-x * x - 1.26551223 + 1.00002368 * t
//                    + 0.37409196 * t2
//                    + 0.09678418 * t3 - 0.18628806 * t4 + 0.27886807
//                    * t5 - 1.13520398 * t6
//                    + 1.48851587 * t7 - 0.82215223 * t8 + 0.17087277
//                    * t9) - 1;
//        }
    }
}
