/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import static java.lang.Math.*;
import org.apache.commons.math3.special.Erf;

/**
 *
 * @author boehm
 */
public class WeightedSigmoid {

    public double sigma;
    private static final double isq2 = 0.7071067814, sq2ipi = 0.7978845605, is2spil2 = 0.5755520493, ln2 = 0.6931471806;
    private static final double minWeight = 1e-9, maxWeight = 2., minSigma = 1E-6;
    //isq2pi2 = 0.1994711401, isqpi = 0.5641895835,
    public final double mu = 1.5;

    public WeightedSigmoid(double s) {
        sigma = s;
    }

    //15.5.15: determine also mu for the comparison methods
    public WeightedSigmoid(PairOwn[] pair, boolean comparison) {
        //get the largest distance of edge as initial guess
        double maxDistEdge = -Double.MAX_VALUE;
        double meanDistEdge = 0.0;
        int edgeCount = 0;
        for (int i = 0; i < pair.length; i++) {
            if (pair[i].isEdge) {
                meanDistEdge += pair[i].dist;
                edgeCount++;
                if (pair[i].dist < maxDistEdge) {
                    maxDistEdge = pair[i].dist;
                }
            }
        }
        meanDistEdge /= edgeCount;
        double muMin = meanDistEdge;
        double muMax = 10 * maxDistEdge;
        int numSteps = (int) Math.floor((muMax - muMin) / 0.1);
        double[] costs = new double[numSteps];

        for (int s = 0; s < numSteps; s++) {
            double aktMu = muMin + numSteps * 0.1;
            PairOwn[] p = pair.clone();
            for (int i = 0; i < p.length; i++) {
                p[i] = new PairOwn((p[i].dist - aktMu) * isq2, p[i].isEdge);
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

        }
    }

    public WeightedSigmoid(WeightedPair[] pair) {
        WeightedPair[] p = pair.clone();
        for (int i = 0; i < p.length; i++) {
            p[i] = new WeightedPair((p[i].dist - mu) * isq2, p[i].alpha, p[i].beta);
        }
        double upper = 20.;
        double lower = 0.;
        for (int i = 0; i < 30; i++) {
            sigma = (upper + lower) / 2.;
            double sum1 = 0;
            // double bla = 0.0; //TEST
            int count = 0;
            for (WeightedPair x : p) {
                double h = x.dist / sigma;
                //    sum1 += x.alpha * (x.dist * exp(-h * h) / (1. - erf(h)));
                //    sum1 -= x.beta * (x.dist * exp(-h * h) / (1. + erf(h)));
                //TEST
                if (x.alpha != 0) {
                    sum1 += x.alpha * (x.dist * exp(-h * h) / (1. - erf(h)));
                }
                if (x.beta != 0) {
                    sum1 -= x.beta * (x.dist * exp(-h * h) / (1. + erf(h)));
                }

//                    if(x.alpha == 1 && x.beta == 0){
//                        bla += x.dist * exp(-h * h) / (1. - erf(h));
//                    }
//                    if(x.beta == 1 && x.alpha == 0){
//                        bla -= x.dist * exp(-h * h) / (1. + erf(h));
//
//                    }
                count++;
//                    if(Double.isNaN(sum1))
//                        System.out.println(count);
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
            for (WeightedPair x : p) {
                double h = x.dist / sigma;
                double hh = x.alpha * (x.dist * exp(-h * h) / (1. - erf(h)));
                sum1 += hh;
                sum2 += hh * (x.dist * x.dist * s3inv - x.dist * s2invc);
                hh = x.beta * (x.dist * exp(-h * h) / (1. + erf(h)));
                sum1 -= hh;
                sum2 -= hh * (x.dist * x.dist * s3inv + x.dist * s2invc);
            }
            sigma -= sum1 / sum2;
            if (sigma < minSigma || Double.isNaN(sigma)) {
                sigma = minSigma;
                break;
            }
            // System.out.println(sigma+", "+log(abs(sum1))) ;
        }
    }

    public static double[] dwFormula(double delta, double mu, double sigma, double alpha, double beta) {
        // Parameters:
        // delta: Euclidean distance between node embeddings
        // mu, sigma: parameters of the sigmoid
        // alpha: number of edges between nodes
        // beta: number of non-edges between nodes
        // return: Parabola array where return[0] is aimed distance and return[1] is weight.
        // sanity checks: alpha=1, beta=0 equals old edge-parabola
        //                alpha=0, beta=1 equals old non-edge-parabola
        // return[0] is invariant to uniform scaling of alpha and beta
        // return[1] is not
        final double sqrt2 = Math.sqrt(2);
        final double sqrt2pi = Math.PI * sqrt2;
        final double sqrtpi = Math.sqrt(Math.PI);
        double delmu = delta - mu;
        double ERF = erf(delmu / sqrt2 / sigma);
        double EXP = Math.exp(-delmu * delmu / (2 * sigma * sigma));
        double COMMON = (Math.sqrt(0.2e1) * Math.PI * alpha * delta
                - Math.sqrt(0.2e1) * Math.PI * alpha * mu + Math.sqrt(0.2e1) * Math.PI * beta * delta - Math.sqrt(0.2e1) * Math.PI * beta * mu) * Math.pow(ERF,
                0.3e1) + ((0.2e1 * Math.sqrt(Math.PI) * alpha * sigma + 0.2e1
                * Math.sqrt(Math.PI) * beta * sigma) * EXP + Math.sqrt(0.2e1) * Math.PI * alpha * delta - Math.sqrt(0.2e1) * Math.PI * alpha * mu
                - Math.sqrt(0.2e1) * Math.PI * beta * delta + Math.sqrt(0.2e1) * Math.PI * beta * mu) * ERF * ERF + ((0.4e1 * Math.sqrt(Math.PI) * alpha * sigma
                - 0.4e1 * Math.sqrt(Math.PI) * beta * sigma) * EXP - Math.sqrt(0.2e1) * Math.PI * alpha * delta + Math.sqrt(0.2e1) * Math.PI * alpha * mu
                - Math.sqrt(0.2e1) * Math.PI * beta * delta + Math.sqrt(0.2e1) * Math.PI * beta * mu) * ERF + (0.2e1 * Math.sqrt(Math.PI) * alpha * sigma + 0.2e1
                * Math.sqrt(Math.PI) * beta * sigma) * EXP - Math.sqrt(0.2e1) * Math.PI * alpha * delta + Math.sqrt(0.2e1) * Math.PI * alpha * mu
                + Math.sqrt(0.2e1) * Math.PI * beta * delta - Math.sqrt(0.2e1) * Math.PI * beta * mu;

        double w = EXP * COMMON * Math.pow(Math.PI, -0.3e1 / 0.2e1) * Math.pow(sigma, -0.3e1) / Math.log(0.2e1) / (Math.pow(ERF, 0.4e1)
                - 0.2e1 * ERF * ERF + 0.1e1) / 0.2e1;
//        if (Double.isInfinite(w))
//            w = 0.001 ;
//        if (w < 1e-9 || Double.isNaN(w) || Double.isInfinite(w)) {
//            return new double[]{delta, 0};
//        }

        double d = ((Math.sqrt(0.2e1) * Math.PI * alpha * delta * delta
                - Math.sqrt(0.2e1) * Math.PI * alpha * delta * mu + alpha
                * Math.sqrt(0.2e1) * Math.PI * sigma * sigma + Math.sqrt(0.2e1) * Math.PI
                * beta * delta * delta - Math.sqrt(0.2e1) * Math.PI * beta * delta * mu
                + beta * Math.sqrt(0.2e1) * Math.PI * sigma * sigma) * Math.pow(ERF,
                0.3e1) + ((0.2e1 * Math.sqrt(Math.PI) * alpha * delta * sigma + 0.2e1
                * Math.sqrt(Math.PI) * beta * delta * sigma) * EXP + Math.sqrt(0.2e1) * Math.PI * alpha * delta * delta - Math.sqrt(0.2e1) * Math.PI * alpha * delta * mu + alpha * Math.sqrt(0.2e1) * Math.PI * sigma * sigma
                - Math.sqrt(0.2e1) * Math.PI * beta * delta * delta + Math.sqrt(0.2e1) * Math.PI * beta * delta * mu - beta * Math.sqrt(0.2e1) * Math.PI * sigma
                * sigma) * ERF * ERF + ((0.4e1 * Math.sqrt(Math.PI) * alpha * delta * sigma - 0.4e1 * Math.sqrt(Math.PI) * beta * delta * sigma) * EXP
                - Math.sqrt(0.2e1) * Math.PI * alpha * delta * delta + Math.sqrt(0.2e1) * Math.PI * alpha * delta * mu - alpha * Math.sqrt(0.2e1) * Math.PI * sigma * sigma - Math.sqrt(0.2e1) * Math.PI * beta * delta * delta
                + Math.sqrt(0.2e1) * Math.PI * beta * delta * mu - beta * Math.sqrt(0.2e1)
                * Math.PI * sigma * sigma) * ERF + (0.2e1 * Math.sqrt(Math.PI) * alpha * delta * sigma + 0.2e1 * Math.sqrt(Math.PI) * beta * delta * sigma) * EXP
                - Math.sqrt(0.2e1) * Math.PI * alpha * delta * delta + Math.sqrt(0.2e1)
                * Math.PI * alpha * delta * mu - alpha * Math.sqrt(0.2e1) * Math.PI * sigma * sigma + Math.sqrt(0.2e1) * Math.PI * beta * delta * delta
                - Math.sqrt(0.2e1) * Math.PI * beta * delta * mu + beta * Math.sqrt(0.2e1)
                * Math.PI * sigma * sigma) / COMMON;

        boolean caseCommonZero = false;
        boolean caseErf4Zero = false;
        boolean caseOutside = false;

        if(COMMON == 0)
            caseCommonZero = true;

        double term1 = Math.pow(ERF, 0.4e1)- 0.2e1 * ERF * ERF + 0.1e1;
        if(term1 == 0)
            caseErf4Zero = true;

        double term2 = (delta-mu)/sigma;

        if(term2 > 5)
            caseOutside = true;

        if(caseCommonZero || caseErf4Zero|| caseOutside){
            w = alpha/(2*Math.log(2)*Math.pow(sigma, 2.0));
            if(alpha != 0){
                d = mu;
            }
            else{
                d = delta;
            }

        }

        double[] res = new double[2];
        res[0] = d;
        res[1] = w;

//        return new double[]{
//            ((Math.sqrt(0.2e1) * Math.PI * alpha * delta * delta
//            - Math.sqrt(0.2e1) * Math.PI * alpha * delta * mu + alpha
//            * Math.sqrt(0.2e1) * Math.PI * sigma * sigma + Math.sqrt(0.2e1) * Math.PI
//            * beta * delta * delta - Math.sqrt(0.2e1) * Math.PI * beta * delta * mu
//            + beta * Math.sqrt(0.2e1) * Math.PI * sigma * sigma) * Math.pow(ERF,
//            0.3e1) + ((0.2e1 * Math.sqrt(Math.PI) * alpha * delta * sigma + 0.2e1
//            * Math.sqrt(Math.PI) * beta * delta * sigma) * EXP + Math.sqrt(0.2e1) * Math.PI * alpha * delta * delta - Math.sqrt(0.2e1) * Math.PI * alpha * delta * mu + alpha * Math.sqrt(0.2e1) * Math.PI * sigma * sigma
//            - Math.sqrt(0.2e1) * Math.PI * beta * delta * delta + Math.sqrt(0.2e1) * Math.PI * beta * delta * mu - beta * Math.sqrt(0.2e1) * Math.PI * sigma
//            * sigma) * ERF * ERF + ((0.4e1 * Math.sqrt(Math.PI) * alpha * delta * sigma - 0.4e1 * Math.sqrt(Math.PI) * beta * delta * sigma) * EXP
//            - Math.sqrt(0.2e1) * Math.PI * alpha * delta * delta + Math.sqrt(0.2e1) * Math.PI * alpha * delta * mu - alpha * Math.sqrt(0.2e1) * Math.PI * sigma * sigma - Math.sqrt(0.2e1) * Math.PI * beta * delta * delta
//            + Math.sqrt(0.2e1) * Math.PI * beta * delta * mu - beta * Math.sqrt(0.2e1)
//            * Math.PI * sigma * sigma) * ERF + (0.2e1 * Math.sqrt(Math.PI) * alpha * delta * sigma + 0.2e1 * Math.sqrt(Math.PI) * beta * delta * sigma) * EXP
//            - Math.sqrt(0.2e1) * Math.PI * alpha * delta * delta + Math.sqrt(0.2e1)
//            * Math.PI * alpha * delta * mu - alpha * Math.sqrt(0.2e1) * Math.PI * sigma * sigma + Math.sqrt(0.2e1) * Math.PI * beta * delta * delta
//            - Math.sqrt(0.2e1) * Math.PI * beta * delta * mu + beta * Math.sqrt(0.2e1)
//            * Math.PI * sigma * sigma) / COMMON,
//            w
//        };
        return res;

    }

    private double log2(double d){
        return Math.log(d)/Math.log(2);
    }

    public static double[] dwFormulaOld(double delta, double mu, double sigma, double alpha, double beta) {
        // Parameters:
        // delta: Euclidean distance between node embeddings
        // mu, sigma: parameters of the sigmoid
        // alpha: number of edges between nodes
        // beta: number of non-edges between nodes
        // return: Parabola array where return[0] is aimed distance and return[1] is weight.
        // sanity checks: alpha=1, beta=0 equals old edge-parabola
        //                alpha=0, beta=1 equals old non-edge-parabola
        // return[0] is invariant to uniform scaling of alpha and beta
        // return[1] is not
        final double sqrt2 = Math.sqrt(2);
        final double sqrt2pi = Math.PI * sqrt2;
        final double sqrtpi = Math.sqrt(Math.PI);
        double delmu = delta - mu;
        double r = erf(delmu / sqrt2 / sigma);
        double x = Math.exp(delmu * delmu / (2 * sigma * sigma));
        double common
                = (sqrt2pi * r * delmu * (r * r - 1) + 2 * sqrtpi * x * sigma * (r * r + 1)) * (alpha + beta)
                + (sqrt2pi * delmu * (r * r - 1) + 4 * sqrtpi * x * r * sigma) * (alpha - beta);
        double w = x * common / (2 * sqrtpi * Math.PI * sigma * sigma * sigma * Math.log(2) * (r * r - 1) * (r * r - 1));
        double d = ((sqrt2pi * r * (delmu * delmu + sigma * sigma) * (r * r - 1) + 2 * sqrtpi * x * delmu * sigma * (r * r + 1)) * (alpha + beta)
                + (sqrt2pi * (delmu * delmu + sigma * sigma) * (r * r - 1) + 4 * sqrtpi * r * x * delmu * sigma) * (alpha - beta)) / common;
        if (w <= 0) {
            d = delta;
        }
        return new double[]{d, w};
//            ((sqrt2pi * r * (delmu * delmu + sigma * sigma) * (r * r - 1) + 2 * sqrtpi * x * delmu * sigma * (r * r + 1)) * (alpha + beta)
//            + (sqrt2pi * (delmu * delmu + sigma * sigma) * (r * r - 1) + 4 * sqrtpi * r * x * delmu * sigma) * (alpha - beta)) / common, x * common
//            / (2 * sqrtpi * Math.PI * sigma * sigma * sigma * Math.log(2) * (r * r - 1) * (r * r - 1))};
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

    public double costAllPairs(WeightedPair[] p) {
        double sum = 0.;
        for (int i = 0; i < p.length; i++) //DEBUG
        {

            sum += p[i].alpha * costEdge(p[i].dist);

            sum += p[i].beta * costNoEdge(p[i].dist);

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
