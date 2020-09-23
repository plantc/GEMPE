package nature;

public class Sigmoid {

    public double mu ;      // Location Parameter
    public double sigma ;   // Scale Parameter
    public double a ;       // Upper plateau probability
    public double b ;       // Lower plateau probability
    private java.util.Random r ;

    public Sigmoid (double mu, double sigma, double a, double b) {
        this.mu = mu ;
        this.sigma = sigma ;
        this.a = a ;
        this.b = b ;
    }

    public Sigmoid (double edgeprob) {
        mu = sigma = 1.0 ;
        a = b = edgeprob ;
    }

    public Sigmoid (PairOwn[] pairs) {
        int n = pairs.length ;
        java.util.Arrays.sort(pairs) ;
//        //DEBUG
//        double[][] dist = new double[pairs.length][1];
//        double[][] isEdge = new double[pairs.length][1];
//        for(int i = 0; i < pairs.length; i++){
//            dist[i][0] = pairs[i].dist;
//            if(pairs[i].isEdge){
//                isEdge[i][0] = 1.0;
//            }
//        }
//        DataUtils du = new DataUtils();
//        du.saveAsMatlab2(dist, isEdge, "dist", "isEdge", "pairCheck.mat");
//        //DEBUG
        double sum = 0 ;
        double sqrsum = 0 ;
        double weight = 0 ;
        int numEdges = 0 ;
        for (int i=0 ; i<n ; i++)
            if(pairs[i].isEdge)
                numEdges++ ;
        double edgeprob = (double) numEdges / n;
        int start = 0 ;
        int stop = n ;
        while (start < n && !pairs[start].isEdge)
            start ++ ;
        while (stop > start && pairs[stop-1].isEdge)
            stop -- ;
        if (start == stop) {
            a=b=edgeprob ;
            mu=sigma=1.0 ;
            return;
        }
        for (int i = start+1 ; i < stop ; i++) {
            if (pairs[i - 1].isEdge && !pairs[i].isEdge) {
                weight += 1;
                double avgdist = (pairs[i - 1].dist + pairs[i].dist) / 2;
                sum += avgdist;
                sqrsum += avgdist * avgdist;
            } else if (!pairs[i - 1].isEdge && pairs[i].isEdge) {
                weight -= 1;
                double avgdist = (pairs[i - 1].dist + pairs[i].dist) / 2;
                sum -= avgdist;
                sqrsum -= avgdist * avgdist;
            }
        }
        if (weight > 1.1 || weight < 0.9) {
            a=b=edgeprob ;
            mu=sigma=1.0 ;
            return;
        }
        // mu = sum/weight ;
        // sigma = Math.sqrt(sqrsum / weight - sum * sum / weight * weight);
        mu = sum ;
        sigma = Math.sqrt(sqrsum - sum*sum) ;
        // and now determine the scaling parameters a and b of the sigmoid

        double sigavg = 0.0;
        for (int i = 0; i < n ; i++) {
            sigavg += sigmoid(pairs[i].dist, mu, sigma, 1, 0);
        }
        sigavg /= n ;
        double nominator = 0;
        double denominator = 0;
        for (int i = 0; i < n ; i++) {
            double h = sigmoid(pairs[i].dist, sum / weight, sigma, 1, 0) - sigavg;
            nominator += h * ((pairs[i].isEdge ? 1 : 0) - edgeprob);
            denominator += h * h;
        }
        b = Math.max (edgeprob - nominator / denominator * sigavg, 0.0001) ;
        a = Math.min (nominator / denominator + b, 0.9999) ;
        if (a<b)
            a=b=edgeprob ;
    }

    public PairOwn generateRandom(java.util.Random r, double minx, double maxx) {
        double h = r.nextDouble()*(maxx-minx)+minx ;
        return new PairOwn(h, f(h) >= r.nextDouble()) ;
    }

    public PairOwn generateRandom (double minx, double maxx) {
        if (r==null)
            r = new java.util.Random(10) ;
        return generateRandom (r, minx, maxx) ;
    }

    public double f(double x) {
        return sigmoid(x,mu,sigma,a,b) ;
    }

    public double diffF(double x) {
        // Gaussian, scaled by (a-b)
        double h = (x-mu)/sigma ;
        return Math.exp(-0.5*h*h)/Math.sqrt(2*Math.PI*sigma*sigma)*(a-b) ;
    }

    public double costEdge(double x) {
        return -Math.log(sigmoid (x,mu,sigma,a,b) ) / Math.log(2.0) ;
    }

    public double costNoEdge(double x) {
        return -Math.log(1.0-sigmoid (x,mu,sigma,a,b) ) / Math.log(2.0) ;
    }

    public double diffCostEdge(double x) {
        // differentiation of the cost for an actual edge
        // x is the distance of the embedding points
        return (a-b)*Math.exp(-Math.pow(mu-x,2.0)/(sigma*sigma)/2)*Math.sqrt(2.0)/Math.sqrt(
                0.3141592653589793E1)/sigma/(b+a+a*normp(Math.sqrt(2.0)*(mu-x)/sigma/2)-b*normp(Math.sqrt(2.0
        )*(mu-x)/sigma/2))/Math.log(2.0);
    }

    public double diffCostNoEdge(double x) {
        // differentiation of the cost of a node pair with no edge between
        return (a-b)*Math.exp(-Math.pow(mu-x,2.0)/(sigma*sigma)/2)*Math.sqrt(2.0)/Math.sqrt(
                0.3141592653589793E1)/sigma/(-b+a*normp(Math.sqrt(2.0)*(mu-x)/sigma/2)-a-b*normp(Math.sqrt(
                2.0)*(mu-x)/sigma/2))/Math.log(2.0);
    }

    private static double sigmoid(double x, double mu, double sigma, double a, double b) {
        if (a==b)
            return a ;
        if (sigma <= 0)
            return x<mu ? a : b ;
        return (a - b) * (1 - normp((x - mu) / sigma)) + b;
    }

    public String toString() {
        return "(mu="+mu+"; sigma="+sigma+"; a="+a+"; b="+b+")" ;
    }

    private static double normp(double z) {
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
            p = expntl * ((((((p6 * zabs + p5) * zabs + p4) * zabs + p3) * zabs
                    + p2) * zabs + p1) * zabs + p0) / (((((((q7 * zabs + q6) * zabs
                    + q5) * zabs + q4) * zabs + q3) * zabs + q2) * zabs + q1) * zabs
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

}

