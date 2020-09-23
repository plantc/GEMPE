/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.util.Pair;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import utils.DataUtils;
import visu.Visualization;

/**
 *
 *
 */
public class Parallel {

    Graph g;
    SimpleSigmoid s;
    double[][] distances; //TO DO: save only upper triangle matrix
    double[][] weights;
    double[][] coord;
    int[][][] pId;
    CoordUpdateResult[] cr;
    DWUpdateResult[] dr;
    double sumWeights;
    double cutOff; //grid
    int numBinsX; //grid
    int numBinsY; //grid
    double[][] minMax; //grid
    static double tau = 1E-12;
    static double mu = 1.5;
    static double fixSigma = 1.0;
    int n;
    int m;
    Random r;
    int numT;
    int iter; //iteration counter
    int d; //dimensionality
    static int iterSigmoid = 10; //after how many iterations the sigmoid is recomputed - 100 is the default value
    static int localMinEscape = 3; //number of tries to escape a local cost minimum
    static int numInits = 10; //number of inits
    double bestCost; //written in initializationPhase, stabilizationPhase
    static boolean verbose = true;
    static boolean save = true;
    boolean saveIntermediateResults = true;
    boolean display = false;
    boolean useGrid;
    boolean improvementInCurrentPhase;
    static double convConst = 0.1;
    static int maxIteration = 100000;

    //7.8.2017: this constructor should be used
    public Parallel(Graph g, int numT, int d) {
        this.g = g;
        this.numT = numT;
        this.d = d;
        n = g.getVertexCount();
        m = n * (n - 1) / 2;
        //r = new Random(1);
        iter = 0;
        coord = new double[n][d];
        cr = new CoordUpdateResult[numT];
        dr = new DWUpdateResult[numT];
//        for (int i = 0; i < coord.length; i++) {
//            for (int j = 0; j < d; j++) {
//                coord[i][j] = r.nextDouble();
//            }
//        }
//        distances = new double[n][n];
//        for (int i = 0; i < n; i++) {
//            for (int j = 0; j < i; j++) {
//                distances[i][j] = distances[j][i] = dist(coord[i], coord[j]);
//            }
//        }
        weights = new double[n][n];
        bestCost = Double.MAX_VALUE;
        //useGrid = true;
//        numBinsX = 0;
//        numBinsY = 0;
        // computeSigmoid();
    }

    //7.8.2017: 3D without grid
    public void threeD() {
        double[][] init = getBestInitialization();
        DataUtils du = new DataUtils();
        Visualization v = new Visualization(g);
        if (save) {
            //du.saveAsMatlab(init, "init", "init.mat");
            du.saveResult(init, "init.txt");
        }
        if (display) {
            v.displayCoordNew(init, " ");
        }
        setCoord(init);
        double[][] res = finishPhase();
        if (save) {
            du.saveResult(res, "res.txt");
            //du.saveAsMatlab(res, "res", "res.mat");
        }
        if (display) {
            v.displayCoordNew(res, " ");
        }
    }

    //26.10.17: run only the finishing phase
    public void twoDFinishOnly(double[][] init) {
        setCoord(init);
        startWithGrid();
        DataUtils du = new DataUtils();
        double[][] res = finishPhase();
        if (improvementInCurrentPhase) {
            if (save) {
                du.saveResult(res, "res.txt");
                //du.saveAsMatlab(res, "res", "res.mat");
            }

        } else {
            System.out.println("No improvement in polishing phase.");
        }
    }

    //8.8.2017: 2D with grid
    public void twoD() {
        double[][] init = getBestInitialization();
        DataUtils du = new DataUtils();
        Visualization v = new Visualization(g);
        if (save) {
            //du.saveAsMatlab(init, "init", "init.mat");
            du.saveResult(init, "init.txt");

        }
        if (display) {
            v.displayCoordNew(init, " ");
        }
        //System.out.println(init[0][0] + " " + init[0][1]);
        double[][] initC = init.clone();
        for (int i = 0; i < initC.length; i++) {
            initC[i] = init[i].clone();
        }
        setCoord(initC);
        //  System.out.println(init[0][0] + " " + init[0][1]); - ok
        startWithGrid();
        //System.out.println(init[0][0] + " " + init[0][1]); //- ok
        double[][] stab = stabilizationPhase();
        if (improvementInCurrentPhase) {
            if (save) {
                //du.saveAsMatlab(stab, "stab", "stab.mat");
                du.saveResult(stab, "stab.txt");
            }
            if (display) {
                v.displayCoordNew(stab, " ");
            }
            double[][] stabC = stab.clone();
            for (int i = 0; i < stabC.length; i++) {
                stabC[i] = stabC[i].clone();
            }
            setCoord(stabC);
        } else {
            System.out.println("No improvement in stabilization phase.");
            //System.out.println(init[0][0] + " " + init[0][1]); // modified
            initC = init.clone();
            for (int i = 0; i < initC.length; i++) {
                initC[i] = init[i].clone();
            }
            setCoord(initC);
        }

        startWithGrid();
        double[][] res = finishPhase();
        if (improvementInCurrentPhase) {
            if (save) {
                //du.saveAsMatlab(res, "res", "res.mat");
                du.saveResult(res, "res.txt");
            }
            if (display) {
                v.displayCoordNew(res, " ");
            }
        } else {
            System.out.println("No improvement in polishing phase.");
        }

    }

    //8.8.2017: 2D with grid
    public void twoDGridStart() {
        double[][] init = getBestInitializationGrid();
        DataUtils du = new DataUtils();
        Visualization v = new Visualization(g);
        if (save) {
            //du.saveAsMatlab(init, "init", "init.mat");
            du.saveResult(init, "init.txt");
        }
        if (display) {
            v.displayCoordNew(init, "init");
        }
        //System.out.println(init[0][0] + " " + init[0][1]);
        double[][] initC = init.clone();
        for (int i = 0; i < initC.length; i++) {
            initC[i] = init[i].clone();
        }
        setCoord(initC);
        //  System.out.println(init[0][0] + " " + init[0][1]); - ok
        startWithGrid();
        //System.out.println(init[0][0] + " " + init[0][1]); //- ok
//        double[][] stab = stabilizationPhase();
//        if (improvementInCurrentPhase) {
//            if (save) {
//                du.saveAsMatlab(stab, "stab", "stab.mat");
//            }
//            if (display) {
//                v.displayCoordNew(stab, " ");
//            }
//            double[][] stabC = stab.clone();
//            for (int i = 0; i < stabC.length; i++) {
//                stabC[i] = stabC[i].clone();
//            }
//            setCoord(stabC);
//        } else {
//            System.out.println("No improvement in stabilization phase.");
//            //System.out.println(init[0][0] + " " + init[0][1]); // modified
//            initC = init.clone();
//            for (int i = 0; i < initC.length; i++) {
//                initC[i] = init[i].clone();
//            }
//            setCoord(initC);
//        }
//
//        startWithGrid();
        double[][] res = finishPhase();
        if (improvementInCurrentPhase) {
            if (save) {
                du.saveResult(res, "res.txt");
                //du.saveAsMatlab(res, "res", "res.mat");
            }
            if (display) {
                v.displayCoordNew(res, "res");
            }
        } else {
            System.out.println("No improvement in polishing phase.");
        }

    }

    public Parallel(Graph g, int numT, double[][] coord, boolean useGrid) {
        this.g = g;
        this.numT = numT;
        n = g.getVertexCount();
        m = n * (n - 1) / 2;
        //r = new Random(1);
        iter = 0;
        this.coord = coord;
        this.useGrid = useGrid;
        d = coord[0].length;
        cr = new CoordUpdateResult[numT];
        dr = new DWUpdateResult[numT];
//        for (int i = 0; i < coord.length; i++) {
//            for (int j = 0; j < d; j++) {
//                coord[i][j] = r.nextDouble();
//            }
//        }
//        distances = new double[n][n];
//        for (int i = 0; i < n; i++) {
//            for (int j = 0; j < i; j++) {
//                distances[i][j] = distances[j][i] = dist(coord[i], coord[j]);
//            }
//        }
        weights = new double[n][n];
        //useGrid = true;
//        numBinsX = 0;
//        numBinsY = 0;
        computeSigmoid(true);
    }

    public void setCoord(double[][] coord) {
        this.coord = coord;
        iter = 0;
        cr = new CoordUpdateResult[numT];
        dr = new DWUpdateResult[numT];
        weights = new double[n][n];
    }

    public void setUseGrid(boolean useGrid) {
        this.useGrid = useGrid;
    }

    public void startWithGrid() {
        iter = 0;
        cr = new CoordUpdateResult[numT];
        dr = new DWUpdateResult[numT];
        weights = new double[n][n];
        useGrid = true;
        computeSigmoid(true);
        //setUpGrid();
    }

    public void reset() {
        iter = 0;
        coord = new double[n][d];
        cr = new CoordUpdateResult[numT];
        dr = new DWUpdateResult[numT];
        weights = new double[n][n];
        bestCost = Double.MAX_VALUE;

    }

    public void initializeCoords(int seed) {
        r = new Random(seed);
        for (int i = 0; i < coord.length; i++) {
            for (int j = 0; j < d; j++) {
                coord[i][j] = r.nextDouble();
            }
        }
        distances = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                distances[i][j] = distances[j][i] = dist(coord[i], coord[j]);
            }
        }
        s = new SimpleSigmoid(fixSigma);
        computeCost();
    }

    public void initializeCoordsFreeSigma(int seed) {
        r = new Random(seed);
        for (int i = 0; i < coord.length; i++) {
            for (int j = 0; j < d; j++) {
                coord[i][j] = r.nextDouble();
            }
        }
        PairOwn[] p = new PairOwn[m];
        distances = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                distances[i][j] = distances[j][i] = dist(coord[i], coord[j]);
                p[i * (i - 1) / 2 + j] = new PairOwn(dist(coord[i], coord[j]), g.isNeighbor(i, j));
            }
        }
        s = new SimpleSigmoid(p);
        double aktCost = s.costAllPairs(p);
        bestCost = aktCost;
        setUpGrid();

    }

    private double computeCost() {
        //compute sigmoid
        PairOwn[] p = new PairOwn[m];
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < i; j++) {
                p[i * (i - 1) / 2 + j] = new PairOwn(dist(coord[i], coord[j]), g.isNeighbor(i, j));
            }
        }
        double aktCost = s.costAllPairs(p);
        if (verbose) {
            System.out.println(iter + " " + s.sigma + " " + aktCost);
        }
        return aktCost;

    }

    private double computeSigmoid(boolean initCost) {
        //compute sigmoid
        PairOwn[] p = new PairOwn[m];
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < i; j++) {
                p[i * (i - 1) / 2 + j] = new PairOwn(dist(coord[i], coord[j]), g.isNeighbor(i, j));
            }
        }
        s = new SimpleSigmoid(p);
        double aktCost = s.costAllPairs(p);
        if (initCost) {
            bestCost = aktCost;
        }
//        if (s.sigma < 1.0) {
//            useGrid = true;
//            setUpGrid();
//        } else {
//            s.sigma = 1.0;
//        }
        //s.sigma = Math.min(s.sigma, 3.0);
        if (useGrid) {
            setUpGrid();
        }
//        //setUpGrid();
        if (verbose) {
            System.out.println(iter + " " + s.sigma + " " + aktCost);
        }

        return aktCost;

    }

    private double computeSigmoid(boolean initCost, double maxSigma) {
        //compute sigmoid
        PairOwn[] p = new PairOwn[m];
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < i; j++) {
                p[i * (i - 1) / 2 + j] = new PairOwn(dist(coord[i], coord[j]), g.isNeighbor(i, j));
            }
        }
        s = new SimpleSigmoid(p);
        double aktCost = s.costAllPairs(p);
        if (initCost) {
            bestCost = aktCost;
        }
//        if (s.sigma < 1.0) {
//            useGrid = true;
//            setUpGrid();
//        } else {
//            s.sigma = 1.0;
//        }
        s.sigma = Math.min(s.sigma, maxSigma);

        if (useGrid) {
            setUpGrid();
        }
//        //setUpGrid();
        if (verbose) {
            System.out.println(iter + " " + s.sigma + " " + aktCost);
        }

        return aktCost;

    }

    private double gridThreshold() {
        if (Math.pow(tau, 2) * Math.pow(s.sigma, 4) > 0.05) {
            return 1.33;

        }
        double xiOld = 1.0;
        double xiNew = 0.0;
        boolean converged = false;
        double cc = 1E-3;
        double xi = 1.0;
        int ii = 0;
        while (!converged) {
            ii++;
            double xiTerm = Math.pow(tau, 2) / Math.pow(xiOld, 2);
            double product = 4 * Math.PI * xiTerm * Math.pow(Math.log(2), 2) * Math.pow(s.sigma, 4);
            xiNew = Math.sqrt(-Math.log(product));
            if (Math.abs(xiOld - xiNew) < cc || ii > 2000) {
                converged = true;
            }
            xiOld = xiNew;
        }

        double dist = Math.sqrt(2.0) * s.sigma * xiNew + mu;
        return dist;
    }

    private void setUpGrid() {
        DataUtils du = new DataUtils();
        cutOff = gridThreshold();
        //cutOff = 100;
        minMax = du.minMax(coord);
        numBinsX = (int) Math.ceil((minMax[0][1] - minMax[0][0]) / cutOff);
        numBinsY = (int) Math.ceil((minMax[1][1] - minMax[1][0]) / cutOff);
        //debug
//        numBinsX = 1;
//        numBinsY = 1;
        useGrid = numBinsX > 1 && numBinsY > 1;

        //System.out.println(cutOff + " " + numBinsX + " " + numBinsY);
//        numBinsX = Math.min(numBinsX, 1);
//        numBinsY = Math.min(numBinsY, 1);
        int[][] numPoints = new int[numBinsX][numBinsY];
        for (int i = 0; i < n; i++) {
            //compute grid cell of point
            //first bin: [0, 0.01[, last bin [0.9, 1.0]
            int iIndex = Math.min((int) Math.floor((coord[i][0] - minMax[0][0]) / cutOff), numBinsX - 1);
            int jIndex = Math.min((int) Math.floor((coord[i][1] - minMax[1][0]) / cutOff), numBinsY - 1);
//
            numPoints[iIndex][jIndex]++;
        }

        pId = new int[numBinsX][numBinsY][];
        for (int i = 0; i < numBinsX; i++) {
            for (int j = 0; j < numBinsY; j++) {
                pId[i][j] = new int[numPoints[i][j]];
            }
        }

        for (int i = 0; i < numBinsX; i++) {
            for (int j = 0; j < numBinsY; j++) {
                numPoints[i][j] = 0;
            }
        }

        for (int i = 0; i < n; i++) {
            //compute grid cell of point
            //first bin: [0, 0.01[, last bin [0.9, 1.0]
            int iIndex = Math.min((int) Math.floor((coord[i][0] - minMax[0][0]) / cutOff), numBinsX - 1);
            int jIndex = Math.min((int) Math.floor((coord[i][1] - minMax[1][0]) / cutOff), numBinsY - 1);
//
            pId[iIndex][jIndex][numPoints[iIndex][jIndex]] = i;
            numPoints[iIndex][jIndex]++;
        }
    }

    //4.6.2020. do specific initialization for movie of final verson.
    public double[][] initForMovie() {
        int seed = 304908421;
        initializeCoords(seed);
        DataUtils du = new DataUtils();
        String fName = "coordi_0.mat";
        du.saveResult(coord, fName);
        //du.saveAsMatlab(coord, "coord", fName);
        System.out.println("NEW INITIALIZATION: SEED " + seed);
        double[][] aktCoord = initializationPhaseForMovie();
        return aktCoord;
    }

    //nunTry: number of random inits, iterSigmoid: number of iterations before re-computation of sigmoid and cost; actually this is only to check if finished, sigmoid does not need to be determined
    public double[][] getBestInitialization() {
        System.out.println("--------Initialization Phase------------");
        //Random rr = new Random(10); //standard
        Random rr = new Random(3); //other try
        double[][] bestCoord = new double[n][d];
        double bc = Double.MAX_VALUE;
        for (int i = 0; i < numInits; i++) {
            reset();
            int s = rr.nextInt();
            System.out.println("NEW INITIALIZATION: SEED " + s);
            initializeCoords(s);
            double[][] aktCoord = initializationPhase();
            if (this.bestCost < bc) {
                bestCoord = aktCoord;
                bc = this.bestCost;
            }

        }
        return bestCoord;
    }

    //nunTry: number of random inits, iterSigmoid: number of iterations before re-computation of sigmoid and cost; actually this is only to check if finished, sigmoid does not need to be determined
    public double[][] getBestInitializationGrid() {
        System.out.println("--------Initialization Phase------------");
        Random rr = new Random(10);
        double[][] bestCoord = new double[n][d];
        double bc = Double.MAX_VALUE;
        useGrid = true;
        for (int i = 0; i < numInits; i++) {
            reset();
            int ss = rr.nextInt();
            System.out.println("NEW INITIALIZATION: SEED " + ss);
            //initializeCoords(s); //max sigma 1
            initializeCoordsFreeSigma(ss);
            s.sigma = Math.min(s.sigma, 1.0);

            double[][] aktCoord = stabilizationPhase();
            if (this.bestCost < bc) {
                bestCoord = aktCoord;
                bc = this.bestCost;
            }

        }
        return bestCoord;
    }

    public double[][] finishPhase() {
        System.out.println("--------Polishing Phase------------");
        DataUtils du = new DataUtils();
        boolean improvement = true;
        improvementInCurrentPhase = false;
        double lastCost = Double.MAX_VALUE;
        int count = 0;
        int start = 4590;
        int step = 100;
        int writeCount = 0;
        double[][] bestCoord = new double[n][d];

        //s.sigma = 0.2; //changed for fake display
        bestCost = Double.MAX_VALUE; //changed for fake display

        while (iter < maxIteration && improvement == true) {
            iter++;
            finalPhaseOne();
            updateDW();
            phaseTwo();
            updateCoords();
            if (iter % iterSigmoid == 0) {
//                if(iter == 400)
//                    System.out.println("m");
                double aktCost = computeSigmoid(false); //this is the original version

                //double aktCost = computeSigmoid(false, 0.3);
                //double aktCost = computeCost(); //no update of sigma

                if (saveIntermediateResults) {
                    int name = start + writeCount * step;
                    //writeCount++;
                    //String fName = "coordi_" + name + ".mat";
                    String fName = "coord_" + iter + ".txt"; //changed for fake movie
                   // du.saveAsMatlab(coord, "coord", fName);
                    du.saveResult(coord, fName);
                }
//                if (verbose && display) {
//                    Visualization v = new Visualization(g);
//                    String ss = Integer.toString(iter);
//                    v.displayCoordNew(coord, ss);
//                }

                if ((aktCost + convConst) < bestCost) {
                    improvement = true;
                    count = 0;
                    bestCoord = coord.clone();
                    bestCost = aktCost;
                    lastCost = aktCost;
                    improvementInCurrentPhase = true;
                } else {
                    count++;
                    if (count == localMinEscape) {
                        improvement = false;
                    }
                }

            }
        }
        return bestCoord;
    }

    //can only be used together with grid, i.e. only 2D and if more than one gridcell. Check this before calling
    public double[][] stabilizationPhase() {
        System.out.println("--------Stabilization Phase------------");
        DataUtils du = new DataUtils();
        boolean improvement = true;
        improvementInCurrentPhase = false;
        double lastCost = Double.MAX_VALUE;
        int count = 0;
        double[][] bestCoord = new double[n][d];
        while (iter < maxIteration && improvement == true) {
            iter++;
            stabilizationPhaseOne(iter);
            updateDW();
            phaseTwo();
            updateCoords();
            if (iter % iterSigmoid == 0) {
                double aktCost = computeSigmoid(false);
                s.sigma = Math.min(s.sigma, 1.0);
                if ((aktCost + convConst) < bestCost) {
                    improvement = true;
                    improvementInCurrentPhase = true;
                    count = 0;
                    bestCoord = coord.clone();
                    bestCost = aktCost;
                    lastCost = aktCost;
                } else {
                    count++;
                    if (count == localMinEscape) {
                        improvement = false;
                    }
                }

            }
        }
        return bestCoord;
    }

    private double[][] copyDouble(double[][] d) {
        double[][] res = new double[d.length][d[0].length];
        for (int i = 0; i < d.length; i++) {
            for (int j = 0; j < d[i].length; j++) {
                res[i][j] = d[i][j];
            }
        }
        return res;
    }

    //iterationsSigmoid: number of loops over all edges before recomputation of the sigmoid, int numTry: number of trys used to escape local cost minimum
    public double[][] initializationPhaseForMovie() {
        improvementInCurrentPhase = false;
        DataUtils du = new DataUtils();
        boolean improvement = true;
        double lastCost = Double.MAX_VALUE;
        int count = 0;
        double[][] bestCoord = new double[n][d];
        while (iter <= 4500) {
            iter++;
            initializationPhaseOne(iter);
            updateDW();
            phaseTwo();
            updateCoords();
            if (iter < 10 || iter % 10 == 0) {
                String fName = "coordi_" + iter + ".mat";
                //du.saveAsMatlab(coord, "coord", fName);
                du.saveResult(coord, fName);
            }

            if (iter % iterSigmoid == 0 || iter < 10) {
                double aktCost = computeCost();
                //s.sigma = Math.min(s.sigma, 1.0);
                if ((aktCost + convConst) < bestCost) {
                    improvement = true;
                    count = 0;
                    bestCoord = coord.clone();
                    bestCost = aktCost;
                    lastCost = aktCost;
                    improvementInCurrentPhase = true;
                } else {
                    count++;
                    if (count == localMinEscape) {
                        improvement = false;
                    }
                }

            }
        }
        return bestCoord;
    }

    //iterationsSigmoid: number of loops over all edges before recomputation of the sigmoid, int numTry: number of trys used to escape local cost minimum
    public double[][] initializationPhase() {
        improvementInCurrentPhase = false;
        DataUtils du = new DataUtils();
        boolean improvement = true;
        double lastCost = Double.MAX_VALUE;
        int count = 0;
        double[][] bestCoord = new double[n][d];
        while (iter < maxIteration && improvement == true) {
            iter++;
            initializationPhaseOne(iter);
            updateDW();
            phaseTwo();
            updateCoords();
            if (iter % iterSigmoid == 0) {
                double aktCost = computeCost();
                //s.sigma = Math.min(s.sigma, 1.0);
                if ((aktCost + convConst) < bestCost) {
                    improvement = true;
                    count = 0;
                    bestCoord = coord.clone();
                    bestCost = aktCost;
                    lastCost = aktCost;
                    improvementInCurrentPhase = true;
                } else {
                    count++;
                    if (count == localMinEscape) {
                        improvement = false;
                    }
                }

            }
        }
        return bestCoord;
    }

    public void run() {
        DataUtils du = new DataUtils();
        while (iter < maxIteration) {
            iter++;
            //System.out.println(iter);
//            if (iter == 27) {
//                System.out.println(iter);
//            }
            initializationPhaseOne(iter);
//            if(iter == 102)
//                System.out.println("m");
            updateDW();
            phaseTwo();
            updateCoords();
//            if (iter % 100 == 0) {
//                computeSigmoid();
//            }
            if (iter % 1000 == 0) {
                computeSigmoid(false);
                du.saveResult(coord, "res_" + iter + ".txt");
                //du.saveAsMatlab(coord, "coord", "res_" + iter + ".mat");
//                Visualization v = new Visualization(g);
//                v.displayCoordNew(coord, " ");

            }
//            if (iter % 1000 == 0) {
//
//                //du.saveAsMatlab(coord, "coord", "res_" + iter + ".mat");
//            }

        }
//        Visualization v = new Visualization(g);
//        v.displayCoordNew(coord, " ");
    }

    public void updateCoords() {
        for (CoordUpdateResult cr1 : cr) {
            for (int j = 0; j < cr1.index.length; j++) {
                for (int k = 0; k < d; k++) {
                    if (!Double.isNaN(cr1.coord[j][k])) //System.out.println("update Coords: NaN: " + cr[i].index[j]) ;
                    {
                        coord[cr1.index[j]][k] = cr1.coord[j][k];
                    }
                }
            }
        }
        if (useGrid) {
            setUpGrid();
        }
    }

    public void updateDW() {
        distances = new double[n][n];
        weights = new double[n][n];
        for (int i = 0; i < dr.length; i++) {
            for (int j = 0; j < dr[i].d.length; j++) {
                if (dr[i].ij[j][0] != -1) {
                    distances[dr[i].ij[j][0]][dr[i].ij[j][1]] = distances[dr[i].ij[j][1]][dr[i].ij[j][0]] = dr[i].d[j];
                    weights[dr[i].ij[j][0]][dr[i].ij[j][1]] = weights[dr[i].ij[j][1]][dr[i].ij[j][0]] = dr[i].w[j];
                }
            }

        }
    }

    public void updateSumWeights() {
        sumWeights = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                sumWeights += weights[i][j];
            }
        }
    }

    public void stabilizationPhaseOne(int iter) {
        int numEdgesPerThread = g.getEdgeCount() / numT;
        int remain = g.getEdgeCount() - (numT * numEdgesPerThread);
        int counter = 0;
        Random rr = new Random(iter + 5);
        ExecutorService executor = Executors.newFixedThreadPool(numT);
        Collection collection = new ArrayList();
        for (int i = 0; i < numT - 1; i++) {
            Pair[] edges = new Pair[numEdgesPerThread];
            for (int j = 0; j < edges.length; j++) {
                edges[j] = g.getEndpoints(counter);
                counter++;
            }
            SamplingGridDWThread gt = new SamplingGridDWThread(edges, this, new Random(rr.nextInt()));
            collection.add(gt);
        }

        Pair[] edges = new Pair[numEdgesPerThread + remain];
        for (int j = 0; j < edges.length; j++) {
            edges[j] = g.getEndpoints(counter);
            counter++;
        }
        SamplingGridDWThread gt = new SamplingGridDWThread(edges, this, new Random(rr.nextInt()));
        collection.add(gt);

        try {
            List<Future<DWUpdateResult>> futures = executor.invokeAll(collection);
            boolean allDone = true;
            int cc = 0;
            for (Future<DWUpdateResult> future : futures) {
                DWUpdateResult dd = future.get();
                dr[cc] = dd;
                cc++;

                if (!future.isDone()) {
                    allDone = false;
                }
            }
            if (allDone) {

            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        executor.shutdown();
    }

    public void initializationPhaseOne(int iter) {
        int numEdgesPerThread = g.getEdgeCount() / numT;
        int remain = g.getEdgeCount() - (numT * numEdgesPerThread);
        int counter = 0;
        Random rr = new Random(iter);
        ExecutorService executor = Executors.newFixedThreadPool(numT);
        Collection collection = new ArrayList();

        for (int i = 0; i < numT - 1; i++) {
            Pair[] edges = new Pair[numEdgesPerThread];
            for (int j = 0; j < edges.length; j++) {
                edges[j] = g.getEndpoints(counter);
                counter++;
            }
            SamplingDWThread t = new SamplingDWThread(edges, this, new Random(rr.nextInt()));
            collection.add(t);
        }

        Pair[] edges = new Pair[numEdgesPerThread + remain];
        for (int j = 0; j < edges.length; j++) {
            edges[j] = g.getEndpoints(counter);
            counter++;
        }
        SamplingDWThread t = new SamplingDWThread(edges, this, new Random(rr.nextInt()));
        collection.add(t);

        try {
            List<Future<DWUpdateResult>> futures = executor.invokeAll(collection);
            boolean allDone = true;
            int cc = 0;
            for (Future<DWUpdateResult> future : futures) {
                DWUpdateResult dd = future.get();
                dr[cc] = dd;
                cc++;

                if (!future.isDone()) {
                    allDone = false;
                }
            }
            if (allDone) {

            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        executor.shutdown();
    }

    public void phaseTwo() {
        int numNodesPerThread = n / numT;
        int remain = n - (numT * numNodesPerThread);
        int counter = 0;
        ExecutorService executor = Executors.newFixedThreadPool(numT);
        Collection collection = new ArrayList();
        for (int i = 0; i < numT - 1; i++) {
            int[] nodes = new int[numNodesPerThread];
            for (int j = 0; j < nodes.length; j++) {
                nodes[j] = counter;
                counter++;
            }
            CoordThread t = new CoordThread(nodes, this);
            collection.add(t);
        }

        int[] nodes = new int[numNodesPerThread + remain];
        for (int j = 0; j < nodes.length; j++) {
            nodes[j] = counter;
            counter++;
        }
        CoordThread t = new CoordThread(nodes, this);
        collection.add(t);

        try {
            List<Future<CoordUpdateResult>> futures = executor.invokeAll(collection);
            boolean allDone = true;
            int cc = 0;
            for (Future<CoordUpdateResult> future : futures) {
                CoordUpdateResult cu = future.get();
                cr[cc] = cu;
                cc++;
                if (!future.isDone()) {
                    allDone = false;
                }
            }
            if (allDone) {

            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        executor.shutdown();
    }

    public static double dist(double[] x, double[] y) {
        int d = x.length;
        double result = 0;
        for (int i = 0; i < d; i++) {
            result += (x[i] - y[i]) * (x[i] - y[i]);
        }
        return Math.sqrt(result);
    }

    //optimize objective function in the final stage
    private void finalPhaseOne() {
        int numEdgesPerThread = g.getEdgeCount() / numT;
        int numNodesPerThread = g.getVertexCount() / numT;
        int remainE = g.getEdgeCount() - (numT * numEdgesPerThread);
        int remainN = g.getVertexCount() - (numT * numNodesPerThread);
        int counterE = 0;
        int counterN = 0;
        ExecutorService executor = Executors.newFixedThreadPool(numT);
        Collection collection = new ArrayList();
        for (int i = 0; i < numT - 1; i++) {
            Pair[] edges = new Pair[numEdgesPerThread];
            int[] nodes = new int[numNodesPerThread];
            for (int j = 0; j < edges.length; j++) {
                edges[j] = g.getEndpoints(counterE);
                counterE++;
            }
            for (int j = 0; j < nodes.length; j++) {
                nodes[j] = counterN;
                counterN++;
            }
            if (useGrid) {
                FullGridDWThread gt = new FullGridDWThread(edges, nodes, this);
                collection.add(gt);
            } else {
                FullDWThread gt = new FullDWThread(edges, nodes, this);
                collection.add(gt);
            }
        }

        Pair[] edges = new Pair[numEdgesPerThread + remainE];
        for (int j = 0; j < edges.length; j++) {
            edges[j] = g.getEndpoints(counterE);
            counterE++;
        }
        int[] nodes = new int[numNodesPerThread + remainN];
        for (int j = 0; j < nodes.length; j++) {
            nodes[j] = counterN;
            counterN++;
        }
        if (useGrid) {
            FullGridDWThread gt = new FullGridDWThread(edges, nodes, this);
            collection.add(gt);
        } else {
            FullDWThread gt = new FullDWThread(edges, nodes, this);
            collection.add(gt);
        }
        try {
            List<Future<DWUpdateResult>> futures = executor.invokeAll(collection);
            boolean allDone = true;
            int cc = 0;
            for (Future<DWUpdateResult> future : futures) {
                DWUpdateResult dd = future.get();
                dr[cc] = dd;
                cc++;

                if (!future.isDone()) {
                    allDone = false;
                }
            }
            if (allDone) {

            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        executor.shutdown();
    }

}
