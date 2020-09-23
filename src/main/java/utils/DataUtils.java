/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package utils;

import com.jmatio.io.MatFileReader;
import com.jmatio.io.MatFileWriter;
import com.jmatio.types.MLArray;
import com.jmatio.types.MLDouble;
import edu.uci.ics.jung.algorithms.cluster.WeakComponentClusterer;
import edu.uci.ics.jung.algorithms.filters.FilterUtils;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import edu.uci.ics.jung.graph.util.Pair;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.PrintStream;
import java.io.StreamTokenizer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Random;
import java.util.Scanner;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
//import visu.VisuImage;

/**
 * @author claudia
 */
public class DataUtils {

    HashMap<String, Integer> i3;

    public double[][] readMatlabMatrix(String dir, String variableName) {
        double[][] dataM = new double[1][1];

        MatFileReader mfr = null;
        try {
            mfr = new MatFileReader(dir);
        } catch (IOException e) {
        }

        if (mfr != null) {
            dataM = ((MLDouble) mfr.getMLArray(variableName)).getArray();
        }
        return dataM;
    }

    public void removeDegreeOneNodes(Graph g, double[][] labels) {
        Graph<Integer, Integer> gg = new UndirectedSparseGraph<Integer, Integer>();
        HashMap<Integer, Integer> oldNameToNewName = new HashMap<Integer, Integer>();
        HashMap<Integer, Integer> newNameToOldName = new HashMap<Integer, Integer>();
        int counter = 0;
        for (int i = 0; i < g.getVertexCount(); i++) {
            if (g.degree(i) > 1) {
                gg.addVertex(counter);
                oldNameToNewName.put(i, counter);
                newNameToOldName.put(counter, i);
                counter++;
            }
        }
        System.out.println(counter + " nodes remaining.");
        double[][] labelsNew = new double[counter][1];
        int edgeCounter = 0;
        for (int i = 0; i < gg.getVertexCount(); i++) {
            labelsNew[i][0] = labels[newNameToOldName.get(i)][0];
        }
        //add all Edges
        for (int i = 0; i < gg.getVertexCount(); i++) {
            for (int j = 0; j < gg.getVertexCount(); j++) {

                if (i > j && g.isNeighbor(newNameToOldName.get(i), newNameToOldName.get(j))) {
                    gg.addEdge(edgeCounter, i, j);
                    edgeCounter++;
                }
            }
        }
        IO ea = new IO();
        ea.writeGraphToMatlab(gg, "reducedGraph");
        saveAsMatlab(labelsNew, "labels", "reducedLabels.mat");

    }

    public void writeCSVFile(double[][] d, String filename) {
        FileOutputStream output = null;
        try {
            File target = new File(filename);
            output = new FileOutputStream(target);
            PrintStream writer = new PrintStream(output);
            //ObjectOutputStream myStream = new ObjectOutputStream(output);

            for (int i = 0; i < d.length; i++) {
                for (int j = 0; j < d[i].length - 1; j++) {
                    writer.print(d[i][j] + ",");
                }
                writer.println(d[i][d[i].length - 1]);
            }
            //System.out.println("done.");
//            for (int i = 0; i < g.getEdgeCount(); i++) {
//                Pair bla = g.getEndpoints(i);
//                writer.println((Integer) bla.getSecond());
//            }

        } catch (Exception ex) {
            System.out.println(ex.getMessage());
//             Logger.getLogger(DataUtils.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                output.close();
            } catch (IOException ex) {
                System.out.println(ex.getMessage());
//                   Logger.getLogger(DataUtils.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    public void writeLabel(double[] d, String filename) {
        FileOutputStream output = null;
        try {
            File target = new File(filename);
            output = new FileOutputStream(target);
            PrintStream writer = new PrintStream(output);

            for (int i = 0; i < d.length; i++) {
                //Pair bla = g.getEndpoints(i);
                //System.out.println(i+" "+ bla.getFirst() + " " + bla.getSecond());
                writer.println(d[i]);
            }
            //System.out.println("done.");
//            for (int i = 0; i < g.getEdgeCount(); i++) {
//                Pair bla = g.getEndpoints(i);
//                writer.println((Integer) bla.getSecond());
//            }

        } catch (Exception ex) {
            System.out.println(ex.getMessage());
//             Logger.getLogger(DataUtils.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                output.close();
            } catch (IOException ex) {
                System.out.println(ex.getMessage());
//                   Logger.getLogger(DataUtils.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    public void writeAdjacencyList(Graph g, String filename) {
        FileOutputStream output = null;
        System.out.println(g.getVertexCount() + " " + g.getEdgeCount());
        try {
            File target = new File(filename);
            output = new FileOutputStream(target);
            PrintStream writer = new PrintStream(output);
            //ObjectOutputStream myStream = new ObjectOutputStream(output);

            writer.println(g.getVertexCount());
            writer.println(g.getEdgeCount());

            for (int i = 0; i < g.getEdgeCount(); i++) {
                Pair bla = g.getEndpoints(i);
                //System.out.println(i+" "+ bla.getFirst() + " " + bla.getSecond());
                writer.println((Integer) bla.getFirst() + " " + (Integer) bla.getSecond());
            }
            //System.out.println("done.");
//            for (int i = 0; i < g.getEdgeCount(); i++) {
//                Pair bla = g.getEndpoints(i);
//                writer.println((Integer) bla.getSecond());
//            }

        } catch (Exception ex) {
            System.out.println(ex.getMessage());
//             Logger.getLogger(DataUtils.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                output.close();
            } catch (IOException ex) {
                System.out.println(ex.getMessage());
//                   Logger.getLogger(DataUtils.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    //    public void writeAdjacencyList(Graph g, String filename) {
//        FileOutputStream output = null;
//        System.out.println(g.getVertexCount() + " " + g.getEdgeCount());
//        try {
//            File target = new File(filename);
//            output = new FileOutputStream(target);
//            PrintStream writer = new PrintStream(output);
//            //ObjectOutputStream myStream = new ObjectOutputStream(output);
//
//            writer.println(g.getVertexCount());
//            writer.println(g.getEdgeCount());
//
//            for (int i = 0; i < g.getEdgeCount(); i++) {
//                Pair bla = g.getEndpoints(i);
//                writer.println((Integer) bla.getFirst());
//            }
//            for (int i = 0; i < g.getEdgeCount(); i++) {
//                Pair bla = g.getEndpoints(i);
//                writer.println((Integer) bla.getSecond());
//            }
//
//        } catch (Exception ex) {
//            Logger.getLogger(DataUtils.class.getName()).log(Level.SEVERE, null, ex);
//        } finally {
//            try {
//                output.close();
//            } catch (IOException ex) {
//                Logger.getLogger(DataUtils.class.getName()).log(Level.SEVERE, null, ex);
//            }
//        }
//    }
    //set up hashmap to index all potential triangles
    public void setUp3IndexTable(int numObj) {
        i3 = new HashMap<String, Integer>();
        int counter = 0;
        for (int index1 = 0; index1 <= numObj - 3; index1++) {
            for (int index2 = index1 + 1; index2 <= numObj - 2; index2++) {
                for (int index3 = index2 + 1; index3 <= numObj - 1; index3++) {
                    String key = ((Integer) index1).toString() + ((Integer) index2).toString() + ((Integer) index3).toString();
                    //System.out.println(key);
                    i3.put(key, counter);
                    //test
                    //System.out.println(i3.get(key));
                    //test
                    counter++;
                }
            }
        }
    }

    public int getIndex3(int index1, int index2, int index3, int numObj) {
        Integer[] indices = {index1, index2, index3};
        Arrays.sort(indices);
        int first = indices[0];
        int second = indices[1];
        int third = indices[2];
        String key = ((Integer) first).toString() + ((Integer) second).toString() + ((Integer) third).toString();
        int bla = -1;
        try {
            bla = i3.get(key);
        } catch (java.lang.NullPointerException e) {
            System.out.println("m");

        }
        return bla;
    }

    public void saveAsMatlab2(double[][] d, double[][] e, String vName_d, String vName_e, String fName) {
        MLDouble v = new MLDouble(vName_d, d);
        MLDouble v1 = new MLDouble(vName_e, e);
        ArrayList ll = new ArrayList();
        ll.add(v);
        ll.add(v1);
        MatFileWriter mw = new MatFileWriter();
        try {
            mw.write(fName, ll);
        } catch (IOException ex) {
        }
        //System.out.println("Matlab matrix saved.");
    }

    public HashMap nodenameToIndexAmazon() {
        HashMap nodeNameToIndex = new HashMap<Integer, Integer>();
        int numNodes = 334863;
        int numEdges = 925872;
        String path = "com-amazon.ungraph.txt";
        int counter = 0;
        try {
            File source = new File(path);
            StreamTokenizer st1 = new StreamTokenizer(new BufferedReader(new FileReader(source)));
            st1.eolIsSignificant(false);
            st1.parseNumbers();
            while (st1.ttype != StreamTokenizer.TT_EOF) {
                try {
                    st1.nextToken();

                } catch (IOException ex) {
                    Logger.getLogger(DataUtils.class
                            .getName()).log(Level.SEVERE, null, ex);
                }
                int nodeName = (int) st1.nval;
                if (!nodeNameToIndex.containsKey(nodeName)) {
                    nodeNameToIndex.put(nodeName, counter);
                    counter++;
                }

            }
            System.out.println(counter + " " + numNodes);

        } catch (FileNotFoundException ex) {
            Logger.getLogger(DataUtils.class
                    .getName()).log(Level.SEVERE, null, ex);
        }
        return nodeNameToIndex;
    }

    public Graph readAmazon() {
        int numNodes = 334863;
        int numEdges = 925872;
        String path = "com-amazon.ungraph.txt";
        String communities = "com-amazon.top5000.cmty.txt";
        Graph g = new UndirectedSparseGraph<Integer, Integer>();
        for (int i = 0; i < numNodes; i++) {
            g.addVertex(i);
        }
        HashMap<Integer, Integer> nodenameToIndex = nodenameToIndexAmazon();
        System.out.println(g.getVertexCount());
        int edgeCounter = 0;
        try {
            File source = new File(path);
            StreamTokenizer st1 = new StreamTokenizer(new BufferedReader(new FileReader(source)));
            st1.eolIsSignificant(false);
            st1.parseNumbers();
            while (st1.ttype != StreamTokenizer.TT_EOF) {
                try {
                    st1.nextToken();

                } catch (IOException ex) {
                    Logger.getLogger(DataUtils.class
                            .getName()).log(Level.SEVERE, null, ex);
                }
                int vertex1 = (int) st1.nval;
                try {
                    st1.nextToken();

                } catch (IOException ex) {
                    Logger.getLogger(DataUtils.class
                            .getName()).log(Level.SEVERE, null, ex);
                }
                int vertex2 = (int) st1.nval;
                if (edgeCounter < numEdges) {
                    if (g.getVertexCount() > numNodes) {
                        System.out.println(g.getVertexCount() + " " + edgeCounter);
                    }
                    if (vertex1 != vertex2) {
                        g.addEdge(edgeCounter, nodenameToIndex.get(vertex1), nodenameToIndex.get(vertex2));
                        edgeCounter++;
                    }
                }

            }
            System.out.println(g.getVertexCount() + " " + g.getEdgeCount());

        } catch (FileNotFoundException ex) {
            Logger.getLogger(DataUtils.class
                    .getName()).log(Level.SEVERE, null, ex);
        }
        double[] communityLabel = new double[numNodes];
        String filename = "amazonCommunities.txt";
        String filename1 = "amazonLabels.txt";
        int row = 0;
        try {
            File target = new File(filename);
            FileOutputStream output = new FileOutputStream(target);
            output = new FileOutputStream(target);
            PrintStream writer = new PrintStream(output);
            File source = new File(communities);
            StreamTokenizer st1 = new StreamTokenizer(new BufferedReader(new FileReader(source)));
            st1.eolIsSignificant(true);
            st1.parseNumbers();
            int counter = 0;
            boolean newRow = false;
            while (st1.ttype != StreamTokenizer.TT_EOF) {
                while (st1.ttype != StreamTokenizer.TT_EOL || newRow == true) {
                    newRow = false;
                    try {
                        st1.nextToken();
                        if (st1.ttype == StreamTokenizer.TT_NUMBER) {

                            int nodeIndex = nodenameToIndex.get((int) st1.nval);
                            writer.print(nodeIndex + " ");
                            communityLabel[nodeIndex] = row;
                            counter++;
                        }
                    } catch (IOException ex) {
                        Logger.getLogger(DataUtils.class
                                .getName()).log(Level.SEVERE, null, ex);
                    }
                }
                row++;
                writer.println();
                newRow = true;
            }
            System.out.println(g.getVertexCount() + " " + g.getEdgeCount());

        } catch (FileNotFoundException ex) {
            Logger.getLogger(DataUtils.class
                    .getName()).log(Level.SEVERE, null, ex);
        }
        System.out.println("communities written");
        //saveAsMatlab(communityLabel, "label", "communityLabel.mat");
        writeLabel(communityLabel, "amazonLabels.txt");
        System.out.println("labels written.");

        return g;
    }

    public void writeLineInputFile(Graph g, String name) {
        try {
            File target = new File(name + "lineInput.txt");
            FileOutputStream output = new FileOutputStream(target);
            PrintStream writer = new PrintStream(output);
            for (int i = 0; i < g.getEdgeCount(); i++) {
                Pair<Integer> akt = g.getEndpoints(i);
                String name1 = akt.getFirst().toString();
                String name2 = akt.getSecond().toString();
                writer.println(name1 + " " + name2 + " " + 1);
                writer.println(name2 + " " + name1 + " " + 1);
            }

        } catch (FileNotFoundException e) {
            System.err.println("Error: " + e);
        }
    }

    public void writeArff(double[][] d, int[] labels, int numCl, String dir) {
        int dim = d[0].length;
        try {
            File target = new File(dir + ".arff");
            FileOutputStream output = new FileOutputStream(target);
            PrintStream writer = new PrintStream(output);
            writer.println("@relation test");
            writer.println();
            for (int i = 0; i < dim; i++) {
                writer.println("@attribute " + "a" + i + " real");
            }
            writer.print("@attribute class {");
            for (int i = 0; i < numCl; i++) {
                writer.print((i + 1) + ",");
            }
            //writer.print(numCl);
            writer.println("}");
            writer.println();
            writer.println("@data");
            //int numV = 5;
            //for(int k = 0; k < numV; k++){
            for (int i = 0; i < d.length; i++) {
                for (int j = 0; j < dim; j++) {
                    writer.print(d[i][j] + ",");
                }
                //writer.println(d[i][dim - 1]);
                writer.println(labels[i]);
            }
            //}
        } catch (FileNotFoundException e) {            //Es gibt kein solches File
            System.err.println("Error: " + e);
        }
    }

    public void procrustesEvaluation(double[][] gt, double[][] res) {
        for (int i = 0; i < gt.length; i++) {
            System.out.println(gt[i][0] + " " + gt[i][1]);
            System.out.println(res[i][0] + " " + res[i][1]);
            System.out.println();
        }
    }

    public void order1Evaluation(Graph g, double[][] coord) {
        int numDist = (g.getVertexCount() * (g.getVertexCount() - 1)) / 2;
        double[][] distEdges = new double[g.getEdgeCount()][1];
        double[][] distNonEdges = new double[numDist - g.getEdgeCount()][1];
        //double[][] isEdge = new double[numDist][1];
        int counterE = 0;
        int counterNE = 0;
        for (int i = 0; i < g.getVertexCount(); i++) {
            for (int j = i + 1; j < g.getVertexCount(); j++) {
                double aktDist = euclideanDistance(coord[i], coord[j]);
                if (g.isNeighbor(i, j)) {
                    distEdges[counterE][0] = aktDist;
                    counterE++;

                } else {
                    distNonEdges[counterNE][0] = aktDist;
                    counterNE++;

                }

            }
        }

//           // dist = scale(dist);
//            double[][] dd = new double[dist.length][1];
//            for(int i = 0; i < dd.length; i++)
//                dd[i][0] = dist[i];
//
        saveAsMatlab2(distEdges, distNonEdges, "distEdges", "distNonEdges", "order1.mat");

    }

    public double euclideanDistance(double[] o, double[] p) {
        double sum = 0.0;
        int dim = o.length;
        for (int i = 0; i < dim; i++) {
            double dist = ((o[i] - p[i]) * (o[i] - p[i]));
            sum = sum + dist;
        }
        return Math.sqrt(sum);
    }

    public Graph getLargestComponent(Graph g) {
        FilterUtils filt = new FilterUtils();
        WeakComponentClusterer<Number, Number> clusterer = new WeakComponentClusterer<Number, Number>();
        Set<Set<Number>> clusterset = clusterer.transform(g);
        Set<Number> largest = Collections.EMPTY_SET;
        for (Set<Number> cluster : clusterset) {
            if (cluster.size() > largest.size()) {
                largest = cluster;
            }
        }
        Graph<Integer, Integer> res = new UndirectedSparseGraph<Integer, Integer>();
        boolean[] isIncluded = new boolean[g.getVertexCount()];
        Object[] nodeNames = largest.toArray();
        int[] index = new int[nodeNames.length];
        //add all nodes
        for (int i = 0; i < index.length; i++) {
            Integer name = (Integer) nodeNames[i];
            isIncluded[name] = true;
            res.addVertex(i);
            index[i] = i;
        }
        for (int i = 0; i < isIncluded.length; i++) {
            if (!isIncluded[i]) {
                System.out.println("missing node: " + i);
            }
        }
        //add all Edges
        int edgeCounter = 0; //changed on 27.05.2020 for display of corona graph

        for (int i = 0; i < nodeNames.length; i++) {
            for (int j = 0; j < nodeNames.length; j++) {
                if (i > j && g.isNeighbor(nodeNames[i], nodeNames[j])) {
                    //res.addEdge(getIndex(i, j, nodeNames.length), i, j);
                    res.addEdge(edgeCounter, i, j);
                    edgeCounter++;
                }
            }
        }
//        //check if graph is connected
//        for(int i = 0; i < res.getVertexCount(); i++){
//            int bla = res.getNeighborCount(i);
//           // if(bla == 0)
//                System.out.println(i + " " + bla);
//        }

        return res;
    }

    //31.07.2017: write only the ground truth coordinates corresponding to the connected nodes
    public void getLargestComponent(Graph g, double[][] gt) {
        FilterUtils filt = new FilterUtils();
        WeakComponentClusterer<Number, Number> clusterer = new WeakComponentClusterer<Number, Number>();
        Set<Set<Number>> clusterset = clusterer.transform(g);
        Set<Number> largest = Collections.EMPTY_SET;
        for (Set<Number> cluster : clusterset) {
            if (cluster.size() > largest.size()) {
                largest = cluster;
            }
        }
        Graph<Integer, Integer> res = new UndirectedSparseGraph<Integer, Integer>();
        boolean[] isIncluded = new boolean[g.getVertexCount()];
        Object[] nodeNames = largest.toArray();
        int[] index = new int[nodeNames.length];
        //add all nodes
        for (int i = 0; i < index.length; i++) {
            Integer name = (Integer) nodeNames[i];
            isIncluded[name] = true;
            res.addVertex(i);
            index[i] = i;
        }

        int counter = 0;

        for (int i = 0; i < isIncluded.length; i++) {
            if (!isIncluded[i]) {
                System.out.println("missing node: " + i);
            } else {
//                gtCoord[counter][0] = gt[counter][0];
//                gtCoord[counter][1] = gt[counter][1];
                counter++;
            }
        }
        double[][] gtCoord = new double[counter][2];
        int cc = 0;
        int cI = 0;
        for (int i = 0; i < isIncluded.length; i++) {
            if (isIncluded[i]) {
                gtCoord[cc][0] = gt[cI][0];
                gtCoord[cc][1] = gt[cI][1];
                //  gtCoord[cc][2] = gt[cI][2];
                cc++;
                cI++;
            } else {
                cI++;
            }

        }

        saveAsMatlab(gtCoord, "gt", "gt.mat");

        //add all Edges
        for (int i = 0; i < nodeNames.length; i++) {
            for (int j = 0; j < nodeNames.length; j++) {
                if (i > j && g.isNeighbor(nodeNames[i], nodeNames[j])) {
                    res.addEdge(getIndex(i, j, nodeNames.length), i, j);
                }
            }
        }
//        //check if graph is connected
//        for(int i = 0; i < res.getVertexCount(); i++){
//            int bla = res.getNeighborCount(i);
//           // if(bla == 0)
//                System.out.println(i + " " + bla);
//        }

//        return res;
    }

    public double[][] createRandomMatrix(int numRows, int numCols) {
        double[][] res = new double[numRows][numCols];
        Random r = new Random(1);
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                res[i][j] = r.nextDouble() * 100;
                if (res[i][j] == 0) {
                    res[i][j] = 0.0001;
                }
            }
        }
        return res;

    }

    public void saveLMLabels(double[][] labels, double[][] lmIndex) {
        double[][] lmLabels = new double[lmIndex.length][1];
        for (int i = 0; i < lmIndex.length; i++) {
            lmLabels[i][0] = labels[(int) lmIndex[i][0]][0];
        }
        saveAsMatlab(lmLabels, "lmLabels", "lmLabels.mat");
    }

    //17.6.14: evaluate accuracy of K-NN query
    public void accuracy(double[] basicSpace, double[][] label, int k, int numObj) {
        double meanAccuracy = 0.0;
        for (int i = 0; i < numObj; i++) {
            int[] neighbors = getKNN(basicSpace, i, k, numObj);
            int countSame = 0;
            for (int j = 0; j < neighbors.length; j++) {
                if (label[neighbors[j]][0] == label[i][0]) {
                    countSame++;
                }
            }
            meanAccuracy += (double) countSame / (double) k;
        }
        meanAccuracy /= (double) numObj;
        System.out.println("K: " + k + " meanAcc: " + meanAccuracy);
    }

    private int[] getKNN(double[] basicSpace, int obj, int k, int numObj) {
        double[] erg = new double[k];      //Array der Grösse MinPts anlegen
        int[] ind = new int[k];

        for (int i = 0; i < k; i++) {
            erg[i] = Double.MAX_VALUE;  //alle auf Max_Value setzen.
            ind[i] = -1;
        }
        for (int i = 0; i < numObj; i++) {   //Durchlauf alle Datenobjekte
            if (i != obj) {
                int c = k - 1;
                double d = basicSpace[getIndex(obj, i, numObj)]; //Distanz zum i-ten Datenobjekt

                while (c >= 0 && d < erg[c]) { //MinPts >= 0 && Distanz zum i-ten Datenobjekt < erg[MinPts-1]
                    if (c < k - 1) {
                        erg[c + 1] = erg[c];
                        ind[c + 1] = ind[c];
                    }
                    erg[c] = d;
                    ind[c] = i;
                    c--;                  //MinPts um 1 reduzieren
                }
            }
        }
        return ind;
    }

    //22.5.14: create small ALOI DS
    public double[][] readCSVFileAloiFour(String path, int numObj, int dim) {
        String greenTin = "615/";
        String orangeTin = "616/";
        String greenBox = "235/";
        String redBox = "249/";
        double[][] res = new double[numObj][dim];
        int counter = 0;
        int counterO = 0;
        try {
            Scanner scanner = new Scanner(new File(path));

            //Set the delimiter used in file
            scanner.useDelimiter(";");
            //scanner.useDelimiter("\n");

            //Get all tokens and store them in some data structure
            //I am just printing them
            while (scanner.hasNext()) {
                String s = scanner.next();
                //System.out.println(s);
                if (s.startsWith(greenTin) || s.startsWith(orangeTin) || s.startsWith(greenBox) || s.startsWith(redBox)) {
                    //String t = scanner.next();
                    String part = "_l";
                    if (s.contains(part)) {
                        counterO++;
                        System.out.println(s);
                        //String tt = scanner.next(); //do not use first column in RGB
                        // System.out.println(tt);
                        for (int i = 0; i < dim; i++) {
                            String ss = scanner.next();
//                            if(i == 7)
//                                System.out.println("m"); //last column of csv file is not correctly interpeted
//                            System.out.println(i + " " + ss);
//                            if (counter == 96) {
//                                System.out.println("m");
//                            }
                            res[counter][i] = Double.valueOf(ss);

                        }
                        counter++;
                    }

                }

            }
            System.out.println("coutnerO" + counterO);
            //Do not forget to close the scanner
            scanner.close();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(DataUtils.class.getName()).log(Level.SEVERE, null, ex);
        }
        return res;
    }

    public void outputAloiLabels(String path, String ending, int numObj, int dim, double[][] l) {
        double[][] res = new double[numObj][dim];
        int counter = 0;

        int[] index = new int[l.length];
        for (int i = 0; i < index.length; i++) {
            index[i] = (int) l[i][0];
        }
        Arrays.sort(index);
        int countL = 0;
        String[] start = new String[index.length];
        for (int i = 0; i < start.length; i++) {
            start[i] = new Integer(index[i]).toString() + "/";
        }

        System.out.println("m");

        try {
            Scanner scanner = new Scanner(new File(path));

            //Set the delimiter used in file
            scanner.useDelimiter(";");
            //scanner.useDelimiter("\n");

            //Get all tokens and store them in some data structure
            //I am just printing them
            while (scanner.hasNext()) {
                String s = scanner.next();
//                if (countL == 100) {
//                    System.out.println("m");
//                }
                //if (s.endsWith(ending) && s.startsWith(start[countL])) {
                if (s.endsWith(ending)) {
                    System.out.println(s);
                    countL++;
                    for (int i = 0; i < dim; i++) {
                        String ss = scanner.next();
//if(i == 7)
//    System.out.println("m");
//System.out.println(i + " " + ss);
                        res[counter][i] = Double.valueOf(ss);

                    }
                    counter++;
                }
            }

            //Do not forget to close the scanner
            scanner.close();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(DataUtils.class.getName()).log(Level.SEVERE, null, ex);
        }
        //return res;
    }

    public double[][] readCoordFile(String path, int numObj, int dim) {
        //int numObj = 114599; //luxemburg
        //int numObj = 9528; //dblp

        double[][] res = new double[numObj][dim];
        int counter = 0;
        int counterd = 0;
        System.out.println("reading file " + path);
        try {
            File source = new File(path);
            StreamTokenizer st1 = new StreamTokenizer(new BufferedReader(new FileReader(path)));
            st1.parseNumbers();
            st1.eolIsSignificant(false);
            while (st1.ttype != StreamTokenizer.TT_EOF) {  //Daten aus File in Vector einlesen
                st1.nextToken();
                //System.out.println(st1.ttype);
                if (st1.ttype == StreamTokenizer.TT_NUMBER) //  System.out.println(st1.nval);
                {
                    res[counter][counterd] = st1.nval;
                    counterd++;
                    // st1.nextToken();
                    if (counterd == dim) {
                        //System.out.println(st1.nval);
                        //res[counter][counterd] = st1.nval;
                        counter++;
                        //System.out.println(counter);
                        counterd = 0;
                    }
                }

            }

        } catch (IOException e) {            //Es gibt kein solches File
            System.err.println("Error: " + e);
        }
        //luxemburg
//        double[][] bla = new double[numObj][dim];
//        for (int i = 0; i < (bla.length - 1); i++) {
//            bla[i][0] = res[i + 1][0];
//            bla[i][1] = res[i + 1][1];
//
//        }
//        bla[bla.length-1][0] = res[0][0];
//        bla[bla.length-1][1] = res[0][1];
        return res;

    }

    public double[][] removeLastCol(double[][] d) {
        double[][] res = new double[d.length][d[0].length - 1];
        for (int i = 0; i < d.length; i++) {
            for (int j = 0; j < res[0].length; j++) {
                res[i][j] = d[i][j];
            }
        }
        return res;
    }


    public double[][] txtToDouble(String path, int dim, int numObj) {
        //int numObj = 114599; //luxemburg
        //int numObj = 9528; //dblp
        //int numObj = 17278; //dblpLarge
        //int numObj = 334863; //amazon
        //int numObj = 115; //football
        //int numObj = 332; //airflights
        //int numObj = 547; //eppstein
        // int numObj = 1024; //tapir
//         int numObj = 2640; //minnesota
//int numObj = 19717;
//int d = 2;

//        int dim = 3;
//
//        String path = "GEM.txt";


        double[][] res = new double[numObj][dim];
        int counter = 0;
        int counterd = 0;
        System.out.println("reading file " + path);
        try {
            File source = new File(path);
            StreamTokenizer st1 = new StreamTokenizer(new BufferedReader(new FileReader(path)));
            st1.parseNumbers();
            st1.eolIsSignificant(false);
            while (st1.ttype != StreamTokenizer.TT_EOF) {  //Daten aus File in Vector einlesen
                st1.nextToken();
                //System.out.println(st1.ttype);
                if (st1.ttype == StreamTokenizer.TT_NUMBER) //  System.out.println(st1.nval);
                {
                    res[counter][counterd] = st1.nval;
                    counterd++;
                    // st1.nextToken();

                    if (counterd == dim) {
                        //System.out.println(st1.nval);
                        //res[counter][counterd] = st1.nval;
                        counter++;
                        if (counter == 17228)
                            System.out.println(counter);
                        counterd = 0;
                    }
                }

            }

        } catch (IOException e) {            //Es gibt kein solches File
            System.err.println("Error: " + e);
        }
        //luxemburg
//        double[][] bla = new double[numObj][dim];
//        for (int i = 0; i < (bla.length - 1); i++) {
//            bla[i][0] = res[i + 1][0];
//            bla[i][1] = res[i + 1][1];
//
//        }
//        bla[bla.length-1][0] = res[0][0];
//        bla[bla.length-1][1] = res[0][1];
        return res;

    }

    public double[][] txtToDouble() {
        //int numObj = 114599; //luxemburg
        //int numObj = 9528; //dblp
        //int numObj = 17278; //dblpLarge
        //int numObj = 334863; //amazon
        //int numObj = 115; //football
        //int numObj = 332; //airflights
        //int numObj = 547; //eppstein
        // int numObj = 1024; //tapir
        // int numObj = 2640; //minnesota
        int numObj = 19719;
        int dim = 128;

        String path = "pubmedLine.txt";
        double[][] res = new double[numObj][dim];
        int counter = 0;
        int counterd = 0;
        System.out.println("reading file " + path);
        try {
            File source = new File(path);
            StreamTokenizer st1 = new StreamTokenizer(new BufferedReader(new FileReader(path)));
            st1.parseNumbers();
            st1.eolIsSignificant(false);
            while (st1.ttype != StreamTokenizer.TT_EOF) {  //Daten aus File in Vector einlesen
                st1.nextToken();
                //System.out.println(st1.ttype);
                if (st1.ttype == StreamTokenizer.TT_NUMBER) //  System.out.println(st1.nval);
                {
                    res[counter][counterd] = st1.nval;
                    counterd++;
                    // st1.nextToken();

                    if (counterd == dim) {
                        //System.out.println(st1.nval);
                        //res[counter][counterd] = st1.nval;
                        counter++;
                        //System.out.println(counter);
                        counterd = 0;
                    }
                }

            }

        } catch (IOException e) {            //Es gibt kein solches File
            System.err.println("Error: " + e);
        }
        //luxemburg
//        double[][] bla = new double[numObj][dim];
//        for (int i = 0; i < (bla.length - 1); i++) {
//            bla[i][0] = res[i + 1][0];
//            bla[i][1] = res[i + 1][1];
//
//        }
//        bla[bla.length-1][0] = res[0][0];
//        bla[bla.length-1][1] = res[0][1];
        return res;

    }

    //    public Graph readDiabetes(){
//        int numNodes = 19717;
//        int numEdges = 44338;
//        double[] classLabel = new double[numNodes];
//        String nodeFile = "Pubmed-Diabetes.NODE.paper";
//    }
    public void graphToTulipAdj(Graph g, String filename) {
        FileOutputStream output = null;
        System.out.println(g.getVertexCount() + " " + g.getEdgeCount());
        try {
            File target = new File(filename);
            output = new FileOutputStream(target);
            PrintStream writer = new PrintStream(output);

            for (int i = 0; i < g.getVertexCount(); i++) {
                for (int j = 0; j < i; j++) {
                    if (g.isNeighbor(i, j)) {
                        writer.print("@ ");
                    } else {
                        writer.print("# ");
                    }
                }
                writer.println(i);
            }

        } catch (Exception ex) {
            Logger.getLogger(DataUtils.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                output.close();
            } catch (IOException ex) {
                Logger.getLogger(DataUtils.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    public Graph readAdjList(String path) {
        Graph g = new UndirectedSparseGraph<Integer, Integer>();
        int numNodes = 0;
        int numEdges = 0;
        try {
            File source = new File(path);
            StreamTokenizer st1 = new StreamTokenizer(new BufferedReader(new FileReader(source)));
            st1.eolIsSignificant(false);
            st1.parseNumbers();
            try {
                st1.nextToken();
                numNodes = (int) st1.nval;
                st1.nextToken();
                numEdges = (int) st1.nval;
            } catch (IOException ex) {
                Logger.getLogger(DataUtils.class.getName()).log(Level.SEVERE, null, ex);
            }
            for (int i = 0; i < numNodes; i++) {
                g.addVertex(i);
            }
            int edgeCounter = 0;
            while (st1.ttype != StreamTokenizer.TT_EOF) {
                try {
                    st1.nextToken();
                } catch (IOException ex) {
                    Logger.getLogger(DataUtils.class.getName()).log(Level.SEVERE, null, ex);
                }
                int vertex1 = (int) st1.nval;
                try {
                    st1.nextToken();
                } catch (IOException ex) {
                    Logger.getLogger(DataUtils.class.getName()).log(Level.SEVERE, null, ex);
                }
                int vertex2 = (int) st1.nval;

                g.addEdge(edgeCounter, vertex1, vertex2);
                edgeCounter++;
            }

            System.out.println(g.getVertexCount() + " " + g.getEdgeCount());

        } catch (FileNotFoundException ex) {
            Logger.getLogger(DataUtils.class.getName()).log(Level.SEVERE, null, ex);
        }
        return g;
    }

    public Graph readLuxemburg() {
        int numNodes = 114599;
        int numEdges = 119666;
        Graph g = new UndirectedSparseGraph<Integer, Integer>();
        for (int i = 0; i < numNodes; i++) {
            g.addVertex(i);
        }
        System.out.println(g.getVertexCount());
        int edgeCounter = 0;
        try {
            String path = "luxembourg_osm.mtx";
            File source = new File(path);
            StreamTokenizer st1 = new StreamTokenizer(new BufferedReader(new FileReader(source)));
            st1.eolIsSignificant(false);
            st1.parseNumbers();
            while (st1.ttype != StreamTokenizer.TT_EOF) {
                try {
                    st1.nextToken();

                } catch (IOException ex) {
                    Logger.getLogger(DataUtils.class
                            .getName()).log(Level.SEVERE, null, ex);
                }
                int vertex1 = (int) st1.nval;
                try {
                    st1.nextToken();

                } catch (IOException ex) {
                    Logger.getLogger(DataUtils.class
                            .getName()).log(Level.SEVERE, null, ex);
                }
                int vertex2 = (int) st1.nval;
                if (edgeCounter < numEdges) {
                    if (g.getVertexCount() > numNodes) {
                        System.out.println(g.getVertexCount() + " " + edgeCounter);
                    }
//                    if(edgeCounter == 72919)
//                        System.out.println("m");
//                    if (!g.containsVertex(vertex1 - 1) || !g.containsVertex(vertex2 - 1)) {
//                        System.out.println("m");
//                    }
                    g.addEdge(edgeCounter, (vertex1 - 1), (vertex2 - 1));
                    edgeCounter++;
                }

            }
            System.out.println(g.getVertexCount() + " " + g.getEdgeCount());

        } catch (FileNotFoundException ex) {
            Logger.getLogger(DataUtils.class
                    .getName()).log(Level.SEVERE, null, ex);
        }
        return g;
    }

    public Graph readLuxemburgChris() {
        int numNodes = 114599;
        int numEdges = 119666;
        Graph g = new UndirectedSparseGraph<Integer, Integer>();
        for (int i = 0; i < numNodes; i++) {
            g.addVertex(i + 1);
        }
        System.out.println(g.getVertexCount());
        int edgeCounter = 0;
        try {
            String path = "luxembourg_osm.mtx";
            File source = new File(path);
            StreamTokenizer st1 = new StreamTokenizer(new BufferedReader(new FileReader(source)));
            st1.eolIsSignificant(false);
            st1.parseNumbers();
            while (st1.ttype != StreamTokenizer.TT_EOF) {
                try {
                    st1.nextToken();

                } catch (IOException ex) {
                    Logger.getLogger(DataUtils.class
                            .getName()).log(Level.SEVERE, null, ex);
                }
                int vertex1 = (int) st1.nval;
                try {
                    st1.nextToken();

                } catch (IOException ex) {
                    Logger.getLogger(DataUtils.class
                            .getName()).log(Level.SEVERE, null, ex);
                }
                int vertex2 = (int) st1.nval;
                if (edgeCounter < numEdges) {
                    if (g.getVertexCount() > numNodes) {
                        System.out.println(g.getVertexCount() + " " + edgeCounter);
                    }
//                    if(edgeCounter == 72919)
//                        System.out.println("m");
//                    if (!g.containsVertex(vertex1 - 1) || !g.containsVertex(vertex2 - 1)) {
//                        System.out.println("m");
//                    }
                    g.addEdge(edgeCounter, vertex1, vertex2);
                    edgeCounter++;
                }

            }
            System.out.println(g.getVertexCount() + " " + g.getEdgeCount());

        } catch (FileNotFoundException ex) {
            Logger.getLogger(DataUtils.class
                    .getName()).log(Level.SEVERE, null, ex);
        }
        return g;
    }

    public double[][] readCSVFile(String path, String ending, int numObj, int dim) {
        double[][] res = new double[numObj][dim];
        int counter = 0;
        try {
            Scanner scanner = new Scanner(new File(path));

            //Set the delimiter used in file
            scanner.useDelimiter(";");
            //scanner.useDelimiter("\n");

            //Get all tokens and store them in some data structure
            //I am just printing them
            while (scanner.hasNext()) {
                String s = scanner.next();

                if (s.endsWith(ending)) {
                    System.out.println(s);
                    for (int i = 0; i < dim; i++) {
                        String ss = scanner.next();
//if(i == 7)
//    System.out.println("m");
//System.out.println(i + " " + ss);
                        res[counter][i] = Double.valueOf(ss);

                    }
                    counter++;
                }
            }

            //Do not forget to close the scanner
            scanner.close();

        } catch (FileNotFoundException ex) {
            Logger.getLogger(DataUtils.class
                    .getName()).log(Level.SEVERE, null, ex);
        }
        return res;
    }

    //27.05.14: read labels for displaying images.
    public String[] readLabels(String path, int numObj) {
        String[] l = new String[numObj];
        try {
            int counter = 0;
            Scanner scanner = new Scanner(new File(path));
            scanner.useDelimiter("\n");
            //I am just printing them
            while (scanner.hasNext()) {
                String s = scanner.next();
                l[counter] = s;
                counter++;

            }
            //        try {
            //            File source = new File(path);
            //            StreamTokenizer st1 = new StreamTokenizer(new BufferedReader(new FileReader(source)));
            //            st1.eolIsSignificant(false);
            //            while (st1.ttype != StreamTokenizer.TT_EOF) {
            //                try {
            //                    //Daten aus File in Vector einlesen
            //                    st1.nextToken();
            //                    System.out.println(st1.ttype);
            //                    System.out.println(st1.nval);
            //                } catch (IOException ex) {
            //                    Logger.getLogger(VisuImage.class.getName()).log(Level.SEVERE, null, ex);
            //                }
            //                if (st1.ttype == StreamTokenizer.TT_WORD) {
            //                    String nn = (String) (st1.sval);
            //                    l[counter] = nn;
            //                    counter++;
            //                }
            //            }
            //
            //        } catch (FileNotFoundException ex) {
            //            Logger.getLogger(VisuImage.class.getName()).log(Level.SEVERE, null, ex);
            //        }

        } catch (FileNotFoundException ex) {
            Logger.getLogger(DataUtils.class
                    .getName()).log(Level.SEVERE, null, ex);
        }
        return l;
    }

    //27.05.14: scale for plotting images
    public int[][] scaleforVisu(double[][] coord, int range) {
        int n = coord.length;
        int d = coord[0].length;
        int[][] res = new int[n][d];

        //determine min and max from data
        double[] min = new double[d];
        double[] max = new double[d];
        for (int i = 0; i < d; i++) {
            min[i] = Double.MAX_VALUE;
            max[i] = -Double.MAX_VALUE;
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < d; j++) {
                if (coord[i][j] < min[j]) {
                    min[j] = coord[i][j];
                }
                if (coord[i][j] > max[j]) {
                    max[j] = coord[i][j];
                }

            }
        }
        int minScaled = Integer.MAX_VALUE;
        int maxScaled = -Integer.MIN_VALUE;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < d; j++) {
                double scaled = (coord[i][j] - min[j]) / (max[j] - min[j]);
                if (scaled > 1 || scaled < 0) {
                    System.out.println("m");
                }

                scaled *= range;
                int iscaled = (int) scaled;
//                double result = scaled - (double) iscaled;
//                if (result > 0.5) {
//                    iscaled++;
//
//                }
                // if (j == 0) {
                res[i][j] = iscaled;
                //}
                //if (j == 1) {
                //    res[i][j] = 1000 - iscaled;
                //}
                if (res[i][j] < minScaled) {
                    minScaled = res[i][j];
                }
                if (res[i][j] > maxScaled) {
                    maxScaled = res[i][j];
                }

            }
        }
        System.out.println("minScaled: " + minScaled);
        System.out.println("maxScaled: " + maxScaled);

        return res;

    }

    //18.02.14: for reading the aloi data
    public double[][] readFile(String path, int dim, int numObj, String ending) {
        double[][] res = new double[numObj][dim];
        int counter = 0;
        System.out.println("reading file " + path);
        try {
            File source = new File(path);
            StreamTokenizer st1 = new StreamTokenizer(new BufferedReader(new FileReader(path)));
            st1.parseNumbers();
            st1.eolIsSignificant(false);
            while (st1.ttype != StreamTokenizer.TT_EOF) {  //Daten aus File in Vector einlesen
                st1.nextToken();
                //System.out.println(st1.ttype);
                if (st1.ttype == StreamTokenizer.TT_NUMBER) //  System.out.println(st1.nval);
                {
                    if (st1.ttype == StreamTokenizer.TT_WORD) {
                        //System.out.println(st1.sval);
                        //check substring
                        String s = st1.sval;
                        //add coordinates
                        if (s.endsWith(ending)) {
                            for (int i = 0; i < dim; i++) {
                                st1.nextToken();
                                res[counter][i] = st1.nval;

                            }
                        }
                    }
                }

            }
        } catch (IOException e) {            //Es gibt kein solches File
            System.err.println("Error: " + e);
        }
        return res;

    }

    //29.4.14: Scale each row of the data to the precision. MinDist: 10^-p, MaxDist: 1
    public double[][] scaleToPrecision(double[][] data, int p) {
        int numDist = data[0].length;
        int numMetrics = data.length;
        double[][] res = new double[data.length][data[0].length];
        //determine min and max from data
        double[] min = new double[numMetrics];
        double[] max = new double[numMetrics];
        for (int i = 0; i < numMetrics; i++) {
            min[i] = Double.MAX_VALUE;
            max[i] = -Double.MAX_VALUE;
        }
        for (int i = 0; i < numMetrics; i++) {
            for (int j = 0; j < numDist; j++) {
                if (data[i][j] < min[i]) {
                    min[i] = data[i][j];
                }
                if (data[i][j] > max[i]) {
                    max[i] = data[i][j];
                }

            }
        }

        //scale between 0 and 1 and set all distances smaller than minDist to minDist
        double minDist = Math.pow(10, -p);
        for (int i = 0; i < numMetrics; i++) {
            for (int j = 0; j < numDist; j++) {
                res[i][j] = (data[i][j] - min[i]) / max[i];
                res[i][j] = Math.round((res[i][j] * Math.pow(10, p))) / Math.pow(10, p);
                if (res[i][j] < minDist) {
                    res[i][j] = minDist;
                }
            }
        }
        return res;
    }

    public double[][] avoidZeros(double[][] data) {
        double minDist = 0.000001;

        for (int i = 0; i < data.length; i++) {
            int countZeros = 0;
            for (int j = 0; j < data[i].length; j++) {
                if (data[i][j] < minDist) {
                    data[i][j] = minDist;
                    countZeros++;
                }
            }
            System.out.println(i + " " + countZeros);
        }
        return data;
    }

    //scale data in [0,..1]
    public double[][] rowScaleData(double[][] data) {
        int metrics = data.length;
        int dist = data[0].length;
        double[] mins = new double[metrics];
        double[] maxs = new double[metrics];
        for (int i = 0; i < metrics; i++) {
            mins[i] = Double.MAX_VALUE;
            maxs[i] = -Double.MAX_VALUE;
        }
        for (int i = 0; i < metrics; i++) {
            for (int j = 0; j < dist; j++) {
                if (data[i][j] < mins[i]) {
                    mins[i] = data[i][j];
                }
                if (data[i][j] > maxs[i]) {
                    maxs[i] = data[i][j];
                }
            }
        }
        for (int i = 0; i < metrics; i++) {
            for (int j = 0; j < dist; j++) {
                data[i][j] = (data[i][j] - mins[i]) / (maxs[i] - mins[i]);
                //data[i][j] *= 100;
                //data[i][j] += Double.MIN_VALUE;

            }
        }

        return data;
    }

    //scale data in [0,..1]
    public double[][] colScaleData(double[][] data) {
        int rows = data.length;
        int cols = data[0].length;
        double[] mins = new double[cols];
        double[] maxs = new double[cols];
        for (int i = 0; i < cols; i++) {
            mins[i] = Double.MAX_VALUE;
            maxs[i] = -Double.MAX_VALUE;
        }
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (data[i][j] < mins[j]) {
                    mins[j] = data[i][j];
                }
                if (data[i][j] > maxs[j]) {
                    maxs[j] = data[i][j];
                }
            }
        }
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                data[i][j] = (data[i][j] - mins[j]) / (maxs[j] - mins[j]);

            }
        }

        return data;
    }

    public double[][] minMax(double[][] data) {
        int rows = data.length;
        int cols = data[0].length;
        double[][] minMax = new double[cols][2];

        for (int i = 0; i < cols; i++) {
            minMax[i][0] = Double.MAX_VALUE;
            minMax[i][1] = -Double.MAX_VALUE;
        }
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (data[i][j] < minMax[j][0]) {
                    minMax[j][0] = data[i][j];
                }
                if (data[i][j] > minMax[j][1]) {
                    minMax[j][1] = data[i][j];
                }
            }
        }

        return minMax;
    }

    //scale data in [0,..1]
    public void colValueRangeCheck(double[][] data) {
        int rows = data.length;
        int cols = data[0].length;
        double[] mins = new double[cols];
        double[] maxs = new double[cols];
        for (int i = 0; i < cols; i++) {
            mins[i] = Double.MAX_VALUE;
            maxs[i] = -Double.MAX_VALUE;
        }
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (data[i][j] < mins[j]) {
                    mins[j] = data[i][j];
                }
                if (data[i][j] > maxs[j]) {
                    maxs[j] = data[i][j];
                }
            }
        }
        // for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            System.out.println("dim: " + j + " min: " + mins[j] + " max: " + maxs[j]);

        }
        // }

    }

    public int countZeros(double[] data) {
        int counter = 0;
        for (int i = 0; i < data.length; i++) {
            if (data[i] < Double.MIN_VALUE) {
                counter++;
            }
        }
        return counter;

    }

    public double[] correctForZeros(double[] data) {
        double[] res = new double[data.length];
        for (int i = 0; i < data.length; i++) {
            res[i] = data[i] + Double.MIN_VALUE;
        }
        return res;

    }

    //without diagonal
    public int getIndex(int index1, int index2, int numObj) {
        int i1 = Math.min(index1, index2);
        int i2 = Math.max(index1, index2);
        int count = i1;
        int dec = 1;
        int placesBefore = 0;
        while (count > 0) {
            placesBefore = placesBefore + (numObj - dec);
            dec++;
            count--;
        }
        int offset = i2 - i1 - 1;
        return placesBefore + offset;
    }

    public int[] getPair(int index, int numObj) {
        int[] result = new int[2];
        result[0] = 0;
        result[1] = 1;

        for (int i = 0; i < index; i++) {
            if (result[1] < numObj - 1) {
                result[1]++;
            } else {
                result[0]++;
                result[1] = result[0] + 1;
            }
        }

        return result;
    }

    //with diagonal
    public int getIndexDiagonal(int index1, int index2, int numObj) {
        int i1 = Math.min(index1, index2);
        int i2 = Math.max(index1, index2);
        int count = i1;
        int dec = 1;
        int placesBefore = 0;
        while (count > 0) {
            placesBefore = placesBefore + (numObj - dec + 1);
            dec++;
            count--;
        }
        int offset = i2 - i1;
        return placesBefore + offset;
    }

    //    public int getIndex(int index1, int index2, int index3, int numObj) {
//        int[] indices = {index1, index2, index3};
//        Arrays.sort(indices);
//        int first = indices[0];
//        int second = indices[1];
//        int third = indices[2];
//        int placesBefore = 0;
//        int countFirst = 0;
//        int countSecond = 1;
//        int countThird = 2;
//
//        while(countFirst < first && countSecond > second && countThird < third){
//            while (countThird < numObj){
//
//            }
//        }
//        while (countFirst < first){ //count whole blocks
//            int countSecond = countFirst + 1;
//            int countThird = countSecond + 1;
//            while(countSecond < (numObj-1)){
//                placesBefore += numObj - countThird;
//
//            }
//        }
//        while (countFirst < first) {
//            while (numObj - countThird > 0) {
//                placesBefore += (numObj - countThird);
//                countThird++;
//            }
//            countFirst++;
//        }
//        int basic
//        countThird = countFirst + 2;
//        while (countSecond < second) {
//            placesBefore += (numObj - countThird);
//            countThird++;
//            countSecond++;
//        }
//        int offSet = (numObj - 1) - third - 1;
//        return placesBefore + offSet;
    // }
    //21.10.13: check the effect of sorting order of the points for metric coding. Random shuffling.
    public int[] shuffleData(int[] data, int numObj) {
        int[] res = new int[data.length];
        Random r = new Random();
        boolean[] used = new boolean[numObj];
        int[] order = new int[numObj];
        int count = 0;
        while (count < numObj) {
            int next = r.nextInt(numObj);
            while (used[next]) {
                next = r.nextInt(numObj);
            }
            used[next] = true;
            order[count] = next;
            count++;
        }
        count = 0;
        for (int j = 0; j < numObj - 1; j++) {
            for (int i = j + 1; i < numObj; i++) {
                res[count] = data[getIndex(order[j], order[i], numObj)];
                //System.out.println(order[j] + " " + order[i] + " " + getIndex(order[j], order[i], 100));
//            System.out.println(order[j] + " " + order[i] + " " + getIndex(order[j], order[i], 100));
//            System.out.println(order[i] + " " + order[j] + " " + getIndex(order[i], order[j], 100));
                //res[count] = data[getIndex(j, i, numObj)]; //this is correct
                count++;
            }
        }
        return res;
    }

    //5.11.13: reorder data according to specified order
    public int[] reorderData(int[] data, int[] order) {
        int[] res = new int[data.length];
        int numObj = order.length;
        int count = 0;
        for (int j = 0; j < numObj - 1; j++) {
            for (int i = j + 1; i < numObj; i++) {
                res[count] = data[getIndex(order[j], order[i], numObj)];
                //System.out.println(order[j] + " " + order[i] + " " + getIndex(order[j], order[i], 100));
//            System.out.println(order[j] + " " + order[i] + " " + getIndex(order[j], order[i], 100));
//            System.out.println(order[i] + " " + order[j] + " " + getIndex(order[i], order[j], 100));
                //res[count] = data[getIndex(j, i, numObj)]; //this is correct
                count++;
            }
        }
        return res;
    }

    //26.09.2013: create integer values from data for coding. Assumes that data is positive. digits the number of digits used to represent the value range
    public int[][] mapToCodingRange(double[][] data, int digits, boolean[] discrete) {
        // System.out.println("a: " + getIndex(7, 8, 100) + " b: " + getIndex(65, 72, 100));
        int metrics = data.length;
        int dist = data[0].length;
        double min = 1.0;
        double max = Math.pow(10, digits);

        int[][] res = new int[metrics][dist];
        for (int i = 0; i < metrics; i++) {

            double minI = Double.MAX_VALUE;
            double maxI = -Double.MAX_VALUE;
            for (int j = 0; j < dist; j++) {
                if (data[i][j] < minI) {
                    minI = data[i][j];
                }
                if (data[i][j] > maxI) {
                    maxI = data[i][j];
                }
            }
            if (!discrete[i]) {
                for (int j = 0; j < dist; j++) {
                    double sc = ((data[i][j] - minI) / (maxI - minI)) * (max - min) + min;  //scales data between 1 and max but violates the triangle inequality of discrete data
                    int dd = (int) (sc + 0.5);
                    res[i][j] = dd;
                }
            } else {
                double range = maxI - minI;
                double factor = max / range;
                for (int j = 0; j < dist; j++) {
                    //double sc = data[i][j] * factor;
                    double sc = data[i][j];
                    //int dd = (int) (sc + 0.5);
                    int dd = (int) sc;
                    res[i][j] = dd;
                }

            }
        }
        return res;
    }

    //    public boolean checkTriangle(int[][] data, int index){
//        int numObj = data[index].length;
//         for (int j = 1; j < numObj; j++) { //diagonals
//            if (j == 3) {
//                System.out.println("m");
//            }
//            int i = 0;
//            int jj = j;
//            while (i < numObj - j) {
//                //process element
//                int lb = 1;
//                int ub = maxValue;
//                int ja = jj;
//                int ia = jj - 1;
//                int jb = jj - 1;
//                int ib = i;
//                int ind = getIndex(i, jj, numObj);
//                while (jb >= 1 && ia > i) {
//                    double a = data[index][getIndex(ia, ja, numObj)];
//                    double b = data[index][getIndex(ib, jb, numObj)];
//                    double ubN = a + b;
//                    double lbN = Math.max(a, b) - Math.min(a, b);
//                    ub = Math.min(ub, (int) ubN);
//                    lb = Math.max(lb, (int) lbN);
//                    jb--;
//                    ia--;
//                }
//                cc[ind][0] = lg2(ub - lb + 1);
//                if(cc[ind][0] == Double.NEGATIVE_INFINITY)
//                    System.out.println("m");
//                cost += cc[ind][0];
//                i++;
//                jj++;
//            }
//
//        }
//
//    }
    public void rangeCheck(double[][] data) {
        //scale each row of the data between 0 and 1
        int metrics = data.length;
        int dist = data[0].length;
        double[] min = new double[metrics];
        double[] max = new double[metrics];
        for (int i = 0; i < min.length; i++) {
            min[i] = Double.MAX_VALUE;
            max[i] = -Double.MAX_VALUE;
        }
        for (int i = 0; i < metrics; i++) {
            for (int j = 0; j < dist; j++) {
                if (data[i][j] < min[i]) {
                    min[i] = data[i][j];
                }
                if (data[i][j] > max[i]) {
                    max[i] = data[i][j];
                }
            }
        }
        for (int i = 0; i < metrics; i++) {
            System.out.println(i + " min: " + min[i] + " max " + max[i]);
        }
    }

    public double[][] scaleLargestAxis(double[][] coord) {
        double max_x = -java.lang.Double.MAX_VALUE;
        double max_y = -java.lang.Double.MAX_VALUE;
        double min_x = -max_x;
        double min_y = -max_y;
        int numObj = coord.length;
        double[][] coordv = new double[coord.length][coord[0].length];
        for (int i = 0; i < numObj; i++) {

            if (coord[i][0] > max_x) {
                max_x = coord[i][0];
            }
            if (coord[i][0] < min_x) {
                min_x = coord[i][0];
            }
            if (coord[i][1] > max_y) {
                max_y = coord[i][1];
            }
            if (coord[i][1] < min_y) {
                min_y = coord[i][1];
            }
        }
        if (max_x - min_x > max_y - min_y) {
            for (int i = 0; i < numObj; i++) {
                coordv[i][0] = ((coord[i][0] - min_x) / (max_x - min_x));
                coordv[i][1] = ((coord[i][1] - min_y) / (max_x - min_x));
            }
        } else {
            for (int i = 0; i < numObj; i++) {
                coordv[i][0] = ((coord[i][0] - min_x) / (max_y - min_y));
                coordv[i][1] = ((coord[i][1] - min_y) / (max_y - min_y));
            }
        }
        return coordv;

    }

    public double[] scale(double[] data) {
        //scale each row of the data between 0 and 1
        int metrics = data.length;

        double rowSum = 0.0;

//        //avoid zeros
//        //double minValue = 1E-6;
//        for (int i = 0; i < metrics; i++) {
//            for (int j = 0; j < dist; j++) {
//                if (data[i][j] < minValue) {
//                    data[i][j] = minValue;
//                }
//            }
//        }
        double min = Double.MAX_VALUE;
        double max = -Double.MAX_VALUE;
        for (int i = 0; i < metrics; i++) {

            //for (int j = 0; j < dist; j++) {
            if (data[i] < min) {
                min = data[i];
            }
            if (data[i] > max) {
                max = data[i];
            }

        }

        for (int i = 0; i < data.length; i++) {
            // for (int j = 0; j < dist; j++) {
            data[i] = ((data[i] - min) / (max - min));

        }

        //}
        //        //adjust to a row sum of 100
//        for (int i = 0; i < metrics; i++) {
//            double weightFactor = 1.0 / (rowSum[i] / 10000.0);
//            for (int j = 0; j < dist; j++) {
//                data[i][j] = data[i][j] * weightFactor;
//            }
//        }
//
        return data;
    }

    public double[][] scaleCoordinates(double[][] data) {
        int numObj = data.length;
        int dim = data[0].length;
        double[] min = new double[dim];
        double[] max = new double[dim];
        for (int i = 0; i < dim; i++) {
            min[i] = Double.MAX_VALUE;
            max[i] = -Double.MAX_VALUE;
        }
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < dim; j++) {
                if (data[i][j] < min[j]) {
                    min[j] = data[i][j];
                }
                if (data[i][j] > max[j]) {
                    max[j] = data[i][j];
                }
            }
        }
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < dim; j++) {
                data[i][j] = ((data[i][j] - min[j]) / (max[j] - min[j]));
            }
        }
        return data;
    }

    public double[][] scaleAndNormalize(double[][] data) {
        //scale each row of the data between 0 and 1
        int metrics = data.length;
        int dist = data[0].length;
        double[] rowSum = new double[metrics];

//        //avoid zeros
//        //double minValue = 1E-6;
//        for (int i = 0; i < metrics; i++) {
//            for (int j = 0; j < dist; j++) {
//                if (data[i][j] < minValue) {
//                    data[i][j] = minValue;
//                }
//            }
//        }
        for (int i = 0; i < metrics; i++) {
            double min = Double.MAX_VALUE;
            double max = -Double.MAX_VALUE;
            for (int j = 0; j < dist; j++) {
                if (data[i][j] < min) {
                    min = data[i][j];
                }
                if (data[i][j] > max) {
                    max = data[i][j];
                }
            }
            for (int j = 0; j < dist; j++) {
                data[i][j] = ((data[i][j] + min) / (max + min));
                rowSum[i] += data[i][j];
            }
        }
//        //adjust to a row sum of 100
//        for (int i = 0; i < metrics; i++) {
//            double weightFactor = 1.0 / (rowSum[i] / 10000.0);
//            for (int j = 0; j < dist; j++) {
//                data[i][j] = data[i][j] * weightFactor;
//            }
//        }
//
        return data;
    }

    public double mCheck(double[] dist, int numObj) {
        double maxVio = Double.MAX_VALUE;
        for (int i = 0; i <= numObj - 3; i++) {
            for (int j = i + 1; j <= numObj - 2; j++) {
                for (int kk = j + 1; kk <= numObj - 1; kk++) {
                    // int ijk = du.getIndex3(i, j, kk, numObj);
                    int ki = getIndex(kk, i, numObj);
                    int jk = getIndex(j, kk, numObj);
                    int ij = getIndex(i, j, numObj);
                    //using ij as long side
                    double diff_ij = dist[ki] + dist[jk] - dist[ij];
                    if (diff_ij < maxVio) {
                        maxVio = diff_ij;
                    }
                    double diff_jk = dist[ki] + dist[ij] - dist[jk];
                    if (diff_jk < maxVio) {
                        maxVio = diff_jk;
                    }
                    double diff_ki = dist[ij] + dist[jk] - dist[ki];
                    if (diff_ki < maxVio) {
                        maxVio = diff_ki;
                    }
                }
            }
        }
        return maxVio;
    }

    public boolean metricCheck(double[] dist, int numObj) {
        double maxVio = 0.0;
        boolean v = true;
        if (v) {
            System.out.println("----------------metric check--------------");
        }
        boolean passed = true;
        int counter_v = 0;
        int counter_t = 0;
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < numObj; j++) {
                for (int k = 0; k < numObj; k++) {
                    if (k != i && k != j && i != j) {

                        int index_ij = getIndex(i, j, numObj);
                        int index_ik = getIndex(i, k, numObj);
                        int index_kj = getIndex(k, j, numObj);
                        double ij = dist[index_ij];
                        double ik = dist[index_ik];
                        double kj = dist[index_kj];

//                        int index_ki = getIndex(k, i, numObj);
//                        int index_jk = getIndex(j, k, numObj);
//                        int index_ij = getIndex(i, j, numObj);
//
//                        double ij = dist[index_ij];
//                        double ki = dist[index_ki];
//                        double jk = dist[index_jk];
                        double upperbound_ij = ik + kj;
                        double lowerbound_ij = Math.abs(ik - kj);

                        if (ij > upperbound_ij) {
                            double diff = ij - upperbound_ij;
                            if (diff > 0) {
                                passed = false;
                                //System.out.println("ij > upperbound: diff: " + diff + " index_ij: " + index_ij + " index_ik: " + index_ik + " index_kj: " + index_kj);
                                counter_v++;

                            }

                        }
                        if (ij < lowerbound_ij) {
                            double diff = lowerbound_ij - ij;
                            if (diff > 0) {
                                passed = false;
                                if (diff > maxVio) {
                                    maxVio = diff;
                                }
                                //System.out.println("ij < lowerbound: diff: " + diff + " index_ij: " + index_ij + " index_ik: " + index_ik + " index_kj: " + index_kj);
                                counter_v++;
                            }

                        }
                        counter_t++;
                    }
                }
            }
        }
        if (v) {
            System.out.println("metric check passed: " + passed + " number of violations: " + counter_v + " of " + counter_t + " triangle inequalities in total. MaxVio: " + maxVio);
        }
        return passed;
    }


    public void saveAsMatlab3(double[][] d, double[][] e, double[][] t, String vName_d, String vName_e, String vName_t, String fName) {
        MLDouble v = new MLDouble(vName_d, d);
        MLDouble v1 = new MLDouble(vName_e, e);
        MLDouble v2 = new MLDouble(vName_t, t);
        ArrayList ll = new ArrayList();
        ll.add(v);
        ll.add(v1);
        ll.add(v2);
        MatFileWriter mw = new MatFileWriter();
        try {
            mw.write(fName, ll);
        } catch (IOException ex) {
        }
        //System.out.println("Matlab matrix saved.");
    }

    public void saveAsMatlab4(double[][] d, double[][] e, double[][] t, double[][] h, String vName_d, String vName_e, String vName_t, String vName_h, String fName) {
        MLDouble v = new MLDouble(vName_d, d);
        MLDouble v1 = new MLDouble(vName_e, e);
        MLDouble v2 = new MLDouble(vName_t, t);
        MLDouble v3 = new MLDouble(vName_h, h);
        ArrayList ll = new ArrayList();
        ll.add(v);
        ll.add(v1);
        ll.add(v2);
        ll.add(v3);
        MatFileWriter mw = new MatFileWriter();
        try {
            mw.write(fName, ll);
        } catch (IOException ex) {
        }
        //System.out.println("Matlab matrix saved.");
    }

    public void saveAsMatlab5(double[][] d, double[][] e, double[][] t, double[][] h, double[][] k, String vName_d, String vName_e, String vName_t, String vName_h, String vName_k, String fName) {
        MLDouble v = new MLDouble(vName_d, d);
        MLDouble v1 = new MLDouble(vName_e, e);
        MLDouble v2 = new MLDouble(vName_t, t);
        MLDouble v3 = new MLDouble(vName_h, h);
        MLDouble v4 = new MLDouble(vName_k, k);
        ArrayList ll = new ArrayList();
        ll.add(v);
        ll.add(v1);
        ll.add(v2);
        ll.add(v3);
        ll.add(v4);
        MatFileWriter mw = new MatFileWriter();
        try {
            mw.write(fName, ll);
        } catch (IOException ex) {
        }
        //System.out.println("Matlab matrix saved.");
    }

    public void saveAsMatlab6(double[][] d, double[][] e, double[][] t, double[][] h, double[][] k, double[][] m, String vName_d, String vName_e, String vName_t, String vName_h, String vName_k, String vName_m, String fName) {
        MLDouble v = new MLDouble(vName_d, d);
        MLDouble v1 = new MLDouble(vName_e, e);
        MLDouble v2 = new MLDouble(vName_t, t);
        MLDouble v3 = new MLDouble(vName_h, h);
        MLDouble v4 = new MLDouble(vName_k, k);
        MLDouble v5 = new MLDouble(vName_m, m);
        ArrayList ll = new ArrayList();
        ll.add(v);
        ll.add(v1);
        ll.add(v2);
        ll.add(v3);
        ll.add(v4);
        ll.add(v5);
        MatFileWriter mw = new MatFileWriter();
        try {
            mw.write(fName, ll);
        } catch (IOException ex) {
        }
        //System.out.println("Matlab matrix saved.");
    }


    public void saveResult(double[][] d, String fName) {
        FileOutputStream output = null;
        try {
            File target = new File(fName);
            output = new FileOutputStream(target);
            PrintStream writer = new PrintStream(output);


            for (int i = 0; i < d.length; i++) {
                for (int j = 0; j < d[i].length - 1; j++)
                    writer.print(d[i][j] + " ");
                writer.println(d[i][d[i].length - 1]);

            }


        } catch (Exception ex) {
            System.out.println(ex.getMessage());
//             Logger.getLogger(DataUtils.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                output.close();
            } catch (IOException ex) {
                System.out.println(ex.getMessage());
//                   Logger.getLogger(DataUtils.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

    }

    public void saveAsMatlab(double[][] d, String vName, String fName) {
        MLDouble v = new MLDouble(vName, d);
        ArrayList ll = new ArrayList();
        ll.add(v);
        MatFileWriter mw = new MatFileWriter();
        try {
            mw.write(fName, ll);
        } catch (IOException ex) {
        }
        //System.out.println("Matlab matrix saved.");
    }


}
