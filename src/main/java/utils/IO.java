/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package utils;

import com.jmatio.io.MatFileReader;
import com.jmatio.io.MatFileWriter;
import com.jmatio.types.MLDouble;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import edu.uci.ics.jung.graph.util.Pair;
import sun.rmi.server.InactiveGroupException;

import java.io.*;
import java.util.ArrayList;
import java.util.Vector;


/**
 * @author claudia
 */
public class IO {

    public Graph matlabToGraph(String dir, String variableName) {
        Graph<Integer, Integer> g = new UndirectedSparseGraph<Integer, Integer>();
        DataUtils du = new DataUtils();


        MatFileReader mfr = null;
        try {
            mfr = new MatFileReader(dir);
        } catch (IOException e) {
        }

        if (mfr != null) {
            double[][] data = ((MLDouble) mfr.getMLArray(variableName)).getArray();
            for (int i = 0; i < data.length; i++) {
                g.addVertex(i);

            }
            int counter = 0;
            for (int i = 0; i < data.length; i++) {
                for (int j = 0; j < data.length; j++) {
                    if (i > j && data[i][j] == 1.0) {
                        //g.addEdge(du.getIndex(i, j, data.length), i, j);
                        g.addEdge(counter, i, j);
                        counter++;
                    }
                }
            }

        }
        return g;
    }

    public void writeGraphToMatlabSorted(Graph<Integer, Integer> g, int[] labels, String filename) {
        int numObj = g.getVertexCount();
        double[][] dist = new double[numObj][numObj];
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < numObj; j++) {
                if (g.isNeighbor(i, j)) {
                    dist[i][j] = 1.0;
                } else {
                    dist[i][j] = 0.0;
                }
            }
        }
        MLDouble q = new MLDouble("graph", dist);
        ArrayList ll = new ArrayList();
        ll.add(q);
        MatFileWriter mw = new MatFileWriter();
        try {
            String name = filename + ".mat";
            mw.write(name, ll);
        } catch (IOException ex) {
            //     Logger.getLogger(VI.class.getName()).log(Level.SEVERE, null, ex);
        }

    }


    public Graph<Integer, Integer> readEdgeList(String filename) {
        Graph<Integer, Integer> g = new UndirectedSparseGraph<Integer, Integer>();
        File nodeFile = new File(filename);
        if (nodeFile.isFile()) {
            try {
                BufferedReader csvReader = new BufferedReader(new FileReader(nodeFile));
                StreamTokenizer st = new StreamTokenizer(csvReader);
                String row;
             st.nextToken();
                int numNodes = (int) st.nval;
                for (int i = 0; i < numNodes; i++)
                    g.addVertex(i);
                int edgeCounter = 0;
                while (st.nextToken() != StreamTokenizer.TT_EOF) {
                    int n1 = (int) st.nval;
                    st.nextToken();
                    int n2 = (int) st.nval;
                    g.addEdge(edgeCounter, n1, n2);
                    edgeCounter++;
                }


            } catch (IOException e) {
            }
        }
        return g;

    }

    public void writeGraphToMatlab(Graph<Integer, Integer> g, String filename) {
        int numObj = g.getVertexCount();
        double[][] dist = new double[numObj][numObj];
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < numObj; j++) {
                if (g.isNeighbor(i, j)) {
                    dist[i][j] = 1.0;
                } else {
                    dist[i][j] = 0.0;
                }
            }
        }
        MLDouble q = new MLDouble("graph", dist);
        ArrayList ll = new ArrayList();
        ll.add(q);
        MatFileWriter mw = new MatFileWriter();
        try {
            String name = filename + ".mat";
            mw.write(name, ll);
        } catch (IOException ex) {
            //     Logger.getLogger(VI.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    public void writeFigure() {
//        double[][] gt = readMatlabMatrix("pr.mat", "gtCoordRed");
//        double[][] embedding = readMatlabMatrix("pr.mat", "Z_OwnBlowup");


        double[][] gt = readMatlabMatrix("twoMoonsProcrustes.mat", "groundTruth");
        double[][] embedding = readMatlabMatrix("twoMoonsProcrustes.mat", "speP");


        String fn = "resultSPE.txt";
        try {
            FileOutputStream fout = new FileOutputStream(new File(fn));
            for (int i = 0; i < gt.length; i++) {
                for (int j = 0; j < gt[i].length; j++) {
                    fout.write((gt[i][j] + " ").getBytes());
                }
                fout.write(("\n").getBytes());
            }
            fout.write(("\n").getBytes());
            fout.write(("\n").getBytes());
            for (int i = 0; i < embedding.length; i++) {
                for (int j = 0; j < embedding[i].length; j++) {
                    fout.write((embedding[i][j] + " ").getBytes());
                }
                fout.write(("\n").getBytes());
            }
            fout.write(("\n").getBytes());
            fout.write(("\n").getBytes());
            for (int i = 0; i < gt.length; i++) {
                for (int j = 0; j < gt[i].length; j++) {
                    fout.write((gt[i][j] + " ").getBytes());
                }
                fout.write(("\n").getBytes());
                for (int j = 0; j < embedding[i].length; j++) {
                    fout.write((embedding[i][j] + " ").getBytes());
                }
                fout.write(("\n").getBytes());
                fout.write(("\n").getBytes());

            }

            fout.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }


    }

    public void writeDoubleToMatlab(double[][] d, String variablename, String filename) {
        MLDouble q = new MLDouble(variablename, d);
        ArrayList ll = new ArrayList();
        ll.add(q);
        MatFileWriter mw = new MatFileWriter();
        try {
            String name = filename + ".mat";
            mw.write(name, ll);
        } catch (IOException ex) {
            //     Logger.getLogger(VI.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void writeDoubleToMatlab(double[][] d, String filename) {
        MLDouble q = new MLDouble(filename, d);
        ArrayList ll = new ArrayList();
        ll.add(q);
        MatFileWriter mw = new MatFileWriter();
        try {
            String name = filename + ".mat";
            mw.write(name, ll);
        } catch (IOException ex) {
            //     Logger.getLogger(VI.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public double[][] readMatlabMatrix(String dir, String variableName) {
        double[][] data = new double[1][1];

        MatFileReader mfr = null;
        try {
            mfr = new MatFileReader(dir);
        } catch (IOException e) {
        }

        if (mfr != null) {
            data = ((MLDouble) mfr.getMLArray(variableName)).getArray();
        }
        return data;
    }


}
