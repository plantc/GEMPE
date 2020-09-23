/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package visu;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import com.jmatio.io.MatFileWriter;
import com.jmatio.types.MLDouble;
import edu.uci.ics.jung.algorithms.layout.CircleLayout;
import edu.uci.ics.jung.algorithms.layout.FRLayout;
import edu.uci.ics.jung.algorithms.layout.GraphElementAccessor;
import edu.uci.ics.jung.algorithms.layout.ISOMLayout;
import edu.uci.ics.jung.algorithms.layout.KKLayout;
import edu.uci.ics.jung.algorithms.layout.SpringLayout;
import edu.uci.ics.jung.algorithms.layout.StaticLayout;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraDistance;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraShortestPath;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.util.Context;
import edu.uci.ics.jung.visualization.BasicVisualizationServer;
import edu.uci.ics.jung.visualization.VisualizationViewer;
import edu.uci.ics.jung.visualization.control.AbstractPopupGraphMousePlugin;
import edu.uci.ics.jung.visualization.control.DefaultModalGraphMouse;
import edu.uci.ics.jung.visualization.control.GraphMouseListener;
import edu.uci.ics.jung.visualization.decorators.AbstractEdgeShapeTransformer;
import edu.uci.ics.jung.visualization.decorators.EdgeShape;
import edu.uci.ics.jung.visualization.decorators.EllipseVertexShapeTransformer;
import edu.uci.ics.jung.visualization.decorators.ToStringLabeller;
import edu.uci.ics.jung.visualization.renderers.Renderer.VertexLabel.Position;
//import org.jblas.*;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Paint;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.geom.Point2D;
import java.awt.geom.AffineTransform;
import java.awt.event.ActionEvent;
import java.awt.geom.Ellipse2D;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.AbstractAction;
import javax.swing.JFrame;
import javax.swing.JPopupMenu;
import org.apache.commons.collections15.Transformer;
import org.apache.commons.collections15.TransformerUtils;
import utils.DataUtils;
import utils.IO;
import nature.GraphCompression;
//import static org.jblas.Eigen.symmetricGeneralizedEigenvectors;

/**
 *
 * @author
 */
public class Visualization {

    Graph g;
    int numObj;
    double maxDist = 30;
    double[] costs;

    public Visualization(Graph g) {
        this.g = g;
        numObj = g.getVertexCount();
        // The Layout<V, E> is parameterized by the vertex and edge types
//        l = new KKLayout(g);
//        l.setSize(new Dimension(300, 300)); // sets the initial size of the space
//// The BasicVisualizationServer<V,E> is parameterized by the edge types
//        BasicVisualizationServer<Integer, Integer> vv = new BasicVisualizationServer<Integer, Integer>(l);
//        vv.setPreferredSize(new Dimension(350, 350)); //Sets the viewing area size
//        JFrame frame = new JFrame("Simple Graph View");
//        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//        frame.getContentPane().add(vv);
//        frame.pack();
//        frame.setVisible(true);
//        for (Integer v : l.getGraph().getVertices()) {
//            Point2D p = l.transform(v);
//            System.out.println("m");

//        }
    }











    public void scaleToDisplaySize() {
    }





    public synchronized double[][] getCoordinatesStochasticGD() {
        //double[][] coord = new double[g.getVertexCount()][2];
//        Embedding e = new Embedding(g,2);
//        Point2D[] pp = e.getIsoCoord();
        //DisplGrid wi = new DisplGrid(g); //current version
        StochasticOwnDispl wi = new StochasticOwnDispl(g);
        //IncrementalDispl wi = new IncrementalDispl(g);
        //StochasticLockfreeSimple wi = new StochasticLockfreeSimple(g);
        //IncrementalOwnDispl wi = new IncrementalOwnDispl(g); //funktioniert, aber Verfaltungen
        //ActiveStochasticOwnDispl wi = new ActiveStochasticOwnDispl(g);

        //WeightedMajorizationDispl wi = new WeightedMajorizationDispl(g);
//        DataUtils du = new DataUtils();
//        double[][] gt = du.readMatlabMatrix("init.mat", "init");
//        wi.setGroundTruth(gt);
        //HierarchicalWeighted wi = new HierarchicalWeighted(g);
        wi.setSize(new Dimension(600, 600)); // sets the initial size of the space
        long startTime = System.currentTimeMillis();
        wi.initialize();
        //wi.InitializeOld();
//
        VisualizationViewer<Integer, Integer> vv = new VisualizationViewer<Integer, Integer>(wi);
        vv.setBackground(Color.WHITE);

        vv.setPreferredSize(new Dimension(650, 650)); //Sets the viewing area size
        EllipseVertexShapeTransformer<Integer> t = new EllipseVertexShapeTransformer<Integer>();

        t.setSizeTransformer(new Transformer<Integer, Integer>() {
            public Integer transform(Integer i) {
                return 6;
            }
        });
        AbstractEdgeShapeTransformer<Integer, Integer> ast = new EdgeShape.Line<Integer, Integer>();

        vv.getRenderContext().setVertexShapeTransformer(t);
        vv.getRenderContext().setEdgeShapeTransformer(ast);
        // vv.getRenderContext().setVertexFillPaintTransformer(new MyVertexFillPaintFunction(clid));
//        DefaultModalGraphMouse gm = new DefaultModalGraphMouse<Integer, Integer>();
//        vv.setGraphMouse(gm);
//        PopupGraphMousePlugin p = new PopupGraphMousePlugin();
//        p.setNames(names);
//        gm.add(p);
        double ll = wi.getBestllh();
        String tt = java.lang.Double.toString(ll);
        JFrame frame = new JFrame();
        frame.setTitle("");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().add(vv);
        frame.pack();
        frame.setVisible(true);
//        while (!wi.done()) {
////            String s = l.getStatus();
////            System.out.println(s);
////            double[][] cc = l.getCoordinates();
////            GraphCompression gcc = new GraphCompression(g, cc, 30);
////            System.out.println(gcc.mdlFunction());
//            try {
//                wait(100);
//            } catch (InterruptedException ex) {
//                Logger.getLogger(Visualization.class.getName()).log(Level.SEVERE, null, ex);
//            }
//
//        }

        while (!wi.done()) {
//            String s = l.getStatus();
//            System.out.println(s);
//            double[][] cc = l.getCoordinates();
//            GraphCompression gcc = new GraphCompression(g, cc, 30);
//            System.out.println(gcc.mdlFunction());
            try {
                wait(100);
            } catch (InterruptedException ex) {
                Logger.getLogger(Visualization.class.getName()).log(Level.SEVERE, null, ex);
            }

        }
        long endTime = System.currentTimeMillis();
        double runtime = (endTime - startTime) / 1000;
        System.out.println("runtime: " + runtime);
        double[][] cd = wi.getCoordinates();
        DataUtils du = new DataUtils();
        du.saveResult(cd, "result.txt");
       // du.saveAsMatlab(cd, "result", "result.mat");

        // displayCoord(cd, "result");
//        GraphCompression gcd = new GraphCompression(g, cd, 30);
//        System.out.println("final: " + gcd.mdlFunction());
//
//
//         try {
//            wait(500000000);
//        } catch (InterruptedException ex) {
//            Logger.getLogger(Visualization.class.getName()).log(Level.SEVERE, null, ex);
//        }
//
//        // }
//        for (Integer v : wi.getGraph().getVertices()) {
//            Point2D p = wi.transform(v);
//            coord[v][0] = p.getX();
//            coord[v][1] = p.getY();
//        }
        return cd;

    }

    public synchronized void displayMovie(double[][] gtCoord) {
        //double[][] coord = new double[g.getVertexCount()][2];
//        Embedding e = new Embedding(g,2);
//        Point2D[] pp = e.getIsoCoord();
        //DisplGrid wi = new DisplGrid(g); //current version
        MovieDisplay wi = new MovieDisplay(g);
        //IncrementalDispl wi = new IncrementalDispl(g);
        //StochasticLockfreeSimple wi = new StochasticLockfreeSimple(g);
        //IncrementalOwnDispl wi = new IncrementalOwnDispl(g); //funktioniert, aber Verfaltungen
        //ActiveStochasticOwnDispl wi = new ActiveStochasticOwnDispl(g);

        //WeightedMajorizationDispl wi = new WeightedMajorizationDispl(g);
//        DataUtils du = new DataUtils();
//        double[][] gt = du.readMatlabMatrix("init.mat", "init");
//        wi.setGroundTruth(gt);
        //HierarchicalWeighted wi = new HierarchicalWeighted(g);
        DataUtils du = new DataUtils();
        double[][] groundTruth = du.readMatlabMatrix("coordi_0.mat", "coord");
        wi.setGroundTruth(groundTruth);

        wi.setSize(new Dimension(600, 600)); // sets the initial size of the space
        long startTime = System.currentTimeMillis();

        wi.initialize();
        //wi.InitializeOld();
//
        VisualizationViewer<Integer, Integer> vv = new VisualizationViewer<Integer, Integer>(wi);

        vv.setPreferredSize(new Dimension(650, 650)); //Sets the viewing area size
        EllipseVertexShapeTransformer<Integer> t = new EllipseVertexShapeTransformer<Integer>();

        t.setSizeTransformer(new Transformer<Integer, Integer>() {
            public Integer transform(Integer i) {
                return 6;
            }
        });
        AbstractEdgeShapeTransformer<Integer, Integer> ast = new EdgeShape.Line<Integer, Integer>();

        vv.getRenderContext().setVertexShapeTransformer(t);
        vv.getRenderContext().setEdgeShapeTransformer(ast);
        vv.getRenderContext().setVertexFillPaintTransformer(new MyVertexFillPaintFunctionGeo(gtCoord));
        vv.setBackground(Color.WHITE);
        // vv.getRenderContext().setVertexFillPaintTransformer(new MyVertexFillPaintFunction(clid));
//        DefaultModalGraphMouse gm = new DefaultModalGraphMouse<Integer, Integer>();
//        vv.setGraphMouse(gm);
//        PopupGraphMousePlugin p = new PopupGraphMousePlugin();
//        p.setNames(names);
//        gm.add(p);
        //double ll = wi.getBestllh();
        //String tt = java.lang.Double.toString(ll);
        JFrame frame = new JFrame();
        frame.setTitle("");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().add(vv);
        frame.pack();
        frame.setVisible(true);

//        while (!wi.done()) {
////            String s = l.getStatus();
////            System.out.println(s);
////            double[][] cc = l.getCoordinates();
////            GraphCompression gcc = new GraphCompression(g, cc, 30);
////            System.out.println(gcc.mdlFunction());
//            try {
//                wait(100);
//            } catch (InterruptedException ex) {
//                Logger.getLogger(Visualization.class.getName()).log(Level.SEVERE, null, ex);
//            }
//
//        }
        while (!wi.done()) {
//            String s = l.getStatus();
//            System.out.println(s);
//            double[][] cc = l.getCoordinates();
//            GraphCompression gcc = new GraphCompression(g, cc, 30);
//            System.out.println(gcc.mdlFunction());
            try {
                wait(100);
            } catch (InterruptedException ex) {
                Logger.getLogger(Visualization.class.getName()).log(Level.SEVERE, null, ex);
            }

        }
        long endTime = System.currentTimeMillis();
        double runtime = (endTime - startTime) / 1000;

        // displayCoord(cd, "result");
//        GraphCompression gcd = new GraphCompression(g, cd, 30);
//        System.out.println("final: " + gcd.mdlFunction());
//
//
//         try {
//            wait(500000000);
//        } catch (InterruptedException ex) {
//            Logger.getLogger(Visualization.class.getName()).log(Level.SEVERE, null, ex);
//        }
//
//        // }
//        for (Integer v : wi.getGraph().getVertices()) {
//            Point2D p = wi.transform(v);
//            coord[v][0] = p.getX();
//            coord[v][1] = p.getY();
//        }
    }



    public double[][] getLaplacianCoord(int dim) {
        double[][] coord = new double[numObj][dim];
        double[][] w = new double[numObj][numObj];
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < numObj; j++) {
                if (g.isNeighbor(i, j)) {
                    w[i][j] = 1.0;
                } else {
                    w[i][j] = 0.0;
                }
            }
        }
        double[][] d = new double[numObj][numObj];
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < numObj; j++) {
                d[i][i] += w[i][j];
            }
        }

//        //laplacian eigenmaps
        Matrix D = new Matrix(d);
        Matrix LL = D.minus(new Matrix(w));
        //debug generalized eigenvalue problem: write out LL and D
        IO ea = new IO();
        ea.writeDoubleToMatlab(d, "D");
        ea.writeDoubleToMatlab(LL.getArrayCopy(), "L");

        double[][] dd = new double[numObj][numObj];
        for (int i = 0; i < numObj; i++) {
            dd[i][i] = 1.0 / d[i][i];
        }
        Matrix DD = new Matrix(dd);
        Matrix Ln = DD.times(LL);

        //EigenvalueDecomposition eig = new EigenvalueDecomposition(Ln); //normalized
        EigenvalueDecomposition eig = new EigenvalueDecomposition(LL); //not normalized
        Matrix vectors = eig.getV();
        Matrix values = eig.getD();

        //check scaling constraint
        Matrix cc = vectors.times(D).times(vectors.transpose());
        System.out.println("m");

        //debug
        ea.writeDoubleToMatlab(vectors.getArrayCopy(), "vectors");
        ea.writeDoubleToMatlab(values.getArrayCopy(), "values");
        int nonZeroIndex = 1;
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < dim; j++) {
                //2 largest eigenvalues
//            double x = vectors.get((vectors.getRowDimension()-1), i) * values.get(values.getRowDimension()-1, values.getRowDimension()-1);
//            double y = vectors.get((vectors.getRowDimension()-2), i) * values.get(values.getRowDimension()-2, values.getRowDimension()-2);

                //2 smallest eigenvalues > 0
                coord[i][j] = vectors.get(i, nonZeroIndex + j);
            }
        }
        return coord;

    }



    private double[] randomWalk(int numSteps) {
        double[] res = new double[numObj];
        Random r = new Random();
        int start = r.nextInt(numObj);
        res[start] += 1.0 / numSteps;
        Integer curr = new Integer(start);
        for (int i = 0; i < numSteps; i++) {
            Collection<Integer> successors = g.getSuccessors(curr);
            Object[] s = successors.toArray();
            int nextIndex = r.nextInt(s.length);
            Integer next = (Integer) s[nextIndex];
            res[next.intValue()] += 1.0 / numSteps;
            curr = next;
        }
        return res;
    }

    private double[] randomBroadWalk(int numSteps) {
        double[] res = new double[numObj];
        Random r = new Random();
        int start = r.nextInt(numObj);
        res[start] += 1.0 / numSteps;
        Integer curr = new Integer(start);
        for (int i = 0; i < numSteps; i++) {
            Collection<Integer> successors = g.getSuccessors(curr);
            Object[] s = successors.toArray();
            int nextIndex = r.nextInt(s.length);
            Integer next = (Integer) s[nextIndex];
            for (int j = 0; j < s.length; j++) {
                Integer b = (Integer) s[j];
                res[b.intValue()] += 1.0 / numSteps;

            }
            curr = next;
        }
        return res;
    }

    protected class PopupGraphMousePlugin extends AbstractPopupGraphMousePlugin
            implements MouseListener {

        public PopupGraphMousePlugin() {
            this(MouseEvent.BUTTON3_MASK);
        }

        public PopupGraphMousePlugin(int modifiers) {
            super(modifiers);
        }
        String[] names;

        public void setNames(String[] names) {
            this.names = names;
        }

        /**
         * If this event is over a Vertex, pop up a menu to allow the user to
         * increase/decrease the voltage attribute of this Vertex
         *
         * @param e
         */
        @SuppressWarnings("unchecked")
        protected void handlePopup(MouseEvent e) {
            final VisualizationViewer<Integer, Number> vv
                    = (VisualizationViewer<Integer, Number>) e.getSource();
            Point2D p = e.getPoint();//vv.getRenderContext().getBasicTransformer().inverseViewTransform(e.getPoint());

            GraphElementAccessor<Integer, Number> pickSupport = vv.getPickSupport();
            if (pickSupport != null) {
                final Integer v = pickSupport.getVertex(vv.getGraphLayout(), p.getX(), p.getY());
                if (v != null) {
                    JPopupMenu popup = new JPopupMenu();
                    popup.add(new AbstractAction("Display Number") {
                        public void actionPerformed(ActionEvent e) {
                            System.out.println(names[v.intValue()]);

                        }
                    });
                    popup.add(new AbstractAction("Bla") {
                        public void actionPerformed(ActionEvent e) {
//                            Double value = Math.max(0,
//                                    transparency.get(v).doubleValue() - 0.1);
//                            transparency.put(v, value);
//                            MutableDouble value = (MutableDouble)transparency.getNumber(v);
//                            value.setDoubleValue(Math.max(0, value.doubleValue() - 0.1));
                            //vv.repaint();
                        }
                    });
                    popup.show(vv, e.getX(), e.getY());
                } else {
                    final Number edge = pickSupport.getEdge(vv.getGraphLayout(), p.getX(), p.getY());
                    if (edge != null) {
                        JPopupMenu popup = new JPopupMenu();
                        popup.add(new AbstractAction(edge.toString()) {
                            public void actionPerformed(ActionEvent e) {
                                System.err.println("got " + edge);
                            }
                        });
                        popup.show(vv, e.getX(), e.getY());

                    }
                }
            }
        }
    }

    public synchronized double[][] getCoordinatesFR() {
        double[][] coord = new double[g.getVertexCount()][2];
        FRLayout<Integer, Integer> l = new FRLayout(g);
        l.setMaxIterations(5000);
        l.setSize(new Dimension(300, 300)); // sets the initial size of the space
// The BasicVisualizationServer<V,E> is parameterized by the edge types
        BasicVisualizationServer<Integer, Integer> vv = new BasicVisualizationServer<Integer, Integer>(l);
        vv.setPreferredSize(new Dimension(350, 350)); //Sets the viewing area size

        EllipseVertexShapeTransformer<Integer> t = new EllipseVertexShapeTransformer<Integer>();

        t.setSizeTransformer(new Transformer<Integer, Integer>() {
            public Integer transform(Integer i) {
                return 6;
            }
        });
        vv.getRenderContext().setVertexShapeTransformer(t);
        JFrame frame = new JFrame("FR");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().add(vv);
        frame.pack();
        frame.setVisible(true);
//

        while (!l.done()) {
//            String s = l.getStatus();
//            System.out.println(s);

            try {
                wait(100);
            } catch (InterruptedException ex) {
                Logger.getLogger(Visualization.class.getName()).log(Level.SEVERE, null, ex);
            }

        }
        for (Integer v : l.getGraph().getVertices()) {
            Point2D p = l.transform(v);
            coord[v][0] = p.getX();
            coord[v][1] = p.getY();
        }
        GraphCompression gk = new GraphCompression(g, coord);
        double costWithoutEmbedding = gk.codingCostNoEmbedding();
        System.out.println("costWithoutEmbedding: " + costWithoutEmbedding);
        double compressionCost = gk.mdlFunctionSimpleSigmoid();
        System.out.println("compression cost: " + compressionCost);
        return coord;

    }

    public synchronized double[][] getCoordinatesFRRuntime() {
        double[][] coord = new double[g.getVertexCount()][2];
        FRLayout<Integer, Integer> l = new FRLayout(g);
        l.setMaxIterations(5000);
        l.setSize(new Dimension(300, 300)); // sets the initial size of the space
// The BasicVisualizationServer<V,E> is parameterized by the edge types
        BasicVisualizationServer<Integer, Integer> vv = new BasicVisualizationServer<Integer, Integer>(l);
        vv.setPreferredSize(new Dimension(350, 350)); //Sets the viewing area size

        EllipseVertexShapeTransformer<Integer> t = new EllipseVertexShapeTransformer<Integer>();

        t.setSizeTransformer(new Transformer<Integer, Integer>() {
            public Integer transform(Integer i) {
                return 6;
            }
        });
        vv.getRenderContext().setVertexShapeTransformer(t);
//        JFrame frame = new JFrame("FR");
//        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//        frame.getContentPane().add(vv);
//        frame.pack();
//        frame.setVisible(true);
//

        while (!l.done()) {
//            String s = l.getStatus();
//            System.out.println(s);

            try {
                wait(100);
            } catch (InterruptedException ex) {
                Logger.getLogger(Visualization.class.getName()).log(Level.SEVERE, null, ex);
            }

        }
        for (Integer v : l.getGraph().getVertices()) {
            Point2D p = l.transform(v);
            coord[v][0] = p.getX();
            coord[v][1] = p.getY();
        }
        GraphCompression gk = new GraphCompression(g, coord);
        double costWithoutEmbedding = gk.codingCostNoEmbedding();
        System.out.println("costWithoutEmbedding: " + costWithoutEmbedding);
        double compressionCost = gk.mdlFunctionSimpleSigmoid();
        System.out.println("compression cost: " + compressionCost);
        return coord;

    }

    public synchronized double[][] getCoordinatesKKLRuntime() {
        double[][] coord = new double[g.getVertexCount()][2];
        KKLayout<Integer, Integer> l = new KKLayout(g);
        l.setMaxIterations(5000);
//        l.setDisconnectedDistanceMultiplier(0.3);
//        l.setLengthFactor(0.7);
        l.setSize(new Dimension(300, 300)); // sets the initial size of the space
// The BasicVisualizationServer<V,E> is parameterized by the edge types
        BasicVisualizationServer<Integer, Integer> vv = new BasicVisualizationServer<Integer, Integer>(l);
        vv.setPreferredSize(new Dimension(350, 350)); //Sets the viewing area size

        EllipseVertexShapeTransformer<Integer> t = new EllipseVertexShapeTransformer<Integer>();

        t.setSizeTransformer(new Transformer<Integer, Integer>() {
            public Integer transform(Integer i) {
                return 6;
            }
        });
        vv.getRenderContext().setVertexShapeTransformer(t);
//        JFrame frame = new JFrame("KKL");
//        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//        frame.getContentPane().add(vv);
//        frame.pack();
//        frame.setVisible(false);
//

        while (!l.done()) {
//            String s = l.getStatus();
//            System.out.println(s);

            try {
                wait(100);
            } catch (InterruptedException ex) {
                Logger.getLogger(Visualization.class.getName()).log(Level.SEVERE, null, ex);
            }

        }
        for (Integer v : l.getGraph().getVertices()) {
            Point2D p = l.transform(v);
            coord[v][0] = p.getX();
            coord[v][1] = p.getY();
        }
        GraphCompression gk = new GraphCompression(g, coord);
        double costWithoutEmbedding = gk.codingCostNoEmbedding();
        System.out.println("costWithoutEmbedding: " + costWithoutEmbedding);
        double compressionCost = gk.mdlFunctionSimpleSigmoid();
        System.out.println("compression cost: " + compressionCost);
        return coord;

    }

    public synchronized double[][] getCoordinatesKKL() {
        double[][] coord = new double[g.getVertexCount()][2];
        KKLayout<Integer, Integer> l = new KKLayout(g);
        l.setMaxIterations(5000);
//        l.setDisconnectedDistanceMultiplier(0.3);
//        l.setLengthFactor(0.7);
        l.setSize(new Dimension(300, 300)); // sets the initial size of the space
// The BasicVisualizationServer<V,E> is parameterized by the edge types
        BasicVisualizationServer<Integer, Integer> vv = new BasicVisualizationServer<Integer, Integer>(l);
        vv.setPreferredSize(new Dimension(350, 350)); //Sets the viewing area size

        EllipseVertexShapeTransformer<Integer> t = new EllipseVertexShapeTransformer<Integer>();

        t.setSizeTransformer(new Transformer<Integer, Integer>() {
            public Integer transform(Integer i) {
                return 6;
            }
        });
        vv.getRenderContext().setVertexShapeTransformer(t);
        JFrame frame = new JFrame("KKL");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().add(vv);
        frame.pack();
        frame.setVisible(true);
//

        while (!l.done()) {
//            String s = l.getStatus();
//            System.out.println(s);

            try {
                wait(100);
            } catch (InterruptedException ex) {
                Logger.getLogger(Visualization.class.getName()).log(Level.SEVERE, null, ex);
            }

        }
        for (Integer v : l.getGraph().getVertices()) {
            Point2D p = l.transform(v);
            coord[v][0] = p.getX();
            coord[v][1] = p.getY();
        }
        GraphCompression gk = new GraphCompression(g, coord);
        double costWithoutEmbedding = gk.codingCostNoEmbedding();
        System.out.println("costWithoutEmbedding: " + costWithoutEmbedding);
        double compressionCost = gk.mdlFunctionSimpleSigmoid();
        System.out.println("compression cost: " + compressionCost);
        return coord;

    }

    public synchronized double[][] getCoordinatesISOM() {
        double[][] coord = new double[g.getVertexCount()][2];
        ISOMLayout<Integer, Integer> l = new ISOMLayout(g);
        l.setSize(new Dimension(300, 300)); // sets the initial size of the space
// The BasicVisualizationServer<V,E> is parameterized by the edge types
        BasicVisualizationServer<Integer, Integer> vv = new BasicVisualizationServer<Integer, Integer>(l);
        vv.setPreferredSize(new Dimension(350, 350)); //Sets the viewing area size
        EllipseVertexShapeTransformer<Integer> t = new EllipseVertexShapeTransformer<Integer>();

        t.setSizeTransformer(new Transformer<Integer, Integer>() {
            public Integer transform(Integer i) {
                return 6;
            }
        });
        vv.getRenderContext().setVertexShapeTransformer(t);
        JFrame frame = new JFrame("ISOM");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().add(vv);
        frame.pack();
        frame.setVisible(true);
//

        while (!l.done()) {
//            String s = l.getStatus();
//            System.out.println(s);
            try {
                wait(100);
            } catch (InterruptedException ex) {
                Logger.getLogger(Visualization.class.getName()).log(Level.SEVERE, null, ex);
            }

        }
        for (Integer v : l.getGraph().getVertices()) {
            Point2D p = l.transform(v);
            coord[v][0] = p.getX();
            coord[v][1] = p.getY();
        }
        GraphCompression gk = new GraphCompression(g, coord);
        double costWithoutEmbedding = gk.codingCostNoEmbedding();
        System.out.println("costWithoutEmbedding: " + costWithoutEmbedding);
        double compressionCost = gk.mdlFunctionSimpleSigmoid();
        System.out.println("compression cost: " + compressionCost);

        return coord;

    }

    public synchronized double[][] getCoordinatesISOMRuntime() {
        double[][] coord = new double[g.getVertexCount()][2];
        ISOMLayout<Integer, Integer> l = new ISOMLayout(g);
        l.setSize(new Dimension(300, 300)); // sets the initial size of the space
// The BasicVisualizationServer<V,E> is parameterized by the edge types
        BasicVisualizationServer<Integer, Integer> vv = new BasicVisualizationServer<Integer, Integer>(l);
        vv.setPreferredSize(new Dimension(350, 350)); //Sets the viewing area size
        EllipseVertexShapeTransformer<Integer> t = new EllipseVertexShapeTransformer<Integer>();

        t.setSizeTransformer(new Transformer<Integer, Integer>() {
            public Integer transform(Integer i) {
                return 6;
            }
        });
        vv.getRenderContext().setVertexShapeTransformer(t);
//        JFrame frame = new JFrame("ISOM");
//        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//        frame.getContentPane().add(vv);
//        frame.pack();
//        frame.setVisible(true);
//

        while (!l.done()) {
//            String s = l.getStatus();
//            System.out.println(s);
            try {
                wait(100);
            } catch (InterruptedException ex) {
                Logger.getLogger(Visualization.class.getName()).log(Level.SEVERE, null, ex);
            }

        }
        for (Integer v : l.getGraph().getVertices()) {
            Point2D p = l.transform(v);
            coord[v][0] = p.getX();
            coord[v][1] = p.getY();
        }
        GraphCompression gk = new GraphCompression(g, coord);
        double costWithoutEmbedding = gk.codingCostNoEmbedding();
        System.out.println("costWithoutEmbedding: " + costWithoutEmbedding);
        double compressionCost = gk.mdlFunctionSimpleSigmoid();
        System.out.println("compression cost: " + compressionCost);

        return coord;

    }

    public synchronized double[][] getCoordinatesCircle() {
        double[][] coord = new double[g.getVertexCount()][2];
        CircleLayout<Integer, Integer> l = new CircleLayout(g);
        l.setSize(new Dimension(1000, 1000)); // sets the initial size of the space

// The BasicVisualizationServer<V,E> is parameterized by the edge types
        BasicVisualizationServer<Integer, Integer> vv = new BasicVisualizationServer<Integer, Integer>(l);
        vv.setBackground(Color.WHITE);
        vv.setPreferredSize(new Dimension(950, 950)); //Sets the viewing area size
        EllipseVertexShapeTransformer<Integer> t = new EllipseVertexShapeTransformer<Integer>();

        t.setSizeTransformer(new Transformer<Integer, Integer>() {
            public Integer transform(Integer i) {
                return 10;
            }
        });
        vv.getRenderContext().setVertexShapeTransformer(t);
        JFrame frame = new JFrame("Circle");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().add(vv);
        frame.pack();
        frame.setVisible(true);
//

//        while (!l.done()) {
////            String s = l.getStatus();
////            System.out.println(s);
//            try {
//                wait(100);
//            } catch (InterruptedException ex) {
//                Logger.getLogger(Isomap.class.getName()).log(Level.SEVERE, null, ex);
//            }
//
//        }
        for (Integer v : l.getGraph().getVertices()) {
            Point2D p = l.transform(v);
            coord[v][0] = p.getX();
            coord[v][1] = p.getY();
        }
        return coord;

    }

    public synchronized double[][] getCoordinatesSpring() {
        double[][] coord = new double[g.getVertexCount()][2];
        SpringLayout<Integer, Integer> l = new SpringLayout(g);
        l.setSize(new Dimension(300, 300)); // sets the initial size of the space
// The BasicVisualizationServer<V,E> is parameterized by the edge types
        BasicVisualizationServer<Integer, Integer> vv = new BasicVisualizationServer<Integer, Integer>(l);
        vv.setPreferredSize(new Dimension(350, 350)); //Sets the viewing area size
        EllipseVertexShapeTransformer<Integer> t = new EllipseVertexShapeTransformer<Integer>();

        t.setSizeTransformer(new Transformer<Integer, Integer>() {
            public Integer transform(Integer i) {
                return 6;
            }
        });
        vv.getRenderContext().setVertexShapeTransformer(t);
        JFrame frame = new JFrame("Spring");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().add(vv);
        frame.pack();
        frame.setVisible(true);
//

//        while (!l.done()) {
////            String s = l.getStatus();
////            System.out.println(s);
        try {
            wait(50000);
        } catch (InterruptedException ex) {
            Logger.getLogger(Visualization.class.getName()).log(Level.SEVERE, null, ex);
        }

        // }
        for (Integer v : l.getGraph().getVertices()) {
            Point2D p = l.transform(v);
            coord[v][0] = p.getX();
            coord[v][1] = p.getY();
        }
        GraphCompression gk = new GraphCompression(g, coord);
        double costWithoutEmbedding = gk.codingCostNoEmbedding();
        System.out.println("costWithoutEmbedding: " + costWithoutEmbedding);
        double compressionCost = gk.mdlFunctionSimpleSigmoid();
        System.out.println("compression cost: " + compressionCost);
        return coord;

    }












    //    //8.3.16: Laplacian Eigenmaps minimizing squared Error (Ink)
//    public double[][] getCoordinatesLaplacianOnly(int dim) {
//        double[][] res = new double[g.getVertexCount()][dim];
//        double[][] d = new double[g.getVertexCount()][g.getVertexCount()];
//        for (int i = 0; i < g.getVertexCount(); i++) {
//            d[i][i] = g.degree(i);
//        }
//        double[][] l = new double[g.getVertexCount()][g.getVertexCount()];
//        for (int i = 0; i < g.getVertexCount(); i++) {
//            for (int j = 0; j < g.getVertexCount(); j++) {
//                if (i == j) {
//                    l[i][j] = d[i][j] - 1.0;
//                } else if (g.isNeighbor(i, j)) {
//                    l[i][j] = l[j][i] = -1.0;
//                }
//            }
//
//        }
//        DoubleMatrix D = new DoubleMatrix(d);
//        DoubleMatrix L = new DoubleMatrix(l);
//        DoubleMatrix[] V = symmetricGeneralizedEigenvectors(L, D, 1, 2);
//        for (int i = 0; i < res.length; i++) {
//            for (int j = 0; j < dim; j++) {
//                res[i][j] = V[0].get(i, j);
//            }
//        }
//        return res;
//    }




    public void displayCoordNew(double[][] coord, String title) {
        double[][] coordv = new double[coord.length][coord[0].length];
        for (int i = 0; i < coord.length; i++) {
            for (int j = 0; j < coord[0].length; j++) {
                coordv[i][j] = coord[i][j];
            }
        }
        int drawingSize = 800; //original values
        int displaySize = 1000;

//        int drawingSize = 1600;
//        int displaySize = 2000;
        int offSet = (displaySize - drawingSize) / 2;
        double max_x = -java.lang.Double.MAX_VALUE;
        double max_y = -java.lang.Double.MAX_VALUE;
        double min_x = -max_x;
        double min_y = -max_y;

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
                coordv[i][0] = ((coord[i][0] - min_x) / (max_x - min_x)) * drawingSize + offSet;
                coordv[i][1] = ((coord[i][1] - min_y) / (max_x - min_x)) * drawingSize + offSet;
            }
        } else {
            for (int i = 0; i < numObj; i++) {
                coordv[i][0] = ((coord[i][0] - min_x) / (max_y - min_y)) * drawingSize + offSet;
                coordv[i][1] = ((coord[i][1] - min_y) / (max_y - min_y)) * drawingSize + offSet;
            }
        }

        StaticLayout<Integer, Integer> l = new StaticLayout(g);
        l.setSize(new Dimension(displaySize, displaySize)); // sets the initial size of the space
        //BasicVisualizationServer<Integer, Integer> vv = new BasicVisualizationServer<Integer, Integer>(l);
        VisualizationViewer<Integer, Integer> vv = new VisualizationViewer<Integer, Integer>(l);
        vv.setPreferredSize(new Dimension(displaySize, displaySize)); //Sets the viewing area size
        // vv.getRenderContext().setVertexFillPaintTransformer(new MyVertexFillPaintFunction(clid));
        vv.setBackground(Color.WHITE);
        EllipseVertexShapeTransformer<Integer> t = new EllipseVertexShapeTransformer<Integer>();
        t.setSizeTransformer(new Transformer<Integer, Integer>() {
            public Integer transform(Integer i) {
                return 5;
                //return 7;
            }
        });

        // double diffCenter_x =
        HashMap<Integer, Point2D> map = new HashMap<Integer, Point2D>();
        for (int i = 0; i < numObj; i++) {
            Point2D.Double p = new Point2D.Double();
            p.setLocation(coordv[i][0], coordv[i][1]);
            map.put(i, p);
        }

        final Stroke edgeStroke = new BasicStroke(5.0f);
        Transformer<Integer, Stroke> edgeStrokeTransformer = new Transformer<Integer, Stroke>() {
            public Stroke transform(Integer s) {
                return edgeStroke;
            }
        };

        Transformer vertexLocations = TransformerUtils.mapTransformer(map);
        vv.getRenderContext().setVertexShapeTransformer(t);
        // vv.getRenderContext().setEdgeStrokeTransformer(edgeStrokeTransformer);
        //vv.getRenderContext().setVertexLabelTransformer(new ToStringLabeller());
        //vv.getRenderer().getVertexLabelRenderer().setPosition(Position.CNTR);
        vv.addGraphMouseListener(new TestGraphMouseListener<Integer>());
        AbstractEdgeShapeTransformer<Integer, Integer> ast = new EdgeShape.Line<Integer, Integer>();

        vv.getRenderContext().setEdgeShapeTransformer(ast);

        l.setInitializer(vertexLocations);

        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().add(vv);
        frame.pack();
        frame.setVisible(true);

    }

    public void displayDegreeDistribution(double[][] coord, String title) {
        int[] clid = new int[g.getVertexCount()];
        for (int i = 0; i < g.getVertexCount(); i++) {
            clid[i] = g.degree(i);
            System.out.println("i " + i + " " + g.degree(i));
        }
        DataUtils du = new DataUtils();

        double[][] coordv = new double[coord.length][coord[0].length];
        for (int i = 0; i < coord.length; i++) {
            for (int j = 0; j < coord[0].length; j++) {
                coordv[i][j] = coord[i][j];
            }
        }
        int drawingSize = 800;
        int displaySize = 1000;
        int offSet = (displaySize - drawingSize) / 2;
        double max_x = -java.lang.Double.MAX_VALUE;
        double max_y = -java.lang.Double.MAX_VALUE;
        double min_x = -max_x;
        double min_y = -max_y;
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
                coordv[i][0] = ((coord[i][0] - min_x) / (max_x - min_x)) * drawingSize + offSet;
                coordv[i][1] = ((coord[i][1] - min_y) / (max_x - min_x)) * drawingSize + offSet;
            }
        } else {
            for (int i = 0; i < numObj; i++) {
                coordv[i][0] = ((coord[i][0] - min_x) / (max_y - min_y)) * drawingSize + offSet;
                coordv[i][1] = ((coord[i][1] - min_y) / (max_y - min_y)) * drawingSize + offSet;
            }
        }

        StaticLayout<Integer, Integer> l = new StaticLayout(g);
        l.setSize(new Dimension(displaySize, displaySize)); // sets the initial size of the space
        //BasicVisualizationServer<Integer, Integer> vv = new BasicVisualizationServer<Integer, Integer>(l);
        VisualizationViewer<Integer, Integer> vv = new VisualizationViewer<Integer, Integer>(l);
        vv.setPreferredSize(new Dimension(displaySize, displaySize)); //Sets the viewing area size
        vv.getRenderContext().setVertexFillPaintTransformer(new MyVertexFillPaintFunctionDegree(clid));
        vv.setBackground(Color.WHITE);
        EllipseVertexShapeTransformer<Integer> t = new EllipseVertexShapeTransformer<Integer>();
        t.setSizeTransformer(new Transformer<Integer, Integer>() {
            public Integer transform(Integer i) {
                return 20;
            }
        });

        // double diffCenter_x =
        HashMap<Integer, Point2D> map = new HashMap<Integer, Point2D>();
        for (int i = 0; i < numObj; i++) {
            Point2D.Double p = new Point2D.Double();
            p.setLocation(coordv[i][0], coordv[i][1]);
            map.put(i, p);
        }

        Transformer vertexLocations = TransformerUtils.mapTransformer(map);
        vv.getRenderContext().setVertexShapeTransformer(t);
        //vv.getRenderContext().setVertexLabelTransformer(new ToStringLabeller());
        //vv.getRenderer().getVertexLabelRenderer().setPosition(Position.CNTR);
        vv.addGraphMouseListener(new TestGraphMouseListener<Integer>());
        AbstractEdgeShapeTransformer<Integer, Integer> ast = new EdgeShape.Line<Integer, Integer>();

        vv.getRenderContext().setEdgeShapeTransformer(ast);

        l.setInitializer(vertexLocations);

        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().add(vv);
        frame.pack();
        frame.setVisible(true);

    }


    public void displayCoordNew(double[][] coord, String title, int[] clid) {
        double[][] coordv = new double[coord.length][coord[0].length];
        for (int i = 0; i < coord.length; i++) {
            for (int j = 0; j < coord[0].length; j++) {
                coordv[i][j] = coord[i][j];
            }
        }
        int drawingSize = 800;
        int displaySize = 1000;
        int offSet = (displaySize - drawingSize) / 2;
        double max_x = -java.lang.Double.MAX_VALUE;
        double max_y = -java.lang.Double.MAX_VALUE;
        double min_x = -max_x;
        double min_y = -max_y;
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
                coordv[i][0] = ((coord[i][0] - min_x) / (max_x - min_x)) * drawingSize + offSet;
                coordv[i][1] = ((coord[i][1] - min_y) / (max_x - min_x)) * drawingSize + offSet;
            }
        } else {
            for (int i = 0; i < numObj; i++) {
                coordv[i][0] = ((coord[i][0] - min_x) / (max_y - min_y)) * drawingSize + offSet;
                coordv[i][1] = ((coord[i][1] - min_y) / (max_y - min_y)) * drawingSize + offSet;
            }
        }

        StaticLayout<Integer, Integer> l = new StaticLayout(g);
        l.setSize(new Dimension(displaySize, displaySize)); // sets the initial size of the space
        //BasicVisualizationServer<Integer, Integer> vv = new BasicVisualizationServer<Integer, Integer>(l);
        VisualizationViewer<Integer, Integer> vv = new VisualizationViewer<Integer, Integer>(l);
        vv.setPreferredSize(new Dimension(displaySize, displaySize)); //Sets the viewing area size
        vv.getRenderContext().setVertexFillPaintTransformer(new MyVertexFillPaintFunction(clid));
        vv.setBackground(Color.WHITE);
//        Transformer<Integer, Paint> edgePaint = new Transformer<Integer, Paint>() {
//            public Paint transform(Integer i) {
//                return Color.WHITE;
//            }
//        };
        final Stroke edgeStroke = new BasicStroke(2.0f);
        Transformer<Integer, Stroke> edgeStrokeTransformer = new Transformer<Integer, Stroke>() {
            public Stroke transform(Integer s) {
                return edgeStroke;
            }
        };
//        Transformer<MyEdge, Stroke> edgeStrokeTransformerWeighted = new Transformer<MyEdge, Stroke>() {
//            public Stroke transform(MyEdge s) {
//                return edgeStroke;
//            }
//        };

        EllipseVertexShapeTransformer<Integer> t = new EllipseVertexShapeTransformer<Integer>();
        t.setSizeTransformer(new Transformer<Integer, Integer>() {
            public Integer transform(Integer i) {
                //return 20; //small graphs
                return 7;
            }
        });

        // double diffCenter_x =
        HashMap<Integer, Point2D> map = new HashMap<Integer, Point2D>();
        for (int i = 0; i < numObj; i++) {
            Point2D.Double p = new Point2D.Double();
            p.setLocation(coordv[i][0], coordv[i][1]);
            map.put(i, p);
        }

        Transformer vertexLocations = TransformerUtils.mapTransformer(map);
        vv.getRenderContext().setVertexShapeTransformer(t);
        //vv.getRenderContext().setEdgeDrawPaintTransformer(edgePaint);
        vv.getRenderContext().setEdgeStrokeTransformer(edgeStrokeTransformer);
        //vv.getRenderContext().setVertexLabelTransformer(new ToStringLabeller());
        //vv.getRenderer().getVertexLabelRenderer().setPosition(Position.CNTR);
        vv.addGraphMouseListener(new TestGraphMouseListener<Integer>());
        AbstractEdgeShapeTransformer<Integer, Integer> ast = new EdgeShape.Line<Integer, Integer>();

        vv.getRenderContext().setEdgeShapeTransformer(ast);

        l.setInitializer(vertexLocations);

        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().add(vv);
        frame.pack();
        frame.setVisible(true);

    }

    public class MyVertexFillPaintFunctionDegree implements Transformer<Integer, Paint> {

        int[] ids;
        int maxDegree;

        public void setIds(int[] ids) {
            this.ids = ids;
            maxDegree = 0;
            for (int i = 0; i < ids.length; i++) {
                if (ids[i] > maxDegree) {
                    maxDegree = ids[i];
                }
            }
        }

        public MyVertexFillPaintFunctionDegree(int[] ids) {
            this.ids = ids;
            maxDegree = 0;
            for (int i = 0; i < ids.length; i++) {
                if (ids[i] > maxDegree) {
                    maxDegree = ids[i];
                }
            }
        }

        public Paint transform(Integer vv) {
            double brightness = (double) ids[vv.intValue()] / (double) maxDegree;
            int bb = (int) Math.rint(brightness * 255);

            return new Color(bb, bb, bb);
        }
        //int vv = v.intValue();
//            if (ids[vv] == 0) {
//                return Color.BLUE;
//            }
//            if (ids[vv] == 1) {
//                return Color.RED;
//            }
//            if (ids[vv] == 2) {
//                return Color.LIGHT_GRAY;
//
//            }
//        }
    }


    public class MyVertexFillPaintFunctionGeo implements Transformer<Integer, Paint> {

        double[][] coord;

        public void setCoord(double[][] coord) {
            this.coord = coord;
        }

        public MyVertexFillPaintFunctionGeo(double[][] coord) {
            this.coord = coord;
        }

        public Paint transform(Integer vv) {
            double xxx=((coord[vv][0]+68.0652)/5.3781);
            double yyy=((coord[vv][1]-43.499)/(49.001-43.499));
            Color cc = new Color((int)Math.min(255.,yyy*350.), (int)Math.max(0.,(255.-350.*yyy)), (int)(0.*xxx)); //lightred
//                return cc;
//                Color cc = new Color(0, 100, 0); //Sigmod, lightgreen
//                return cc;
            return cc; //football
            //return Color.LIGHT_GRAY;
            //return Color.BLUE;
        }



    }


    public class MyVertexFillPaintFunction implements Transformer<Integer, Paint> {

        int[] ids;

        public void setIds(int[] ids) {
            this.ids = ids;
        }

        public MyVertexFillPaintFunction(int[] ids) {
            this.ids = ids;
        }

        public Paint transform(Integer vv) {
            if (ids[vv.intValue()] == 0) { //KDD
//                Color cc = new Color(255, 200, 200); //lightred
//                return cc;
//                Color cc = new Color(0, 100, 0); //Sigmod, lightgreen
//                return cc;
                return Color.red; //football
                //return Color.LIGHT_GRAY;
                //return Color.BLUE;
            }

            if (ids[vv.intValue()] == 1) { //ICDM
//                Color cc = new Color(255, 0, 0); //red
//                return cc;
                return Color.BLUE; //football
                //return Color.GREEN;
            }
            if (ids[vv.intValue()] == 2) {
//                Color cc = new Color(130, 0, 0); //darkred, SDM
//                return cc;
                //return Color.GREEN; //football
                return Color.YELLOW;
            }
            if (ids[vv.intValue()] == 3) {
//                Color cc = new Color(0, 100, 0); //Sigmod, lightgreen
//                return cc;
//                airflights
                return Color.RED;

                // return Color.ORANGE; //football
            }
            if (ids[vv.intValue()] == 4) {
//                Color cc = new Color(0, 255, 0); //VLDB, green
//                return cc;
                return Color.MAGENTA; //football
                //return Color.red;
            }
            if (ids[vv.intValue()] == 5) {
//                Color cc = new Color(0, 130, 0); //ICDE, darkgreen
//                return cc;

//airflights
                return Color.WHITE;

                //return Color.CYAN; //football
            }
            if (ids[vv.intValue()] == 6) {
//                Color cc = new Color(150, 150, 255); //STOC, lightblue
//                return cc;

//airflights
                return Color.CYAN;
//football
                //return Color.PINK;
            }
            if (ids[vv.intValue()] == 7) {
//                Color cc = new Color(0, 0, 255); //SODA, blue
//                return cc;
                return Color.GREEN;
            }

            if (ids[vv.intValue()] == 8) {
                Color cc = new Color(132, 60, 12);//FOCS, darkblue
                return cc;
                //return Color.DARK_GRAY;
            }

            if (ids[vv.intValue()] == 9) {
                return Color.LIGHT_GRAY;
            }
            if (ids[vv.intValue()] == 10) {
                return Color.YELLOW;
            }
            if (ids[vv.intValue()] == 11) {
                return Color.BLACK;
            }

            return Color.yellow;
        }
        //int vv = v.intValue();
//            if (ids[vv] == 0) {
//                return Color.BLUE;
//            }
//            if (ids[vv] == 1) {
//                return Color.RED;
//            }
//            if (ids[vv] == 2) {
//                return Color.LIGHT_GRAY;
//
//            }
//        }
    }

    public class MyVertexSizePaintFunction implements Transformer<Integer, Shape> {

        public MyVertexSizePaintFunction() {
        }

        public Shape transform(Integer vv) {
            Ellipse2D circle = new Ellipse2D.Double(0, 0, 3, 3);
            // in this case, the vertex is twice as large
            return AffineTransform.getScaleInstance(2, 2).createTransformedShape(circle);

        }
    }





    /**
     * A nested class to demo the GraphMouseListener finding the right vertices
     * after zoom/pan
     */
    static class TestGraphMouseListener<V> implements GraphMouseListener<V> {

        public void graphClicked(V v, MouseEvent me) {
            System.err.println("Vertex " + v + " was clicked at (" + me.getX() + "," + me.getY() + ")");
        }

        public void graphPressed(V v, MouseEvent me) {
            System.err.println("Vertex " + v + " was pressed at (" + me.getX() + "," + me.getY() + ")");
        }

        public void graphReleased(V v, MouseEvent me) {
            System.err.println("Vertex " + v + " was released at (" + me.getX() + "," + me.getY() + ")");
        }
    }

    public void displayCoord(double[][] coord, double[] costs, String title) {
        this.costs = new double[costs.length];
        for (int i = 0; i < costs.length; i++) {
            this.costs[i] = costs[i];
        }
        double[][] coordv = new double[coord.length][coord[0].length];
        for (int i = 0; i < coord.length; i++) {
            for (int j = 0; j < coord[0].length; j++) {
                coordv[i][j] = coord[i][j];
            }
        }
        int displaySize = 1000;
        double max_x = -java.lang.Double.MAX_VALUE;
        double max_y = -java.lang.Double.MAX_VALUE;
        double min_x = -max_x;
        double min_y = -max_y;
        for (int i = 0; i < numObj; i++) {

            if (coord[i][0] > max_x) {
                max_x = coord[i][0];
            }
            if (coord[i][1] < min_x) {
                min_x = coord[i][1];
            }
            if (coord[i][1] > max_y) {
                max_y = coord[i][1];
            }
            if (coord[i][1] < min_y) {
                min_y = coord[i][1];
            }
        }
        for (int i = 0; i < numObj; i++) {
            coordv[i][0] = ((coord[i][0] - min_x) / (max_x - min_x)) * (0.8 * displaySize);
            coordv[i][0] += (0.1 * displaySize);
            coordv[i][1] = ((coord[i][1] - min_y) / (max_y - min_y)) * (0.8 * displaySize);
            coordv[i][1] += (0.1 * displaySize);
        }

        StaticLayout<Integer, Integer> l = new StaticLayout(g);
        l.setSize(new Dimension(displaySize, displaySize)); // sets the initial size of the space
        BasicVisualizationServer<Integer, Integer> vv = new BasicVisualizationServer<Integer, Integer>(l);
        vv.setPreferredSize(new Dimension(displaySize, displaySize)); //Sets the viewing area size
        EllipseVertexShapeTransformer<Integer> t = new EllipseVertexShapeTransformer<Integer>();
        t.setSizeTransformer(new Transformer<Integer, Integer>() {
            public Integer transform(Integer i) {
                return 20;
            }
        });

        // double diffCenter_x =
        HashMap<Integer, Point2D> map = new HashMap<Integer, Point2D>();
        for (int i = 0; i < numObj; i++) {
            Point2D.Double p = new Point2D.Double();
            p.setLocation(coordv[i][0], coordv[i][1]);
            map.put(i, p);
        }

        Transformer vertexLocations = TransformerUtils.mapTransformer(map);
        vv.getRenderContext().setVertexShapeTransformer(t);
        vv.getRenderContext().setVertexLabelTransformer(new ToStringLabeller());
        vv.getRenderer().getVertexLabelRenderer().setPosition(Position.CNTR);
        vv.getRenderContext().setEdgeDrawPaintTransformer(edgePaint);
        l.setInitializer(vertexLocations);

        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().add(vv);
        frame.pack();
        frame.setVisible(true);

    }
    Transformer<Integer, Paint> edgePaint = new Transformer<Integer, Paint>() {
        public Paint transform(Integer s) {
            if (costs[s] > 1) {
                return Color.RED;
            } else {
                return Color.BLACK;
            }
        }
    };

    public void displayCoord(double[][] coord, String title) {
        double[][] coordv = new double[coord.length][coord[0].length];
        for (int i = 0; i < coord.length; i++) {
            for (int j = 0; j < coord[0].length; j++) {
                coordv[i][j] = coord[i][j];
            }
        }
        int displaySize = 1000;
        double max_x = -java.lang.Double.MAX_VALUE;
        double max_y = -java.lang.Double.MAX_VALUE;
        double min_x = -max_x;
        double min_y = -max_y;
        for (int i = 0; i < numObj; i++) {

            if (coord[i][0] > max_x) {
                max_x = coord[i][0];
            }
            if (coord[i][1] < min_x) {
                min_x = coord[i][1];
            }
            if (coord[i][1] > max_y) {
                max_y = coord[i][1];
            }
            if (coord[i][1] < min_y) {
                min_y = coord[i][1];
            }
        }
        for (int i = 0; i < numObj; i++) {
            coordv[i][0] = ((coord[i][0] - min_x) / (max_x - min_x)) * (0.8 * displaySize);
            coordv[i][0] += (0.1 * displaySize);
            coordv[i][1] = ((coord[i][1] - min_y) / (max_y - min_y)) * (0.8 * displaySize);
            coordv[i][1] += (0.1 * displaySize);
        }

        StaticLayout<Integer, Integer> l = new StaticLayout(g);
        l.setSize(new Dimension(displaySize, displaySize)); // sets the initial size of the space
        BasicVisualizationServer<Integer, Integer> vv = new BasicVisualizationServer<Integer, Integer>(l);
        vv.setPreferredSize(new Dimension(displaySize, displaySize)); //Sets the viewing area size
        EllipseVertexShapeTransformer<Integer> t = new EllipseVertexShapeTransformer<Integer>();
        t.setSizeTransformer(new Transformer<Integer, Integer>() {
            public Integer transform(Integer i) {
                return 5;
            }
        });

        // double diffCenter_x =
        HashMap<Integer, Point2D> map = new HashMap<Integer, Point2D>();
        for (int i = 0; i < numObj; i++) {
            Point2D.Double p = new Point2D.Double();
            p.setLocation(coordv[i][0], coordv[i][1]);
            map.put(i, p);
        }

        Transformer vertexLocations = TransformerUtils.mapTransformer(map);
        vv.getRenderContext().setVertexShapeTransformer(t);
        vv.getRenderContext().setVertexLabelTransformer(new ToStringLabeller());
        vv.getRenderer().getVertexLabelRenderer().setPosition(Position.CNTR);
        // vv.getRenderContext().setEdgeDrawPaintTransformer(edgePaint);
        l.setInitializer(vertexLocations);

        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().add(vv);
        frame.pack();
        frame.setVisible(true);

    }

    public void getCoordinatesOwn(int dim, boolean mdlFunct) {
        int numBins = 30;

        double[][] coord = new double[numObj][dim];
        double[][] bestCoord = new double[numObj][dim];
//        double[][] bestCoord = new double[numObj][dim];
//        double[][] worstCoord = new double[numO]
        Random r = new Random(10);
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < dim; j++) {
                coord[i][j] = r.nextDouble();
            }
        }
        int numVertices = g.getVertexCount();
        int numEdges = g.getEdgeCount();
        int numPossEdges = (numVertices * (numVertices - 1)) / 2;
        double pLink_g = (double) (numEdges * 2) / (double) numVertices;
        double pNoLink_g = 1.0 - pLink_g;
        GraphCompression gc = new GraphCompression(g, coord, numBins);
        System.out.println("codingCostNoEmbedding: " + gc.codingCostEmbedding());
        double baseMDL = 0.0;
        if (mdlFunct) {
            baseMDL = gc.mdlFunction();
        } else {
            baseMDL = gc.mdlHisto();
        }
        System.out.println("baseMDL: " + baseMDL);
        int iter = 0;
        double bestMdl = java.lang.Double.MAX_VALUE;
        while (iter < 100) {
            for (int i = 0; i < numObj; i++) {
                // System.out.println(i);
                double[] meanLink = new double[dim];
                double[] meanNotLink = new double[dim];
                int counter_links = 0;
                int counter_noLinks = 0;
                for (int k = 0; k < numObj; k++) {
                    if (g.isNeighbor(i, k)) {
                        counter_links++;
                        for (int j = 0; j < dim; j++) {
                            meanLink[j] += coord[k][j];
                        }
                    } else {
                        counter_noLinks++;
                        for (int j = 0; j < dim; j++) {
                            meanNotLink[j] += coord[k][j];
                        }
                    }
                }
                if (counter_noLinks > 0 && counter_links > 0) {
                    for (int j = 0; j < dim; j++) {
                        meanLink[j] /= (double) counter_links;
                        meanNotLink[j] /= (double) counter_noLinks;
                    }
                    double[] maxDistantPoint = new double[dim];
                    for (int k = 0; k < numObj; k++) {
                        for (int j = 0; j < dim; j++) {
                            if ((Math.abs(coord[k][j] - meanNotLink[j])) > maxDistantPoint[j]) {
                                maxDistantPoint[j] = coord[k][j];
                            }
                        }
                    }
                    //option 1: move to center of links
                    if (iter == 12) {
                        System.out.println("m");
                    }
                    if (iter % 2 == 0) {
                        for (int j = 0; j < dim; j++) {
                            coord[i][j] = (coord[i][j] + meanLink[j]) / 2;
                        }
                    } else {

                        //option 2 : move away from not links
                        for (int j = 0; j < dim; j++) {
                            coord[i][j] = (coord[i][j] + maxDistantPoint[j]) / 2;
                        }
                    }

//                        //weighted sum
//                         for (int j = 0; j < dim; j++) {
//                            coord[i][j] = (((double) (counter_noLinks) * maxDistantPoint[j]) + ((double) (counter_links) * coord[i][j])) / (counter_links + counter_noLinks);
//                        }
                }
//                    if (iter % 2 == 0) {
//                        //System.out.println(iter);
//                        System.arraycopy(meanLink, 0, coord[i], 0, dim);
//                    } else {
//                        double[] maxDistantPoint = new double[dim];
//                        for (int k = 0; k < numObj; k++) {
//                            for (int j = 0; j < dim; j++) {
//                                if ((Math.abs(coord[k][j] - meanNotLink[j])) > maxDistantPoint[j]) {
//                                    maxDistantPoint[j] = coord[i][j];
//                                }
//                            }
//                        }
//                        for (int j = 0; j < dim; j++) {
//                            coord[i][j] = ((double) (counter_noLinks) * maxDistantPoint[j]) + ((double) (counter_links) * coord[i][j]) / (counter_links + counter_noLinks);
//                        }
//                        System.arraycopy(maxDistantPoint, 0, coord[i], 0, dim);
//                    }

                if (counter_links == 0 && counter_noLinks > 0) {
                    double[] maxDistantPoint = new double[dim];
                    for (int k = 0; k < numObj; k++) {
                        for (int j = 0; j < dim; j++) {
                            if ((Math.abs(coord[k][j] - meanNotLink[j])) > maxDistantPoint[j]) {
                                maxDistantPoint[j] = coord[k][j];
                            }
                        }
                    }
                    System.arraycopy(maxDistantPoint, 0, coord[i], 0, dim);

                }
                if (counter_links > 0 && counter_noLinks == 0) {
                    for (int j = 0; j < dim; j++) {
                        meanLink[j] /= (double) counter_links;
                    }
                    System.arraycopy(meanLink, 0, coord[i], 0, dim);
                }
            }

            gc = new GraphCompression(g, coord, numBins);
//            if (iter == 89) {
//                System.out.println("m");
//            }
//            if  (iter == 95){
//                System.out.println("m");
//            }
            double newMdl = 0.0;
            if (mdlFunct) {
                newMdl = gc.mdlFunction();
            } else {
                newMdl = gc.mdlHisto();
            }
            //double newMdl = gc.mdlFunction();
            System.out.println("iter: " + iter + " " + newMdl);
            if (newMdl < bestMdl) {
                bestMdl = newMdl;
                for (int l = 0; l < numObj; l++) {
                    System.arraycopy(coord[l], 0, bestCoord[l], 0, dim);
                }
                //re-scale bestCoordinates between 0 and 1
                double[] min = new double[dim];
                double[] max = new double[dim];
                for (int j = 0; j < dim; j++) {
                    min[j] = java.lang.Double.MAX_VALUE;
                    max[j] = -min[j];
                }
                for (int l = 0; l < numObj; l++) {
                    for (int j = 0; j < dim; j++) {
                        if (bestCoord[l][j] < min[j]) {
                            min[j] = bestCoord[l][j];
                        }
                        if (bestCoord[l][j] > max[j]) {
                            max[j] = bestCoord[l][j];
                        }
                    }
                }
                for (int l = 0; l < numObj; l++) {
                    for (int j = 0; j < dim; j++) {
                        bestCoord[l][j] = (bestCoord[l][j] - min[j]) / (max[j] - min[j]);

                    }
                }

            }
            iter++;

        }

//            if (iter == 89) {
//                dataCheck(coord, "iter_89");
//            }
//            if (iter == 95) {
//                dataCheck(coord, "iter_95");
//            }
//            if(iter == 43)
//                dataCheck(coord, "iter_43");
//            iter++;
//            if(iter == 84)
//                dataCheck(coord, "iter_84");
        System.out.println("bestMdl: " + bestMdl);
        dataCheck(bestCoord, "best");
        int test_260 = g.getNeighborCount(260);
        int test_171 = g.getNeighborCount(171);
        displayCoord(bestCoord, "own");
        //display best coordinates

    }

    //
//    public void getCoordinatesMDSWithoutSp(){
//         for (int i = 0; i < numObj; i++) {
//            for (int j = 0; j < numObj; j++) {
//                if (i > j) {
//                    if(g.)
//
//            }
//        }
//        }
//    }
    private double euclideanDistance(double[] o, double[] p) {
        double dist = 0.0;
        for (int i = 0; i
                < o.length; i++) {
            dist += (o[i] - p[i]) * (o[i] - p[i]);
        }
        return Math.sqrt(dist);

    }


    private void dataCheck(double[][] d, String filename) {
        MLDouble q = new MLDouble("test", d);
        ArrayList ll = new ArrayList();
        ll.add(q);
        MatFileWriter mw = new MatFileWriter();
        try {
            String name = filename + ".mat";
            mw.write(name, ll);
        } catch (IOException ex) {
        }
    }
}
