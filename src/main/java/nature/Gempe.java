/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import edu.uci.ics.jung.graph.Graph;
import utils.IO;

/**
 *
 * @author plantc59cs
 */
public class Gempe {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        if (args.length != 5) {
            System.out.println("Please run GEMPE with arguments: input file name, input variable name, number of threads, display (true/false), save intermediate results (true/false)");
        } else {
            String inputFileName = args[0];
            String inputVarName = args[1];
            int d = 2;
            int numT = Integer.parseInt(args[2]);
            boolean display = Boolean.valueOf(args[3]);
            boolean intermediate = Boolean.valueOf(args[4]);

            IO ea = new IO();
            Graph g = ea.matlabToGraph(inputFileName, inputVarName);
            Parallel gp = new Parallel(g, numT, d);
            gp.saveIntermediateResults = intermediate;
            gp.display = display;
            // gp.twoDGridStart();
            gp.twoD();
            System.out.println("finished");

        }
    }
}
