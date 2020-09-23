/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

/**
 *
 * @author plantc59cs
 */
public class DWUpdateResult {

    double[] d;
    double[] w;
    int[][] ij;

    public DWUpdateResult(int size) {
        d = new double[size];
        w = new double[size];
        ij = new int[size][2];
        for(int i = 0; i < size; i++){
            this.ij[i][0] = -1; //to check if this has been already written
        }
    }

}
