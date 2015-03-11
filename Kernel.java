/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package gaborfilter;

/**
 *
 * @author Kyle
 */
public class Kernel {
    
    static double[][] kernel;
    static int kval;
    
    public Kernel(int k) {
        
        kval = k;
        kernel = new double[kval][kval];
        return;
    }
}
