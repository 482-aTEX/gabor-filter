/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package gaborfilter;

import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.imageio.ImageIO;

/**
 *
 * @author Kyle
 */
public class GaborFilter {

    static int height, width, kval;
    static double[][] kernel;
    
    static void gaborKernel(double lambda, double theta, double psi, double sigma, double gamma) {
        
        kernel = new double[kval][kval];
        
 	theta = (theta * Math.PI) / (double)180;
	psi = (psi * Math.PI) / (double)180;
        
	double gauss, sinusoid, xprime, yprime;
        //for each determined pixel of the kernel, we calculate its value according to equations seen:
        //http://en.wikipedia.org/wiki/Gabor_filter
        int x, y;
	for(y = -(kval/2); y < (kval/2)+1; ++y) {
		for(x = -(kval/2); x < (kval/2)+1; ++x) {
		
			xprime = (x * Math.cos(theta)) + (y * Math.sin(theta));
			yprime = -(x * Math.sin(theta)) + (y * Math.cos(theta));
                            
			gauss = Math.exp(-( (xprime*xprime + lambda*lambda*yprime*yprime)/(double)(2.0*sigma*sigma)));		
                        sinusoid = Math.cos((2 * Math.PI * (xprime / lambda)) + psi);
			
			kernel[y+(kval/2)][x+(kval/2)] = gauss * sinusoid;
		}
	}
        return;
    }
    
    static void normalizeKernel() {
        
        double sum = 0;
        int i, j;
        for(i = 0; i < kval; ++i) {
            for(j = 0; j < kval; ++j) {
                
                sum = sum + kernel[i][j];
            }
        }    
        
        for(i = 0; i < kval; ++i) {
            for(j = 0; j < kval; ++j) {
                
                kernel[i][j] = kernel[i][j] / sum;
            }
        }         
        
        return;
    }
    
    static void printKernel() {
        
        int i, j;
        for(i = 0; i < kval; ++i) {
            for(j = 0; j < kval; ++j) {
                System.out.print(kernel[i][j] + " ");
            }
            System.out.println();
        }
        return;
    }
    
    static BufferedImage applyFilter(BufferedImage src){
        
        int offset = (kval / 2) + 1;
        int filteredDim = (width - kval + 1) - 10;
        
        BufferedImage filtered = new BufferedImage(filteredDim, filteredDim, BufferedImage.TYPE_BYTE_GRAY);
        
        double product = 0;
        int rgb, r, g, b, grayVal = 0;
        
        System.out.println("offset: " + offset + "filteredDim: " + filteredDim);
        //the filtered image is smaller than the original based on kernel size
        //we calculate each pixel of the filtered image in these loops
        int i, j, k, l;
        for(i = offset; i < filteredDim + offset; ++i) {
            for(j = offset; j < filteredDim + offset; ++j) {
                
                //for each pixel, perform several inner products with the Gabor kernel as seen at:
                //http://demonstrations.wolfram.com/ImageKernelsAndConvolutionLinearFiltering/
                for(k = 0; k < kval; ++k) {
                    for(l = 0; l < kval; ++l) {
                    
                        //System.out.println("Getting RGB: " + (i+k) + ", "+ (j+l));
                        
                        if(((i+k) < (width-1)) && ((j+l) < (height-1))) {
                            
                            rgb = src.getRGB((i+k), (j+l));
                            //the following code would be to convert a standard RGB pixel to grayscale
                            //i have already converted into 8bit grayscale image so the claculation is simpler
                            r = (rgb >> 16) & 0xFF;
                            g = (rgb >> 8) & 0xFF;
                            b = (rgb & 0xFF);
                            grayVal = (r + g + b) / 3;
                            //grayVal = rgb & 0xFF;
                            product = product + (grayVal * kernel[k][l]);
                        }
                    }
                }
                
                //System.out.println("Gray: " + (grayVal) + "Pixel: " + product);
                filtered.setRGB(i - offset, j - offset, (int)product);
            }
        }
        
        return filtered;
    }
    
    public static void main(String[] args) {

        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
        System.out.print("Enter a filename to be filtered: ");
        String filename = null;
        try {
            filename = br.readLine();
        } catch (IOException ex) {
            Logger.getLogger(GaborFilter.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        //read in image for filtering
        BufferedImage image = null;
        try {
            image = ImageIO.read(new File(filename));
        } catch (IOException e) {
        }
        
        height = image.getHeight();     //height of source image
        width = image.getWidth();       //width of source image
        
        //convert the image into a grayscale image (not neccessary)
        BufferedImage greyImage = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);
        
        //write the new image into a file
        Graphics2D graphics = greyImage.createGraphics();
        graphics.drawImage(image, 0, 0, null);
        try {
            ImageIO.write(greyImage, "jpg", new File("greyImage.png"));
        } catch (IOException ex) {
            Logger.getLogger(GaborFilter.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        //input parameters fot Gabor kernel
        double lambda, theta, psi, sigma, gamma;
        
        //lambda = 50 takes appx 50 minutes
        lambda = 10;               //pixel range: ( < 1/5 of image length or width to prevent edge effects)
                                   //controls "density" of the band
        
        theta = 0;                 //degree range: ( 0-360, 0 is vertical )
                                   //controls "orientation" of the band
        
        psi = 0;                   //degree range: ( -180-180 ) 
                                   //controls "unknown"
        
        sigma = .56*lambda;        //Used defailt bandwidth of 1
        
        gamma = .5;                //range: (0-1] where 1 is completely circular
                                   //controls "elipticity/height" of the band
        
        //the following code was sampled from:
        //http://patrick-fuller.com/gabor-filter-image-processing-for-scientists-and-engineers-part-6/
        //Getting the required kernel size from the sigma
        //threshold = 0.005, k is minimum odd integer required
        int k = (int)Math.ceil(Math.sqrt(-2 * sigma * sigma * Math.log(0.005)));
        if(k % 2 == 1) k++;
        
        kval = 2*k+1;            //length and width of kernel (tune parameter)
        //kval = 19;
        
        //calculate Gabor kernel (stored globally)

            gaborKernel(lambda, theta, psi, sigma, gamma);
        
        
        normalizeKernel();
        //printKernel();
        
        //convolute the image matrix with the kernel
        BufferedImage filteredImage = applyFilter(image);
        
        //write the resulting filtered image to a file
        Graphics2D printResult = filteredImage.createGraphics();
        printResult.drawImage(filteredImage, 0, 0, null);
        try {
            ImageIO.write(filteredImage, "jpg", new File("filteredImage.png"));
        } catch (IOException ex) {
            Logger.getLogger(GaborFilter.class.getName()).log(Level.SEVERE, null, ex);
        }       
    }
}
