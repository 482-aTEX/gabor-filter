/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package gaborfilter;

import java.awt.Graphics2D;
import java.awt.Image;
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
    static Kernel[][] kernels;
    
    static Kernel gaborKernel(double lambda, double theta, double psi, double sigma, double gamma, Kernel kernel) {
        
 	theta = (theta * Math.PI) / (double)180;
	psi = (psi * Math.PI) / (double)180;
        
	double gauss, sinusoid, sinusoid2, xprime, yprime;
        //for each determined pixel of the kernel, we calculate its value according to equations seen:
        //http://en.wikipedia.org/wiki/Gabor_filter
        int x, y;
	for(y = -(kval/2); y < (kval/2)+1; ++y) {
		for(x = -(kval/2); x < (kval/2)+1; ++x) {
		
			xprime = (x * Math.cos(theta)) + (y * Math.sin(theta));
			yprime = -(x * Math.sin(theta)) + (y * Math.cos(theta));
                        
                                    /*System.out.println("xprime: " + xprime);
                                    System.out.println("yprime: " + yprime);
                                    System.out.println("sigma: " + sigma);
                                    System.out.println("gamma: " + gamma);
                                    System.out.println("lambda: " + lambda);*/ 
                                    
			gauss = Math.exp(-( (xprime*xprime + lambda*lambda*yprime*yprime)/(double)(2.0*sigma*sigma)));		
                        
                        sinusoid = Math.cos((2 * Math.PI * (xprime / lambda)) + psi);
                        sinusoid2 = Math.sin((2 * Math.PI * (xprime / lambda)) + psi);
			
			//kernel.kernel[y+(kval/2)][x+(kval/2)] = Math.sqrt(Math.pow((gauss * sinusoid), 2.0) + Math.pow((gauss * sinusoid2), 2.0));
		
                        double gaborReal = Math.exp(-(Math.pow(xprime/sigma, 2) + Math.pow(y*gamma/sigma, 2))/2)*Math.cos(2*Math.PI*x/lambda + psi);
                        double gaborImage = Math.exp(-(Math.pow(xprime/sigma, 2) + Math.pow(y*gamma/sigma, 2))/2)*Math.sin(2*Math.PI*x/lambda + psi);
                        kernel.kernel[y+(kval/2)][x+(kval/2)] =  Math.sqrt(Math.pow(gaborReal, 2) + Math.pow(gaborImage, 2));
                }
	}
        return kernel;
    }
    
    static Kernel normalizeKernel(Kernel kernel) {
        
        double sum = 0f;
        int i, j;
        for(i = 0; i < kval; ++i) {
            for(j = 0; j < kval; ++j) {
                
                sum = sum + kernel.kernel[i][j];
            }
        }    
        
        //===================
        sum /= kval*kval;
        for(i = 0; i < kval; ++i) {
            for(j = 0; j < kval; ++j) {
                
                kernel.kernel[i][j] = kernel.kernel[i][j] - sum;
            }
        }
        //===================
        
        //original normalization below:
        /*for(i = 0; i < kval; ++i) {
            for(j = 0; j < kval; ++j) {
                
                kernel.kernel[i][j] = kernel.kernel[i][j] / sum;
            }
        }*/         
        
        return kernel;
    }
    
    static void printKernel(Kernel kernel) {
        
        int i, j;
        for(i = 0; i < kval; ++i) {
            for(j = 0; j < kval; ++j) {
                System.out.print(kernel.kernel[i][j] + " ");
            }
            System.out.println();
        }
        return;
    }
    
    static BufferedImage applyFilter(BufferedImage src, Kernel kernel){
        
        //int offset = (kval / 2) + 1;
        int filteredDim = (width - kval + 1);
        
        BufferedImage filtered = new BufferedImage(filteredDim, filteredDim, BufferedImage.TYPE_BYTE_GRAY);
        
        double product = 0;
        int rgb, grayVal = 0;
        
        //the filtered image is smaller than the original based on kernel size
        //we calculate each pixel of the filtered image in these loops
        int i, j, k, l;
        for(i = 0; i < filteredDim; ++i) {
            for(j = 0; j < filteredDim; ++j) {
                
                //for each pixel, perform several inner products with the Gabor kernel as seen at:
                //http://demonstrations.wolfram.com/ImageKernelsAndConvolutionLinearFiltering/
                for(k = 0; k < kval; ++k) {
                    for(l = 0; l < kval; ++l) {
                            
                            rgb = src.getRGB((i+k), (j+l));
                            grayVal = rgb & 0xFF;
                            product = product + (grayVal * kernel.kernel[k][l]);
                    }
                }
                
                //System.out.println("Gray: " + (grayVal) + "Pixel: " + product);
                filtered.setRGB(i, j, (int)product);
                product = 0;
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
        Image image = null;
        try {
            image = ImageIO.read(new File(filename));
        } catch (IOException e) {
        }
        
        height = image.getHeight(null);     //height of source image
        width = image.getWidth(null);       //width of source image
        
        
            //Image image = ImageIO.read(new File("hiresbird.jpg"));
            // Creating buffered image from the given file. NOTE: It's crucial to build the data that way!
            //BufferedImage bufferedImage = new BufferedImage(image.getWidth(null), image.getHeight(null), BufferedImage.TYPE_INT_RGB);
        
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
                                   //controls "unknown" **phase offset
        
        sigma = .56*lambda;        //Used defailt bandwidth of 1
        
        gamma = .5;                //range: (0-1] where 1 is completely circular
                                   //controls "elipticity/height" of the band **aspect ratio
        
        //the following code was sampled from:
        //http://patrick-fuller.com/gabor-filter-image-processing-for-scientists-and-engineers-part-6/
        //Getting the required kernel size from the sigma
        //threshold = 0.005, k is minimum odd integer required
        int k = (int)Math.ceil(Math.sqrt(-2 * sigma * sigma * Math.log(0.005)));
        if(k % 2 == 1) k++;
        
        //kval = 2*k+1;            //length and width of kernel (tune parameter)
        kval = 3;
        
        kernels = new Kernel[4][6];
        //calculate Gabor kernels (stored locally)
       
        for(int z = 0; z < 4; ++z) {            //6 different orientations
            for(int y = 0; y < 6; ++y) {        //4 different scales     
                Kernel ker = new Kernel(kval);
                kernels[z][y] = normalizeKernel(gaborKernel(lambda, (15*y), psi, sigma, gamma, ker)); 
            }
        }
        
        //unify 6 different orientations into one kernal stored at kernels[0][0]
        for(int i = 1; i < 6; ++i) {
            
            for(int j = 0; j < kval; ++j) {
                for(int m = 0; m < kval; ++m) {
                    
                    (kernels[0][0]).kernel[j][m] += (kernels[0][i]).kernel[j][m];
                }
            }   
        }
        
        printKernel(kernels[0][0]);

        //convolute the image matrix with each of the kernels
        /*BufferedImage[][] filtered = new BufferedImage[4][6];
        for(int z = 0; z < 4; ++z) {            //6 different orientations
            for(int y = 0; y < 6; ++y) {        //4 different scales  
                filtered[z][y] = applyFilter(greyImage, kernels[z][y]);
            }
        }
        
        double[][] energyMap = new double[4][6];
        
        int pix = filtered[0][0].getHeight();
        
        
        for(int z = 0; z < 4; ++z) {            //6 different orientations
            for(int y = 0; y < 6; ++y) {        //4 different scales  
                for(int i = 0; i < pix; ++i) {
                    for(int j = 0; j < pix; ++j) {

                        energyMap[z][y] = Math.pow(filtered[z][y].getRGB(i,j), 2.0);
                    }       
                }
            }
        }
        
        double[] orientations = new double[6];
        
        for(int z = 0; z < 4; ++z) {            //6 different orientations
            
            orientations[z] = 0;
            for(int y = 0; y < 6; ++y) {        //4 different scales  
                
                orientations[z] += energyMap[z][y];
            }
        }
        
        double sgoed = 0;
        for(int i = 1; i < 7; ++i) {
            
            sgoed += (orientations[i-1] + orientations[i%6]);
        }
        
        System.out.println("SGOED Value: " + sgoed);*/
        
        BufferedImage filteredImage = applyFilter(greyImage, kernels[0][5]);
        
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
