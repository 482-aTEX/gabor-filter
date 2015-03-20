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
import javax.imageio.ImageIO;

/**
 *
 * @author Kyle
 */
public class GaborFilter {

    static int height, width, kval;
    static Kernel[][] kernels;
    static int[][] convolution;    
    
    //generate a gabor kernel according to the entered parameters
    static Kernel gaborKernel(float lambda, float theta, float psi, float sigma, float gamma, Kernel kernel) {
        
 	theta = (theta * (float)Math.PI) / (float)180;
	psi = (psi * (float)Math.PI) / (float)180;
        
	float gauss, sinusoid, sinusoid2,gaborReal, gaborImage, xprime, yprime;
        //for each determined pixel of the kernel, we calculate its value according to equations seen:
        //http://en.wikipedia.org/wiki/Gabor_filter
        int x, y;
	for(y = -(kval/2); y < (kval/2)+1; ++y) {
		for(x = -(kval/2); x < (kval/2)+1; ++x) {
		
			xprime = (x * (float)Math.cos(theta)) + (y * (float)Math.sin(theta));
			yprime = -(x * (float)Math.sin(theta)) + (y * (float)Math.cos(theta));
                                    
			gauss = (float)Math.exp(-( (xprime*xprime + lambda*lambda*yprime*yprime)/(2.0*sigma*sigma)));		
                        
                        sinusoid = (float)Math.cos((2 * Math.PI * (xprime / lambda)) + psi);
                        sinusoid2 = (float)Math.sin((2 * Math.PI * (xprime / lambda)) + psi);
					
                        gaborReal = gauss * sinusoid;
                        gaborImage = gauss * sinusoid2;
                        
                        //kernel.kernel[y+(kval/2)][x+(kval/2)] =  Math.sqrt(gaborReal*gaborReal + gaborImage*gaborImage);
                        kernel.kernel[y+(kval/2)][x+(kval/2)] =  gaborReal;
                }
	}
        return kernel;
    }
    
    //normalize the generated kernel so its sum is approximately 1
    static Kernel normalizeKernel(Kernel kernel) {
        
        float sum = 0f;
        int i, j;
        for(i = 0; i < kval; ++i) {
            for(j = 0; j < kval; ++j) {
                
                sum = sum + (float)Math.abs(kernel.kernel[i][j]);
            }
        }    
        
        float kcheck = 0;
        for(i = 0; i < kval; ++i) {
            for(j = 0; j < kval; ++j) {
                
                kernel.kernel[i][j] = kernel.kernel[i][j] / sum;
                kcheck += kernel.kernel[i][j];
            }
        }   
        
        return kernel;
    }
    
    //prints the kernel values
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
    
    //convolve the generated kernel with the grayscale src image
    static int applyFilter(BufferedImage src, Kernel kernel){
        
        int filteredDim = (width - kval + 1);
        
        BufferedImage filtered = new BufferedImage(filteredDim, filteredDim, BufferedImage.TYPE_BYTE_GRAY);
        
        float product = 0;
        int rgb, rgb2, grayVal = 0;
        
        convolution = new int[filteredDim][filteredDim];
        
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
                            product = product + (float)((float)grayVal * kernel.kernel[k][l]);
                            
                    }
                }
                
                rgb2 = (int)product<<16 | (int)product << 8 | (int)product;
                convolution[i][j] = (int)product;
                product = 0;
            }
        }
        
        return filteredDim;
    }
    
    //normalize the values obtained from convolution to a 0-255 grayscale range
    static BufferedImage constructResult(int filteredDim) {
        
        BufferedImage result = new BufferedImage(filteredDim, filteredDim, BufferedImage.TYPE_BYTE_GRAY);
        
        int max = convolution[0][0];
        int min = convolution[0][0];
        int grayVal;
                
        for(int i = 1; i < filteredDim; ++i) {
            for(int j = 1; j < filteredDim; ++j) {
            
                grayVal = convolution[i][j];
                if(grayVal < min)
                    min = grayVal;
                if(grayVal > max)
                    max = grayVal;
            }
        }
                
        int newGray, rgb;
        for(int i = 0; i < filteredDim; ++i) {
            for(int j = 0; j < filteredDim; ++j) {
                
                grayVal = convolution[i][j];
                newGray = (grayVal-min) * ((255-0)/(max-min)) + 0;
                rgb = newGray<<16 | newGray << 8 | newGray;
                result.setRGB(i,j,rgb);
            }
        }
        
        return result;
    }
    
    //inport image and apply gabor filter and energy calculations
    public static void main(String[] args) throws IOException {

        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
        System.out.print("Enter a filename to be filtered: ");
        String filename = null;
        filename = br.readLine();

        
        //read in image for filtering
        Image image = null;
        image = ImageIO.read(new File(filename));
        
        height = image.getHeight(null);     //height of source image
        width = image.getWidth(null);       //width of source image
        
        BufferedImage greyImage = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);
        
        //write the new image into a file
        Graphics2D graphics = greyImage.createGraphics();
        graphics.drawImage(image, 0, 0, null);
        ImageIO.write(greyImage, "jpg", new File("greyImage.png"));
        
        //input parameters for Gabor kernel
        float lambda, theta, psi, sigma, gamma;
        
        //lambda = 50 takes appx 50 minutes
        lambda = 20;               //pixel range: ( < 1/5 of image length or width to prevent edge effects)
                                   //controls "density" of the band
        
        theta = 0;                 //degree range: ( 0-360, 0 is vertical )
                                   //controls "orientation" of the band
        
        psi = 0;                   //degree range: ( -180-180 ) 
                                   //controls "unknown" **phase offset
        
        sigma = (float).56*lambda;        //Used defailt bandwidth of 1
        
        gamma = (float).5;                //range: (0-1] where 1 is completely circular
                                   //controls "elipticity/height" of the band **aspect ratio
        
        //the following code was sampled from:
        //http://patrick-fuller.com/gabor-filter-image-processing-for-scientists-and-engineers-part-6/
        //Getting the required kernel size from the sigma
        //threshold = 0.005, k is minimum odd integer required
        int k = (int)Math.ceil(Math.sqrt(-2 * sigma * sigma * Math.log(0.005)));
        if(k % 2 == 1) k++;
        
        //kval = 3;
        kval = 2*k+1;            //length and width of kernel (tune parameter)
        
        kernels = new Kernel[4][6];
        //calculate Gabor kernels (stored locally)
       
        for(int z = 0; z < 4; ++z) {            //4 different scales  
            for(int y = 0; y < 6; ++y) {        //6 different orientations   
                Kernel ker = new Kernel(kval);
                kernels[z][y] = gaborKernel(lambda, (15*y), psi, sigma, gamma, ker); 
            }
        }
        
        //printKernel(kernels[0][0]);
        
        //unify 6 different orientations into one kernal stored at kernels[0][0]
        /*for(int i = 1; i < 6; ++i) {
            
            for(int j = 0; j < kval; ++j) {
                for(int m = 0; m < kval; ++m) {
                    
                    (kernels[0][0]).kernel[j][m] += (kernels[0][i]).kernel[j][m];
                }
            }   
        }*/
        
        //printKernel(kernels[0][0]);

        //convolute the image matrix with each of the kernels
        /*BufferedImage[][] filtered = new BufferedImage[4][6];
        for(int z = 0; z < 4; ++z) {            //6 different orientations
            for(int y = 0; y < 6; ++y) {        //4 different scales  
                filtered[z][y] = applyFilter(greyImage, kernels[z][y]);
            }
        }
        
        float[][] energyMap = new float[4][6];
        
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
        
        float[] orientations = new float[6];
        
        for(int z = 0; z < 4; ++z) {            //6 different orientations
            
            orientations[z] = 0;
            for(int y = 0; y < 6; ++y) {        //4 different scales  
                
                orientations[z] += energyMap[z][y];
            }
        }
        
        float sgoed = 0;
        for(int i = 1; i < 7; ++i) {
            
            sgoed += (orientations[i-1] + orientations[i%6]);
        }
        
        System.out.println("SGOED Value: " + sgoed);*/

        int filteredDim = applyFilter(greyImage, normalizeKernel(kernels[0][0]));
        
        BufferedImage result = constructResult(filteredDim);
        
        //write the resulting filtered image to a file
        ImageIO.write(result, "jpg", new File("filteredImage.png"));
    }
}