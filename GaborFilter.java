/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package gaborfilter;

import java.awt.Color;
import java.awt.color.ColorSpace;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import javax.imageio.ImageIO;

/**
 *
 * @author Kyle
 */
public class GaborFilter {

    static int height, width, kval;
    static Kernel[][] kernels;
    static float[][] convolution;    
    static PrintWriter writer;
    
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

                        gauss = (float)Math.exp(-( (xprime*xprime + gamma*gamma*yprime*yprime)/(2.0*sigma*sigma)));		

                        sinusoid = (float)Math.cos((2 * Math.PI * (xprime / lambda)) + psi);
                        sinusoid2 = (float)Math.sin((2 * Math.PI * (xprime / lambda)) + psi);

                        gaborReal = gauss * sinusoid;
                        gaborImage = gauss * sinusoid2;

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
    static void printKernel(Kernel kernel) throws IOException {
        
        
        int i, j;
        for(i = 0; i < kval; ++i) {
            for(j = 0; j < kval; ++j) {
                writer.print(kernel.kernel[i][j] + ", ");
            }
            
            writer.println();
        }
        writer.println();
        writer.flush();
        return;
    }
    
    //convolve the generated kernel with the grayscale src image
    static int applyFilter(BufferedImage src, Kernel kernel){
        
        int filteredDim2 = (width - kval + 1);
        int filteredDim = (height - kval + 1);
        
        BufferedImage filtered = new BufferedImage(filteredDim, filteredDim, BufferedImage.TYPE_BYTE_GRAY);
        
        float product = 0;
        int rgb, grayVal = 0;
        
        convolution = new float[filteredDim2][filteredDim];
        
        //the filtered image is smaller than the original based on kernel size
        //we calculate each pixel of the filtered image in these loops
        int i, j, k, l;
        for(i = 0; i < filteredDim2; ++i) {
            for(j = 0; j < filteredDim; ++j) {
                
                //for each pixel, perform several inner products with the Gabor kernel as seen at:
                //http://demonstrations.wolfram.com/ImageKernelsAndConvolutionLinearFiltering/
                for(k = 0; k < kval; ++k) {
                    for(l = 0; l < kval; ++l) {
                            
                            rgb = src.getRGB((i+k), (j+l));
                            grayVal = rgb & 0xFF;
                            
                            product = product + (float)grayVal * kernel.kernel[k][l];    
                    }
                }
                
                convolution[i][j] = product;
                //System.out.println("gray: " + grayVal + " compconvol: " + product);
                product = 0;
            }
        }
        
        return filteredDim;
    }
    
    //normalize the values obtained from convolution to a 0-255 grayscale range
    static BufferedImage constructResult(int filteredDim) {
        
        BufferedImage result = new BufferedImage(filteredDim, filteredDim, BufferedImage.TYPE_BYTE_GRAY);
        
        int max = (int)convolution[0][0];
        int min = (int)convolution[0][0];
        int grayVal;
                
        for(int i = 1; i < filteredDim; ++i) {
            for(int j = 1; j < filteredDim; ++j) {
            
                grayVal = (int)convolution[i][j];
                if(grayVal < min)
                    min = grayVal;
                if(grayVal > max)
                    max = grayVal;
            }
        }
        System.out.println("Min: " + min + " Max: " + max);
                
        int newGray, rgb;
        for(int i = 0; i < filteredDim; ++i) {
            for(int j = 0; j < filteredDim; ++j) {
                
                grayVal = (int)convolution[i][j];
                
                //newGray = (grayVal-min) * ((255-0)/(max-min)) + 0;
                //rgb = newGray<<16 | newGray << 8 | newGray;
                rgb = grayVal<<16 | grayVal << 8 | grayVal;
                System.out.println("setting rgb: " + grayVal);
                //if(i < kval && j < kval) {
                    
                    //result.setRGB(i,j, kernels[0][0].kernel[i][j]);
                //}
                result.setRGB(i,j,rgb);
            }
        }
        
        return result;
    }
    
    //inport image and apply gabor filter and energy calculations
    public static void main(String[] args) throws IOException {

        float max;
        float min;
        
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
        ImageIO.write(greyImage, "png", new File("greyImage.png"));
        
        //input parameters for Gabor kernel
        float lambda, theta, psi, sigma, gamma;
        
        lambda = 2;               //pixel range: ( < 1/5 of image length or width to prevent edge effects)
                                   //controls "density" of the band
        
        theta = 0;                 //degree range: ( 0-360, 0 is vertical )
                                   //controls "orientation" of the band
        
        psi = 0;                   //degree range: ( -180-180 ) 
                                   //controls "unknown" **phase offset
        
        sigma = (float).56*lambda;        //Used defailt bandwidth of 1
        
        gamma = (float).05;                //range: (0-1] where 1 is completely circular
                                   //controls "elipticity/height" of the band **aspect ratio
        
        //the following code was sampled from:
        //http://patrick-fuller.com/gabor-filter-image-processing-for-scientists-and-engineers-part-6/
        //Getting the required kernel size from the sigma
        //threshold = 0.005, k is minimum odd integer required
        int k = (int)Math.ceil(Math.sqrt(-2 * sigma * sigma * Math.log(0.005)));
        if(k % 2 == 1) k++;
        
        kval = 2*k+1;            //length and width of kernel (tune parameter)
        //kval = 3;
        
        
        System.out.println("kernel width: " + kval);
        System.out.println("image size: " + width + "x" + height);
        
        kernels = new Kernel[4][6];
        //calculate Gabor kernels (stored locally)
       
        writer = new PrintWriter(new FileWriter("9x9kernels.txt", true));
        writer.println();
        
        for(int z = 0; z < 4; ++z) {            //4 different scales  
            for(int y = 0; y < 6; ++y) {        //6 different orientations   
                Kernel ker = new Kernel(kval);
                kernels[z][y] = gaborKernel(lambda, (30*y), psi, sigma, (gamma*(float)Math.pow(2.0, z)), ker);
                printKernel(normalizeKernel(kernels[z][y]));
                //System.out.println("===================================================================");
            }
        }
        //System.out.println("GENERATED KERNELS");
        
        //convolute the image matrix with each of the kernels
        Kernel[][] filtered = new Kernel[4][6];
        for(int z = 0; z < 4; ++z) {            //4 different scales
            for(int y = 0; y < 6; ++y) {        //6 different orientations  
                
                int newK = applyFilter(greyImage, normalizeKernel(kernels[z][y]));
                filtered[z][y] = new Kernel(newK);
                filtered[z][y].kernel = convolution;
                //printKernel(filtered[z][y]);
            }
        }
        //int pix = filtered[0][0].kval;
        int pix = width - (kval);
        int pix2 = height - (kval);
        
        double[][] energyMap = new double[4][6];
        float[] orientations = new float[6];
        float sgoed[] = new float[pix*pix2];
        float newSgoed = 0;
        float totalSgoed = 0;
        
        BufferedImage result = new BufferedImage(pix, pix, BufferedImage.TYPE_BYTE_GRAY);
        max = (float)0.0;
        min = (float)10000000.0;   
        //System.out.println("Energy Map: ");                
        for(int i = 0; i < pix; ++i) {
            for(int j = 0; j < pix2; ++j) {
                
                sgoed[j+i*pix] = 0;
                
                //===================================================================
                for(int z = 0; z < 4; ++z) {            //4 different scales  
                    for(int y = 0; y < 6; ++y) {        //6 different orientations

                        
                        energyMap[z][y] = (float)Math.pow(filtered[z][y].kernel[i][j], 2.0);
                    }       
                }//System.out.print(energyMap[z][y] + " ");
                
                
                for(int y = 0; y < 6; ++y) {            //4 different scales
            
                    orientations[y] = 0;
                    for(int z = 0; z < 4; ++z) {        //6 different orientations  

                        orientations[y] += energyMap[z][y];
                    }//System.out.print(orientations[y] + " ");
                }//System.out.println();

                for(int n = 1; n < 7; ++n) {

                    newSgoed += Math.abs((orientations[n%6] - orientations[n-1]));
                    if(newSgoed < min)
                       min = newSgoed;
                    if(newSgoed > max)
                        max = newSgoed;
                    
                } 
                //System.out.println("sgoed: " + newSgoed);
                totalSgoed += newSgoed;
                if(newSgoed > 15000) {
                    greyImage.setRGB(i, j, 0xFFFFFFFF);
                    
                    for(int q = 5; q > 0; --q) {
                        for(int t = 5; t > 0; --t) {
                            greyImage.setRGB(i-q, j-t, 0xFFFFFFFF);
                        }
                    }
                    
                }
                if(newSgoed < 15000) {
                    greyImage.setRGB(i, j, 0xFFFFFFFF);
                }
                if(newSgoed < 14000) {
                    greyImage.setRGB(i, j, 0xDDDDDDDD);
                }
                if(newSgoed < 12000) {
                    greyImage.setRGB(i, j, 0xBBBBBBBB);
                }
                if(newSgoed < 10000) {
                    greyImage.setRGB(i, j, 0x99999999);
                }
                if(newSgoed < 8000) {
                    greyImage.setRGB(i, j, 0x88888888);
                }
                if(newSgoed < 6000) {
                    greyImage.setRGB(i, j, 0x77777777);
                }
                if(newSgoed < 4000) {
                    greyImage.setRGB(i, j, 0x55555555);
                }
                if(newSgoed < 2000) {
                    greyImage.setRGB(i, j, 0x33333333);
                }
                if(newSgoed < 500)
                    greyImage.setRGB(i, j, 0x22222222);
                if(newSgoed < 400)
                    greyImage.setRGB(i, j, 0x11111111);
                if(newSgoed < 200)
                    greyImage.setRGB(i, j, 0x00000000);
                
                
                //writer.print(newSgoed + "    ");
                
                sgoed[j+i*pix] = newSgoed;
                newSgoed = 0;
                for(int y = 0; y < 6; ++y) {
                    orientations[y] = 0;
                }
        
            //=======================================================================
            }//writer.println();
        }
        
        
        /*for(int i = 0; i < pix; ++i) {
            for(int j = 0; j < pix; ++j) {
                
                System.out.println("check " + sgoed[j+i*pix]);
                int normalSgoed = ((int)sgoed[j+i*pix]-(int)min) * ((255-0)/((int)max-(int)min)) + 0;
                int rgb = normalSgoed<<16 | normalSgoed << 8 | normalSgoed;
                System.out.println("N " + normalSgoed);
                result.setRGB(i,j,rgb);
            }
        }*/
        System.out.println("min sgoed: " + min + " max sgoed: " + max);
        
        System.out.println("sum of SGOEDs: " + totalSgoed);
        
        ImageIO.write(greyImage, "png", new File("greyImage.png"));
        
        ImageIO.write(result, "jpg", new File("filteredImage.jpg"));
        //System.out.println("Orientation Energies: "); 
        //System.out.println("SGOED Value: " + sgoed);

        /*for(int i = 0; i < 6; ++i) {
            int filteredDim = applyFilter(greyImage, normalizeKernel(kernels[0][i]));
            BufferedImage result = constructResult(filteredDim);
            //write the resulting filtered image to a file
            
        }*/
    }
}