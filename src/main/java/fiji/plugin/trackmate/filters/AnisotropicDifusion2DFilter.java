/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fiji.plugin.trackmate.filters;

import fiji.plugin.trackmate.detection.DetectionUtils;
import ij.ImagePlus;
import ij.ImageStack;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.BenchmarkAlgorithm;
import net.imglib2.algorithm.OutputAlgorithm;
import net.imglib2.img.Img;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import ij.gui.GenericDialog;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;
import java.awt.Rectangle;
import net.imglib2.img.display.imagej.ImageJFunctions;

public class AnisotropicDifusion2DFilter< T extends RealType< T> & NativeType< T>> extends BenchmarkAlgorithm implements OutputAlgorithm< Img< T>>, Filter {

    private int nb_iter = 20;    // Number of iterations
    private int nb_smoothings = 1;     // Number of smoothings per iteration
    private double dt = 20.0;  // Adapting time step
    private double a1 = 0.5f;  // Diffusion limiter along minimal variations
    private double a2 = 0.9f;  // Diffusion limiter along maximal variations
    private int save = 20;    // Iteration saving step
    private boolean sstats = false; // display xdt value in each iteration
    private boolean tstats = false; // measure needed runtime
    private boolean add_labels = false; // add labels to output stack
    private double edgeheight = 5;     // edge threshold

    private Img< T> output;
    final RandomAccessibleInterval< T> input;
    
    
    public AnisotropicDifusion2DFilter(final RandomAccessibleInterval< T> source) {
        input=source;
 
        
    }


    @Override
    public boolean checkInput() {
        return true;//throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public boolean process() {
        if (!DetectionUtils.filterApplied){
             DetectionUtils.filterApplied=true;
           inputParametersGui(); 
          
        }
 
                
        ImagePlus imp = ImageJFunctions.wrap( input,"" );
        
        ImageProcessor ip=imp.getChannelProcessor();
        
         int scount;   // number of stacks
        ImagePlus imp2=null;   //instruccion de inicio
       ImageStack stack,stack2;
    stack=imp.getStack();
        scount=stack.getSize();
    

        int channels = ip instanceof ColorProcessor ? 3 : 1;
        Rectangle r = ip.getRoi();   // genera un rectangulo r
        int rwidth = r.width;         // ancho del rectangulo r
        int rheight = r.height;       // alto del rectangulo r
        int width = ip.getWidth(); // obtiene el ancho e la imagen
        int height = ip.getHeight(); // obtiene la altura de la imagen
        ImageProcessor ip2;
        stack2 = imp.createEmptyStack();
        double[][][][] grad = new double[2][rwidth][rheight][channels];
        double[][][] G = new double[rwidth][rheight][3]; // must be clean for each slice
        double[][][] T = new double[rwidth][rheight][3];
        double[][][] veloc = new double[rwidth][rheight][3];
        double val1, val2;
        double vec1, vec2;
        double xdt;
        double[][][] ipf = new double[rwidth][rheight][channels];
        int[] pixel = new int[channels];  //almacena canales de la image
        double fx, fy;
        double ipfnew;
        long[] time = new long[6];
        boolean breaked = false;
        int iter = 0;
        double average, stddev, drange, drange2, drange0;
        // consts
        final double c1 = (double) (0.25 * (2 - Math.sqrt(2.0))), c2 = (double) (0.5f * (Math.sqrt(2.0) - 1));


        // slices loop // segmentos o porciones del lazo
        for (int s = 0; s++ < scount && !breaked;) {

            // convert image into the double channels
            for (int x = rwidth; x-- > 0;) // recorre en X del ancho maximo hasta 0
            {
                for (int y = rheight; y-- > 0;) // recorre la imagen Y del alto maximo hasta 0
                {
                    pixel = ip.getPixel(x + r.x, y + r.y, pixel);
                    for (int k = channels; k-- > 0;) {
                        ipf[x][y][k] = pixel[k];      //****almacena canales de la imagen*****                  
                    }                                 //*******recorre la imagen y almacena en variable pixel*********
                }
            }
            // get initial stats for later normalizing
            average = 0;
            double initial_max = ipf[0][0][0], initial_min = ipf[0][0][0], pix;
            for (int x = rwidth; x-- > 0;) {
                for (int y = rheight; y-- > 0;) {
                    for (int k = channels; k-- > 0;) // realizar primera normalizacion
                    {                               //ubicacion x,y de cada canal
                        pix = ipf[x][y][k];

                        if (pix > initial_max) {
                            initial_max = pix;
                        }

                        if (pix < initial_min) {
                            initial_min = pix;
                        }

                        average += pix;      //promedio=promedio+pix
                    }
                }
            }
            average /= (rwidth * rheight * channels); // promedio=promeio/(rwidth*rheight*channels)
            // standard deviation
            stddev = 0;
            for (int x = rwidth; x-- > 0;) {
                for (int y = rheight; y-- > 0;) // ******sacar la desviacion********
                {
                    for (int k = channels; k-- > 0;) {
                        pix = ipf[x][y][k];
                        stddev += (pix - average) * (pix - average);  //promedio=promedio+pix
                    }
                }
            }
            stddev = Math.sqrt(stddev / (rwidth * rheight * channels));

            //version 0.3 normalization
            drange = (edgeheight * stddev) / 256.0;
            drange0 = (6 * stddev) / 256.0;         //*******primera norm.0-256
            drange2 = drange * drange;


            // PDE main iteration loop   -- Partial Differential  Equations main iteration loop
            for (iter = 0; (iter < nb_iter) && (!breaked); iter++) {


                double Ipp, Icp, Inp = 0, Ipc, Icc, Inc = 0, Ipn, Icn, Inn = 0;
                //double Ipp,Icp,Inp=0,Ipc,Icc,Inc=0,Ipn,Icn,Inn=0;
                /*/
                for(int k=0; k<channels; k++)
                    for(int y=0,py=0,ny=1; ny<rheight || y==--ny; py=y++,ny++)
                    {
                        Icp=Ipp=ipf[0][py][k];
                        Icc=Ipc=ipf[0][y] [k];
                        Icn=Ipn=ipf[0][ny][k];
                        for(int nx=1,x=0,px=0; nx<rwidth || x==--nx; Ipp=Icp,Ipc=Icc,Ipn=Icn,Icp=Inp,Icc=Inc,Icn=Inn,px=x++,nx++ )
                        {
                            Inp=ipf[nx][py][k];
                            Inc=ipf[nx][y] [k];
                            Inn=ipf[nx][ny][k];
                            grad[0][x][y][k]=(double)(-c1*Ipp-c2*Ipc-c1*Ipn+c1*Inp+c2*Inc+c1*Inn);
                            grad[1][x][y][k]=(double)(-c1*Ipp-c2*Icp-c1*Inp+c1*Ipn+c2*Icn+c1*Inn);
                        }
                    }
                /*/
                // the following seems several times faster
                for (int x = rwidth; x-- > 0;) {
                    int px = x - 1;
                    if (px < 0) {
                        px = 0;
                    }
                    int nx = x + 1;
                    if (nx == rwidth) {
                        nx--;
                    }
                    for (int y = rheight; y-- > 0;) {
                        int py = y - 1;
                        if (py < 0) {
                            py = 0;
                        }
                        int ny = y + 1;
                        if (ny == rheight) {
                            ny--;
                        }
                        for (int k = channels; k-- > 0;) {
                            Ipp = ipf[px][py][k];
                            Ipc = ipf[px][y][k];
                            Ipn = ipf[px][ny][k];
                            Icp = ipf[x][py][k];
                            Icn = ipf[x][ny][k];
                            Inp = ipf[nx][py][k];
                            Inc = ipf[nx][y][k];
                            Inn = ipf[nx][ny][k];
                            double IppInn = c1 * (Inn - Ipp);
                            double IpnInp = c1 * (Ipn - Inp);
                            grad[0][x][y][k] = IppInn - IpnInp - c2 * Ipc + c2 * Inc;
                            grad[1][x][y][k] = IppInn + IpnInp - c2 * Icp + c2 * Icn;
                        }
                    }
                }

                // compute structure tensor field G
//				G=new double[rwidth][rheight][3]; // must be clean for each slice
                for (int x = rwidth; x-- > 0;) {
                    for (int y = rheight; y-- > 0;) {
                        G[x][y][0] = 0.0f;
                        G[x][y][1] = 0.0f;
                        G[x][y][2] = 0.0f;
                        for (int k = channels; k-- > 0;) {   //version 0.2 normalization
                            fx = grad[0][x][y][k];
                            fy = grad[1][x][y][k];
                            G[x][y][0] += fx * fx;
                            G[x][y][1] += fx * fy;
                            G[x][y][2] += fy * fy;
                        }
                    }
                }

                // compute the tensor field T,used to drive the diffusion
                for (int x = rwidth; x-- > 0;) {
                    for (int y = rheight; y-- > 0;) {
                        // eigenvalues:
                        double a = G[x][y][0], b = G[x][y][1], c = G[x][y][1], d = G[x][y][2], e = a + d;
                        double f = Math.sqrt(e * e - 4 * (a * d - b * c));
                        double l1 = 0.5 * (e - f), l2 = 0.5 * (e + f);
                        // more precise computing of quadratic equation
                        if (e > 0) {
                            if (l1 != 0) {
                                l2 = (a * d - b * c) / l1;
                            }
                        } else {
                            if (l2 != 0) {
                                l1 = (a * d - b * c) / l2;
                            }
                        }
                        val1 = l2 / drange2;
                        val2 = l1 / drange2;
                        // slight cheat speedup for default a1 value
                        double f1 = (a1 == .5) ? 1 / Math.sqrt(1.0f + val1 + val2) : Math.pow(1.0f + val1 + val2, -a1);
                        double f2 = Math.pow(1.0f + val1 + val2, -a2);
                        // eigenvectors:
                        double u, v, n;
                        if (Math.abs(b) > Math.abs(a - l1)) {
                            u = 1;
                            v = (l1 - a) / b;
                        } else {
                            if (a - l1 != 0) {
                                u = -b / (a - l1);
                                v = 1;
                            } else {
                                u = 1;
                                v = 0;
                            }
                        }
                        n = Math.sqrt(u * u + v * v);
                        u /= n;
                        v /= n;
                        vec1 = u;
                        vec2 = v;
                        double vec11 = vec1 * vec1, vec12 = vec1 * vec2, vec22 = vec2 * vec2;
                        T[x][y][0] = f1 * vec11 + f2 * vec22;
                        T[x][y][1] = (f1 - f2) * vec12;
                        T[x][y][2] = f1 * vec22 + f2 * vec11;
                    }
                }
 
                xdt = 0.0;
                // multiple smoothings per iteration
                for (int sit = 0; sit < nb_smoothings && !breaked; sit++) {
                    // compute the PDE velocity and update the iterated image
                    //Inp=Inc=Inn=0; 
                    /*/
                    for(int k=0; k<channels; k++)
                        for(int y=0,py=0,ny=1; ny<rheight || y==--ny; py=y++,ny++)
                        {
                            Icp=Ipp=ipf[0][py][k];
                            Icc=Ipc=ipf[0][y] [k];
                            Icn=Ipn=ipf[0][ny][k];
                            for(int nx=1,x=0,px=0; nx<rwidth || x==--nx; Ipp=Icp,Ipc=Icc,Ipn=Icn,Icp=Inp,Icc=Inc,Icn=Inn,px=x++,nx++ )
                            {
                                Inp=ipf[nx][py][k];
                                Inc=ipf[nx][y] [k];
                                Inn=ipf[nx][ny][k];
                                double a=T[x][y][0],
                                    b=T[x][y][1],
                                    c=T[x][y][2],
                                    ixx=Inc+Ipc-2*Icc,
                                    iyy=Icn+Icp-2*Icc,
                                    ixy=0.5f*(Ipp+Inn-Ipn-Inp);
                                veloc[x][y][k]=a*ixx + b*ixy + c*iyy; 
                            }
                        }
                    /*/
                    // the following seems several times faster
                    for (int x = rwidth; x-- > 0;) {
                        int px = x - 1;
                        if (px < 0) {
                            px = 0;
                        }

                        int nx = x + 1;
                        if (nx == rwidth) {
                            nx--;
                        }

                        for (int y = rheight; y-- > 0;) {
                            int py = y - 1;
                            if (py < 0) {
                                py = 0;
                            }

                            int ny = y + 1;
                            if (ny == rheight) {
                                ny--;
                            }

                            for (int k = channels; k-- > 0;) {
                                Ipp = ipf[px][py][k];
                                Ipc = ipf[px][y][k];
                                Ipn = ipf[px][ny][k];
                                Icp = ipf[x][py][k];
                                Icc = ipf[x][y][k];
                                Icn = ipf[x][ny][k];
                                Inp = ipf[nx][py][k];
                                Inc = ipf[nx][y][k];
                                Inn = ipf[nx][ny][k];
                                double ixx = Inc + Ipc - 2 * Icc,
                                        iyy = Icn + Icp - 2 * Icc,
                                        ixy = 0.5f * (Ipp + Inn - Ipn - Inp);
                                veloc[x][y][k] = T[x][y][0] * ixx + T[x][y][1] * ixy + T[x][y][2] * iyy;
                            }
                        }
                    }
                    if (dt > 0) {
                        double max = veloc[0][0][0], min = veloc[0][0][0];
                        for (int x = rwidth; x-- > 0;) {
                            for (int y = rheight; y-- > 0;) {
                                for (int k = channels; k-- > 0;) {
                                    if (veloc[x][y][k] > max) {
                                        max = veloc[x][y][k];
                                    }

                                    if (veloc[x][y][k] < min) {
                                        min = veloc[x][y][k];
                                    }
                                }
                            }
                        }
                        xdt = dt / Math.max(Math.abs(max), Math.abs(min)) * drange0;
                    } else {
                        xdt = -dt;
                    }


                    // update image // actualizar la imagen
                    for (int x = rwidth; x-- > 0;) {
                        for (int y = rheight; y-- > 0;) {
                            for (int k = channels; k-- > 0;) {
                                ipfnew = ipf[x][y][k] + veloc[x][y][k] * xdt;
                                ipf[x][y][k] = ipfnew;
                                // normalize image to the original range
                                if (ipf[x][y][k] < initial_min) {
                                    ipf[x][y][k] = initial_min;
                                }

                                if (ipf[x][y][k] > initial_max) {
                                    ipf[x][y][k] = initial_max;
                                }
                            }
                        }
                    }
          
                }
                // save result
                if ((scount == 1 && ((save > 0 && (iter + 1) % save == 0) || breaked)) || (iter == nb_iter - 1)) // create new image if you break the cycle with Esc,if it is saving step,or if we get below stopping condition
                {
                    ip2 = ip.createProcessor(width, height);
                    for (int x = width; x-- > 0;) {
                        for (int y = height; y-- > 0;) {
                            for (int k = channels; k-- > 0;) {
                                if ((x < r.x) || (x >= r.x + rwidth) || (y < r.y) || (y >= r.y + rheight)) {
                                    ip2.putPixel(x, y, ip.getPixel(x, y));
                                } else {
                                    pixel[k] = (int) ipf[x - r.x][y - r.y][k];
                                    ip2.putPixel(x, y, pixel);
                                }
                            }
                        }
                    }
                    if (scount == 1) {

                        if (imp2 == null) {
                            stack2 = new ImageStack(ip2.getWidth(), ip2.getHeight());
                            stack2.addSlice("iter" + (iter + 1), ip2);
                            imp2 = new ImagePlus(imp.getShortTitle() + "-iter" + (iter + 1) + ((breaked) ? "-interrupted" : ""), stack2);
                        //    imp2.show();
                        } else {
                            stack2.addSlice("iter" + (iter + 1), ip2);
                            imp2.setStack(null, stack2);
                            imp2.setSlice(imp2.getStackSize());
                            imp2.setTitle(imp.getShortTitle() + "-iter" + (iter + 1));
                        }
                    } else {
                        stack2.addSlice(stack.getSliceLabel(s), ip2);
                    }
                }



            }
        }
        if (imp2 != null && imp2.getStackSize() > 1) {
            ip2 = imp.getProcessor().duplicate();

            stack2.addSlice("orig", ip2, 0);
            imp2.setStack(null, stack2);
            imp2.setSlice(imp2.getStackSize());
        }

        if (scount > 1) {
            new ImagePlus(imp.getShortTitle() + "-iter" + nb_iter, stack2).show();
        }

        output=ImageJFunctions.wrap(imp2);

        return true;
    }

    @Override
    public Img< T> getResult() {
        return output;
    }

    @Override
    public void inputParametersGui() {
        GenericDialog gd = new GenericDialog("2D Anisotropic Diffusion Tschumperle-Deriche v");
        gd.addNumericField("Number of iterations", nb_iter, 0);
        gd.addNumericField("Smoothings per iteration", nb_smoothings, 0);
        gd.addNumericField("a1(Diffusion limiter along minimal variations)", a1, 2);
        gd.addNumericField("a2(Diffusion limiter along maximal variations)", a2, 2);
        gd.addNumericField("dt(Time step)", dt, 1);
        gd.addNumericField("edge threshold rheight", edgeheight, 1);
        String[] labels = {"Show_filter stats", "Show_time stats", "Add labels"};
        boolean[] values = {sstats, tstats, add_labels};
        gd.addCheckboxGroup(2, 2, labels, values);
        gd.addMessage("Incorrect values will be replaced by defaults.\nLabels are drawn in the foreground color.\nPress Esc to stop processing.");
        gd.showDialog();
        // the user presses the Cancel button
        if (gd.wasCanceled()) {

        }
        nb_iter = (int) gd.getNextNumber();
        if (nb_iter < 1) {
            nb_iter = 1;
        }
        nb_smoothings = (int) gd.getNextNumber();
        if (nb_smoothings < 1) {
            nb_smoothings = 1;
        }

        a1 = (double) gd.getNextNumber();
        a2 = (double) gd.getNextNumber();
        dt = (double) gd.getNextNumber();
        edgeheight = (double) gd.getNextNumber();
        sstats = (boolean) gd.getNextBoolean();
        tstats = (boolean) gd.getNextBoolean();
        add_labels = (boolean) gd.getNextBoolean();
    }

}
