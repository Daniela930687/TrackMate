/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fiji.plugin.trackmate.filters;

import fiji.plugin.trackmate.detection.DetectionUtils;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.process.ImageProcessor;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.BenchmarkAlgorithm;
import net.imglib2.algorithm.OutputAlgorithm;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

public class ImprovedPropagatedFilter< T extends RealType< T> & NativeType< T>> extends BenchmarkAlgorithm implements OutputAlgorithm< Img< T>>, Filter {

    private double sd = 15;//Sigma espacial
    private int sr = 3;//Sigma de rango
    int height, width;

    private Img< T> output;
    final RandomAccessibleInterval< T> input;

    public ImprovedPropagatedFilter(final RandomAccessibleInterval< T> source) {
               input=source;
        if (!DetectionUtils.filterApplied){
                DetectionUtils.filterApplied=true;
           inputParametersGui(); 
       
        }
    }

    @Override
    public boolean checkInput() {
        return true;  // throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public boolean process() {

        ImagePlus imp = ImageJFunctions.wrap(input, "");
        height=imp.getHeight();
        width=imp.getWidth();
        ImageProcessor ip = imp.getChannelProcessor();

        int a, x, y;
        byte[] pixelsA = (byte[]) ip.getPixels();
        //Imagen Original
        double[][] Img = new double[height][width];
        double[][] Is = new double[height][width];
        for (a = y = 0; y < height; y++) {
            for (x = 0; x < width; x++) {
                Img[y][x] = pixelsA[a] & 0xff;
                Is[y][x] = pixelsA[a++] & 0xff;
            }
        }
        //Variables
        int alm = 0;
        alm = ((sr * 2) + 1) * ((sr * 2) + 1);
        double Zs, sum, proc, dist;
        double[] Wst = new double[alm];
        double[] Dst = new double[alm];
        double[] Rst = new double[alm];
        int b, c, d, xx, yy, m, n, i, j;

        for (y = sr; y < height - sr; y++) {
            for (x = sr; x < width - sr; x++) {
                Zs = 0;
                sum = 0;
                a = 0;
                b = 0;
                c = 0;
                m = 0;
                d = 0;
                //Procedimiento del filtro propagado
                //Obtencion de los pesos de Rst
                for (i = j = n = 0; n < sr; n++, j++, i++) {
                    for (yy = (-1 - i), xx = (-1 - j); xx <= (n); xx++) {
                        proc = (Img[y][x] - Img[y + yy][x + xx]);
                        Rst[a++] = Math.exp(-(proc * proc) / (2 * sd * sd));
                    }
                    for (yy = (-1 - i), xx = (1 + j); yy <= (n); yy++) {
                        proc = (Img[y][x] - Img[y + yy][x + xx]);
                        Rst[a++] = Math.exp(-(proc * proc) / (2 * sd * sd));
                    }
                    for (yy = (1 + i), xx = (1 + j); (-n) <= xx; xx--) {
                        proc = (Img[y][x] - Img[y + yy][x + xx]);
                        Rst[a++] = Math.exp(-(proc * proc) / (2 * sd * sd));
                    }
                    for (yy = (1 + i), xx = (-1 - j); (-n) <= yy; yy--) {
                        proc = (Img[y][x] - Img[y + yy][x + xx]);
                        Rst[a++] = Math.exp(-(proc * proc) / (2 * sd * sd));
                    }
                }
                //Obtencion de los pesos de Dst
                for (b = i = j = n = 0; n < sr; n++, j++, i++) {
                    //Para los pesos de R
                    for (yy = (-1 - i), xx = (-1 - j); xx <= (n); xx++) {
                        dist = Math.sqrt((xx * xx) + (yy * yy));
                        if ((dist % Math.sqrt(2)) == 0) {
                            proc = (Img[y + yy + 1][x + xx + 1] - Img[y + yy][x + xx]);
                            Dst[b++] = Math.exp(-(proc * proc) / (2 * sd * sd));
                        } else if ((dist % 1) == 0) {
                            proc = (Img[y + yy + 1][x + xx] - Img[y + yy][x + xx]);
                            Dst[b++] = Math.exp(-(proc * proc) / (2 * sd * sd));
                        } else if ((dist % 1) != 0 || (dist % Math.sqrt(2)) != 0) {
                            proc = (Img[y + yy + 1][x + xx] - Img[y + yy][x + xx]);
                            Dst[b++] = Math.exp(-(proc * proc) / (2 * sd * sd));
                        }
                    }
                    for (yy = (-1 - i), xx = (1 + j); yy <= (n); yy++) {
                        dist = Math.sqrt((xx * xx) + (yy * yy));
                        if ((dist % Math.sqrt(2)) == 0) {
                            proc = (Img[y + yy + 1][x + xx - 1] - Img[y + yy][x + xx]);
                            Dst[b++] = Math.exp(-(proc * proc) / (2 * sd * sd));
                        } else if ((dist % 1) == 0) {
                            proc = (Img[y + yy][x + xx - 1] - Img[y + yy][x + xx]);
                            Dst[b++] = Math.exp(-(proc * proc) / (2 * sd * sd));
                        } else if ((dist % 1) != 0 || (dist % Math.sqrt(2)) != 0) {
                            proc = (Img[y + yy][x + xx - 1] - Img[y + yy][x + xx]);
                            Dst[b++] = Math.exp(-(proc * proc) / (2 * sd * sd));
                        }
                    }

                    for (yy = (1 + i), xx = (1 + j); (-n) <= xx; xx--) {
                        dist = Math.sqrt((xx * xx) + (yy * yy));
                        if ((dist % Math.sqrt(2)) == 0) {
                            proc = (Img[y + yy - 1][x + xx - 1] - Img[y + yy][x + xx]);
                            Dst[b++] = Math.exp(-(proc * proc) / (2 * sd * sd));
                        } else if ((dist % 1) == 0) {
                            proc = (Img[y + yy - 1][x + xx] - Img[y + yy][x + xx]);
                            Dst[b++] = Math.exp(-(proc * proc) / (2 * sd * sd));
                        } else if ((dist % 1) != 0 || (dist % Math.sqrt(2)) != 0) {
                            proc = (Img[y + yy - 1][x + xx] - Img[y + yy][x + xx]);
                            Dst[b++] = Math.exp(-(proc * proc) / (2 * sd * sd));
                        }
                    }

                    for (yy = (1 + i), xx = (-1 - j); (-n) <= yy; yy--) {
                        dist = Math.sqrt((xx * xx) + (yy * yy));
                        if ((dist % Math.sqrt(2)) == 0) {
                            proc = (Img[y + yy - 1][x + xx + 1] - Img[y + yy][x + xx]);
                            Dst[b++] = Math.exp(-(proc * proc) / (2 * sd * sd));
                        } else if ((dist % 1) == 0) {
                            proc = (Img[y + yy][x + xx + 1] - Img[y + yy][x + xx]);
                            Dst[b++] = Math.exp(-(proc * proc) / (2 * sd * sd));
                        } else if ((dist % 1) != 0 || (dist % Math.sqrt(2)) != 0) {
                            proc = (Img[y + yy][x + xx + 1] - Img[y + yy][x + xx]);
                            Dst[b++] = Math.exp(-(proc * proc) / (2 * sd * sd));
                        }
                    }
                }
                //Obtencio de los pesos Wst
                for (c = 0; c < 8; c++) {
                    Wst[c] = Rst[c] * Dst[c];
                    Zs += Wst[c];
                }
                for (m = 1, n = 1, j = 0, i = d = 8; n < sr; n++, d++, m = m + 2) {
                    for (; i < (d + 2); i++) {
                        Wst[c] = Wst[j] * Rst[i] + Dst[i];
                        Zs += Wst[c++];
                    }
                    d += 2;
                    j++;
                    for (; i < d + m; i++, j++) {
                        Wst[c] = Wst[j] * Rst[i] + Dst[i];
                        Zs += Wst[c++];
                    }
                    d += m;
                    for (; i < (d + 3); i++) {
                        Wst[c] = Wst[j] * Rst[i] + Dst[i];
                        Zs += Wst[c++];
                    }
                    d += 3;
                    j++;
                    for (; i < d + m; i++, j++) {
                        Wst[c] = Wst[j] * Rst[i] + Dst[i];
                        Zs += Wst[c++];
                    }
                    d += m;
                    for (; i < (d + 3); i++) {
                        Wst[c] = Wst[j] * Rst[i] + Dst[i];
                        Zs += Wst[c++];
                    }
                    d += 3;
                    j++;
                    for (; i < d + m; i++, j++) {
                        Wst[c] = Wst[j] * Rst[i] + Dst[i];
                        Zs += Wst[c++];
                    }
                    d += m;
                    for (; i < (d + 3); i++) {
                        Wst[c] = Wst[j] * Rst[i] + Dst[i];
                        Zs += Wst[c++];
                    }
                    d += 3;
                    j++;
                    for (; i < d + m; i++, j++) {
                        Wst[c] = Wst[j] * Rst[i] + Dst[i];
                        Zs += Wst[c++];
                    }
                    d += m;
                    Wst[c] = Wst[j] * Rst[i] + Dst[i];
                    Zs += Wst[c++];
                    j++;
                    i++;
                }
                for (c = i = j = n = 0; n < sr; n++, j++, i++) {
                    for (yy = (-1 - i), xx = (-1 - j); xx <= (n); xx++) {
                        sum += Img[y + yy][x + xx] * Wst[c++];
                    }
                    for (yy = (-1 - i), xx = (1 + j); yy <= (n); yy++) {
                        sum += Img[y + yy][x + xx] * Wst[c++];
                    }
                    for (yy = (1 + i), xx = (1 + j); (-n) <= xx; xx--) {
                        sum += Img[y + yy][x + xx] * Wst[c++];
                    }
                    for (yy = (1 + i), xx = (-1 - j); (-n) <= yy; yy--) {
                        sum += Img[y + yy][x + xx] * Wst[c++];
                    }
                }
                Is[y][x] = sum / Zs;
            }
        }
        //Mostra imagen final
        for (a = y = 0; y < width; y++) {
            for (x = 0; x < height; x++) {
                pixelsA[a++] = (byte) Is[y][x];
            }
        }
           output=ImageJFunctions.wrap(imp);

        return true;//   throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public Img< T> getResult() {
        return output;
    }

    @Override
    public void inputParametersGui() {
        GenericDialog gd = new GenericDialog("Improved Propagated");
        gd.addNumericField("Valor de la desviaci�n(sigma sd):", sd, 0);
        gd.addNumericField("Tamano de la mascara(rando sr):", sr, 0);
        gd.showDialog();
        sd = (double) gd.getNextNumber();
        sr = (int) gd.getNextNumber();
    }

}
