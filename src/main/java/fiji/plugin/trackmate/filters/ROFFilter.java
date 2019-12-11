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

public class ROFFilter< T extends RealType< T> & NativeType< T>> extends BenchmarkAlgorithm implements OutputAlgorithm< Img< T>>, Filter {

    private float lambda = 80;
    private int nbIterations = 10;
    int type;

    private Img< T> output;

    final RandomAccessibleInterval< T> input;

    public ROFFilter(final RandomAccessibleInterval< T> source) {
        input=source;
        if (!DetectionUtils.filterApplied){
             DetectionUtils.filterApplied=true;
           inputParametersGui(); 
          
        }
    }

    @Override
    public boolean checkInput() {
        return true;// throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public boolean process() {
        //  throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.

        ImagePlus imp = ImageJFunctions.wrap(input, "");
        ImageProcessor ip = imp.getChannelProcessor();
        final int height = ip.getHeight();
        final int width = ip.getWidth();
        final int size = width * height;
        final float dt = 0.25f;
        final float g = 1;
        final float[] u = new float[size];
        final float[] p = new float[size * 2];
        final float[] d = new float[size * 2];
        final float[] du = new float[size * 2];
        final float[] div_p = new float[size];
        int a, i, j;
        byte[] pixels8 = null;
        short[] pixels16 = null;
        float[] pixels32 = null;

        switch (type) {
            case 8:
                pixels8 = (byte[]) ip.getPixels();
                pixels32 = (float[]) ip.convertToFloatProcessor().getPixels();
                break;
            case 16:
                pixels16 = (short[]) ip.getPixels();
                pixels32 = (float[]) ip.convertToFloatProcessor().getPixels();
                break;
            case 32:
                pixels32 = (float[]) ip.getPixels();
                break;
        }
        for (int iteration = 0; iteration < nbIterations; iteration++) {
            for (i = 0; i < width; i++) {
                for (j = 1; j < height - 1; j++) {
                    div_p[i + j * width] = p[i + j * width] - p[i + width * (j - 1)];
                }
                //Handle boundaries
                div_p[i] = p[i];
                div_p[i + width * (height - 1)] = -p[i + width * (height - 1)];
            }
            for (j = 0; j < height; j++) {
                a = j * width;
                for (i = 1; i < width - 1; i++) {
                    div_p[i + a] += p[i + width * (j + height)] - p[i - 1 + width * (j + height)];
                }
                //Handle boundaries
                div_p[a] = p[width * (j + height)];
                div_p[width - 1 + a] = -p[width - 1 + width * (j + height)];
            }
            //Update u
            for (a = 0; a < size; a++) {
                u[a] = pixels32[a] - lambda * div_p[a];
            }
            //Calculate forward derivatives
            for (j = 0; j < height; j++) {
                for (i = 0; i < width; i++) {
                    a = j * width;
                    if (i < width - 1) {
                        du[i + width * (j + height)] = u[i + 1 + a] - u[i + a];
                    }
                    if (j < height - 1) {
                        du[i + a] = u[i + width * (j + 1)] - u[i + a];
                    }
                }
            }
            for (j = 0; j < height; j++) {
                a = j * width;
                for (i = 0; i < width; i++) {
                    final float du1 = du[a + i], du2 = du[i + a + width * height];
                    d[i + a] = 1 + dt / lambda / g * Math.abs((float) Math.sqrt(du1 * du1 + du2 * du2));
                    d[i + a + width * height] = 1 + dt / lambda / g * Math.abs((float) Math.sqrt(du1 * du1 + du2 * du2));
                    p[i + a] = (p[i + a] - dt / lambda * du[i + a]) / d[i + a];
                    p[i + a + width * height] = (p[i + a + width * height] - dt / lambda * du[i + a + width * height]) / d[i + a + width * height];
                }
            }
        }
        switch (type) {
            case 8:
                for (i = 0; i < size; i++) {
                    pixels8[i] = (byte) u[i];
                }
                break;
            case 16:
                for (i = 0; i < size; i++) {
                    pixels16[i] = (short) u[i];
                }
                break;
            case 32:
                System.arraycopy(u, 0, pixels32, 0, size);
                break;
        }
        output = ImageJFunctions.wrap(imp);

        return true;
    }

    @Override
    public Img< T> getResult() {
        return output;
    }

    @Override
    public void inputParametersGui() {
        final GenericDialog gd = new GenericDialog("ROF Denoise");
        gd.addNumericField("Lambda", lambda, 2);
        gd.addNumericField("Number of iterations", nbIterations, 2);
        gd.showDialog();
        if (gd.wasCanceled()) {

        }
        lambda = (float) gd.getNextNumber();
        nbIterations = (int) gd.getNextNumber();

    }

}
