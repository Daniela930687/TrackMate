/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fiji.plugin.trackmate.filters;

import fiji.plugin.trackmate.detection.DetectionUtils;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.filter.Convolver;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;
import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;
import java.util.ArrayList;
import java.util.List;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.BenchmarkAlgorithm;
import net.imglib2.algorithm.OutputAlgorithm;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

public class NonLocalMeansFilter< T extends RealType< T> & NativeType< T>> extends BenchmarkAlgorithm implements OutputAlgorithm< Img< T>>, Filter {

    ImagePlus imp;
    private final int weightPrecision = 1000;// Precision of the used Weight. Values>1000 can cause Overflo
    private int width;//Width of the image
    private int height;//Height of the image
    private int size;//Size image
    private int[][] pixels;// Pixels of the Image. First Dimension is for Colour-Channel,second for Pixels. Pixels are arranged in width*height
//    private int[][] pixelsExpand;// Pixels of the expanded Version of the Image. This Version is needed to prevent Walking out of Bounds.
    private int widthE;// Width of the expanded Image
    private int heightE;// Height of the expanded Image
    private int w;// Big Search-Window
    private int n;// Small Search-Window used for Patches
    private double sigma2;// Variance of the Image-Noise
    private double h2;// Smoothing-Parameter
    private int distConst;// Constant Value that the Distance gets Multiplied with. Currently unused.
    private int dim;// Dimension of the Image. (1=Grayscale,3=RGB)
    private int nextdx;// Variable used to store the next Value for dx. 
    private int nextdy;// Variable used to store the next Value for dy. 
    private boolean autoEstimate = false;// Use Auto-Estimation to determine Image-Noise
    private int sigma = 15;// Standard-Value for Sigma
    private int smoothingFactor = 1;
    private static final int BIG_ENOUGH_INT = 16 * 1024;
    private static final double BIG_ENOUGH_ROUND = BIG_ENOUGH_INT + 0.5;

    private Img< T> output;

    final RandomAccessibleInterval< T> input;

    public NonLocalMeansFilter(final RandomAccessibleInterval< T> source) {
        input = source;
        if (!DetectionUtils.filterApplied) {
            DetectionUtils.filterApplied = true;
            inputParametersGui();

        }
    }

    @Override
    public boolean process() {

        ImagePlus imp = ImageJFunctions.wrap(input, "");
        ImageProcessor ip = imp.getChannelProcessor();

              width = imp.getWidth();
        height = imp.getHeight();
        size = height * width;

        if (autoEstimate) {//Implements the gaussian noise level estimation algorithm of Immerkaer,
            //J.,1996. Fast noise variance estimation. Computer Vision and Image
            //Understanding,64(2),pp.300–302.
            FloatProcessor fp = null;
            switch (ip.getBitDepth()) {
                case 8:
                    ByteProcessor bp = (ByteProcessor) ip;
                    fp = bp.duplicate().convertToFloatProcessor();
                    break;
                case 24:
                    ColorProcessor cp = (ColorProcessor) ip;
                    fp = cp.duplicate().convertToFloatProcessor();
                    break;
                case 16:
                    ShortProcessor sp = (ShortProcessor) ip;
                    fp = sp.duplicate().convertToFloatProcessor();
                    break;
                case 32:
                    fp = (FloatProcessor) ip.duplicate();
                    break;
                default:
                    break;
            }
            Convolver convolver = new Convolver();
            //Generate kernel{1,-2,1,-2,4,-2,1,-2,1}
            int k = 6;
            int n = 2 * k + 1;
            float[] kernel = new float[n * n];
            kernel[0] = 1;
            kernel[k] = -2;
            kernel[2 * k] = 1;
            kernel[k * n] = -2;
            kernel[k * n + k] = 4;
            kernel[k * n + 2 * k] = -2;
            kernel[2 * k * n] = 1;
            kernel[2 * k * n + k] = -2;
            kernel[2 * k * n + 2 * k] = 1;
            //------------------------------------------------------------------
            convolver.convolve(fp, kernel, 2 * k + 1, 2 * k + 1);
            double sum = 0;
            double sub = 2 * k;
            float[] pixelsF = (float[]) fp.getPixels();
            for (int q = 0; q < size; q++) {
                sum += Math.abs(pixelsF[q]);
            }
            sigma = (int) (Math.sqrt(Math.PI / 2) * 1.0 / (6.0 * (width - sub) * (height - sub)) * sum);//sigma;
            IJ.log("Sigma=" + sigma);
        }
        sigma = smoothingFactor * sigma;
        //----------------------------------------------------------------------
        //Apply Non local means
        //initSettings(sigma,ip);
        //Initialize needed Settings
        int type = imp.getType();
        // Init recommended Algorithm Settings
        double hfactor;
        if (type == ImagePlus.COLOR_256 || type == ImagePlus.COLOR_RGB)// Color Image
        {
            if (sigma > 0 && sigma <= 25) {
                n = 1;
                w = 10;
//                n=3;
//                w=17;
                hfactor = 0.55;
            } else if (sigma > 25 && sigma <= 55) {
                n = 2;
                w = 17;
                hfactor = 0.4;
            } else {
                n = 3;
                w = 17;
                hfactor = 0.35;
            }
        } else//Gray Image
        if (sigma > 0 && sigma <= 15) {
            n = 1;
            w = 10;
            hfactor = 0.4;
        } else if (sigma > 15 && sigma <= 30) {
            n = 2;
            w = 10;
            hfactor = 0.4;
        } else if (sigma > 30 && sigma <= 45) {
            n = 3;
            w = 17;
            hfactor = 0.35;
        } else if (sigma > 45 && sigma <= 75) {
            n = 4;
            w = 17;
            hfactor = 0.35;
        } else {
            n = 5;
            w = 17;
            hfactor = 0.3;
        }
        widthE = width + 2 * w + 2 * n;
        heightE = height + 2 * w + 2 * n;
        // Read Pixels from ImageProcessor and store them in pixels Array
        //Converts the Image into its proper form so that it can be used by the Algorithm
        switch (type) {
            case ImagePlus.COLOR_256:
                dim = 1;
                byte[] pixelsC = (byte[]) ip.getPixels();
                pixels = new int[dim][size];
                for (int i = 0; i < size; i++) {
                    pixels[0][i] = (int) pixelsC[i] & 0xff;
                }
                break;
            case ImagePlus.COLOR_RGB:
                dim = 3;
                int[] pixelArrayC = (int[]) ip.getPixels();
                pixels = new int[dim][size];
                for (int i = 0; i < size; i++) {
                    int qtemp = pixelArrayC[i];
                    pixels[0][i] = (qtemp & 0xff0000) >> 16;
                    pixels[1][i] = (qtemp & 0x00ff00) >> 8;
                    pixels[2][i] = qtemp & 0x0000ff;
                }
                break;
            case ImagePlus.GRAY16:
                dim = 1;
                short[] pixelsS = (short[]) ip.getPixels();
                pixels = new int[dim][size];
                for (int i = 0; i < size; i++) {
                    pixels[0][i] = (int) (pixelsS[i] & (0xffff));
                }
                break;
            case ImagePlus.GRAY32:
                dim = 1;
                float[] pixelF = (float[]) ip.getPixels();
                pixels = new int[dim][size];
                for (int i = 0; i < size; i++) {
                    pixels[0][i] = (int) pixelF[i];
                }
                break;
            case ImagePlus.GRAY8:
                dim = 1;
                byte[] pixelsB = (byte[]) ip.getPixels();
                pixels = new int[dim][size];
                for (int i = 0; i < size; i++) {
                    pixels[0][i] = (int) (pixelsB[i] & 0xff);
                }
                break;
        }
        double h = hfactor * sigma;
        sigma2 = sigma * sigma * 2 * (dim * (2 * n + 1) * (2 * n + 1));
//        sigma2=2*sigma*sigma;
        distConst = (dim * (2 * n + 1) * (2 * n + 1));
//        h2=(h*h)/(dim*(2*n+1)*(2*n+1));
        h2 = (h * h);
        // Multithreadding related Initializations
        nextdx = -w;
        nextdy = -w;
        //----------------------------------------------------------------------
        try {
            int windowWidth = 512;
            int windowHeight = 512;
            double[][] result = NLMeansDenoising(windowWidth, windowHeight);
            //------------------------------------------------------------------
            //Converts a denoised Picture back to its original Format and saves it in
            //the ImageProcessor
            switch (ip.getBitDepth()) {
                case 8:
                    byte[] pixelsB = (byte[]) ip.getPixels();
                    for (int i = 0; i < size; i++) {
                        pixelsB[i] = (byte) (result[0][i]);
                    }
                    break;
                case 24:
                    int[] pixelsC = (int[]) ip.getPixels();
                    for (int i = 0; i < size; i++) {
                        pixelsC[i] = (((int) result[0][i] & 0xff) << 16) + (((int) result[1][i] & 0xff) << 8) + ((int) result[2][i] & 0xff);
                    }
                    break;
                case 16:
                    short[] pixelsS = (short[]) ip.getPixels();
                    for (int i = 0; i < size; i++) {
                        pixelsS[i] = (short) (result[0][i]);
                    }
                    break;
                case 32:
                    float[] pixelsF = (float[]) ip.getPixels();
                    for (int i = 0; i < size; i++) {
                        pixelsF[i] = (float) (result[0][i]);
                    }
                    break;
            }
        } catch (InterruptedException e) {
            e.printStackTrace();
            //IJ.showMessage("Error while computing Denoised Image.");
        }

        output = ImageJFunctions.wrap(imp);

        return true;
    }

    private double[][] NLMeansDenoising(int windowWidth, int windowHeight) throws InterruptedException {
        double[][] result = new double[dim][size];

        for (int ys = 0; ys < height; ys += windowHeight) {
            for (int xs = 0; xs < width; xs += windowWidth) {
                int imagePartWidth = (windowWidth + xs > width) ? windowWidth - ((windowWidth + xs) - width) : windowWidth;
                int imagePartHeight = (windowHeight + ys > height) ? windowHeight - ((windowHeight + ys) - height) : windowHeight;
                int[][] imagePartE = expandImage(pixels, xs, ys, imagePartWidth,
                        imagePartHeight, width, height, false);

//                double[][] partResult=NLMeansMultithreadInstance(imagePartE,
//                        Runtime.getRuntime().availableProcessors(),imagePartWidth,
//                        imagePartHeight);
                double[][] partResult = NLMeansMultithreadInstance(imagePartE,
                        Runtime.getRuntime().availableProcessors(), imagePartWidth, imagePartHeight);
                // save Partial Result in Image
                nextdx = -w;
                nextdy = -w;
                int ystart = ys;
                int xstart = xs;
                for (int y = ystart; y < ystart + imagePartHeight; y++) {
//                    if (y >=height) continue;
                    for (int x = xstart; x < xstart + imagePartWidth; x++) {
//                        if (x >=width) continue;
                        for (int d = 0; d < dim; d++) {
                            result[d][y * width + x] = partResult[d][(y - ystart) * imagePartWidth + x - xstart];
                        }
                    }
                }
            }
        }
        return result;
    }

    /**
     * Multi Threaded Implementation of the Non-local Means Algorithm. This
     * accelerated Version is based of: Darbon,Jérôme,et al. "Fast nonlocal
     * filtering applied to electron cryomicroscopy." Biomedical Imaging: From
     * Nano to Macro,2008. ISBI 2008. 5th IEEE International Symposium on.
     * IEEE,2008.
     *
     * @param image The image as Integer Array. Colors are stored within first
     * dimension of Array. Gets computed via convertImage()
     * @param threadcount Number of Threads used for Denoising
     * @param ip ImageProcessor for the original Image
     * @throws InterruptedException
     */
    private double[][] NLMeansMultithreadInstance(int[][] image, int threadcount,
            int width, int height) throws InterruptedException {
        int widthE = width + 2 * w + 2 * n;
        int heightE = height + 2 * w + 2 * n;
        long[][] u = new long[dim][widthE * heightE];
        long[] wMaxArr = new long[widthE * heightE];
        long[] wSumArr = new long[widthE * heightE];

        List<Worker> workerList = new ArrayList<Worker>(threadcount);
        for (int i = 0; i < threadcount; i++) {
            Worker worker = new Worker(width, height, image, u, wMaxArr, wSumArr);
            worker.start();
            workerList.add(worker);
        }
        for (Worker worker : workerList) {
            worker.join();
        }
        return finishPicture(u, image, wMaxArr, wSumArr, width, height);
    }

    private synchronized void deliverImagePart(long[][] imagePart, long[][] u, int widthE, int heightE) {
        for (int y = 0; y < heightE; y++) {
            int offset = y * widthE;
            for (int x = 0; x < widthE; x++) {
                for (int d = 0; d < dim; d++) {
                    u[d][offset + x] += imagePart[d][offset + x];
                }
            }
        }
    }

    /**
     * This Method is used to deliver a partial result of the Weight Sum Array.
     * The Weight Sum Array stores the sum of all Weights that are used for each
     * pixel. It is used within finishPicture(...) to properly Scale each Pixel.
     *
     * @param arr Weight Sum Array
     */
    private synchronized void deliverWSumArr(long[] arr, long[] wSumArr, int widthE, int heightE) {
        for (int y = 0; y < heightE; y++) {
            int offset = y * widthE;
            for (int x = 0; x < widthE; x++) {
                wSumArr[offset + x] += arr[offset + x];
            }
        }
    }

    /**
     * This Method is used to deliver a partial result of the Weight Max Array.
     * The Weight Max Array stores the maximum Weight that is used per Pixel.
     * This Weight is used as Weight between the Pixel and itself.
     *
     * @param arr Maximum Weight Array
     */
    private synchronized void deliverWMaxArr(long[] arr, long[] wMaxArr, int widthE, int heightE) {
        for (int y = 0; y < heightE; y++) {
            int offset = y * widthE;
            for (int x = 0; x < widthE; x++) {
                if (wMaxArr[offset + x] < arr[offset + x]) {
                    wMaxArr[offset + x] = arr[offset + x];
                }
            }
        }
    }

    /**
     * Finishes the Picture by dividing every Pixel with the Sum of all Weights
     * for the respective Pixel,and by performing the last denoising step. As
     * last Step,the Pixels get weighted with the maximum Weight for each Pixel.
     *
     * @param picture The Denoised Picture
     * @param wMaxArr Array with highest used Weight for each Pixel
     * @param wSumArr Array with Sum of Weights for each Pixel
     * @return
     */
    private double[][] finishPicture(long[][] picture, int[][] pixelsExpand,
            long[] wMaxArr, long[] wSumArr, int width, int height) {
        double[][] result = new double[dim][width * height];
        int wn = w + n;
        int widthE = width + 2 * wn;

        // x and y coordinates are based off the original Image (NOT the expanded Image)
        for (int y = 0; y < height; y++) {
            int offset = y * width;// y offset for original Image coordinates
            int offset2 = (y + wn) * widthE;// y offset for expanded Image coordinates
            for (int x = 0; x < width; x++) {
                int k = offset + x;// array Position for Pixel x,y
                int kwn = offset2 + x + wn;// same as k,but in expanded Image coordinates
                for (int d = 0; d < result.length; d++) {
                    result[d][k] = picture[d][kwn];
                    if (wMaxArr[kwn] == 0) {
                        // If Sum of all Weights is 0,just copy the original Pixel
                        result[d][k] += pixelsExpand[d][kwn];
                    } else {
                        // Weight the original Pixel with the maximum Weight
                        result[d][k] += pixelsExpand[d][kwn] * wMaxArr[kwn];
                        wSumArr[kwn] += wMaxArr[kwn];

                        // Divide Pixel by sum of all Weights
                        result[d][k] /= wSumArr[kwn];
                    }
                }
            }
        }
        return result;
    }

    private void denoise(long[][] targetArr, int[][] pixelsExpand, long[][] S,
            long[] wMaxArr, long[] wSumArr, int widthE, int heightE, int dx, int dy) {
        int wn = w + n;
        for (int y = wn; y < heightE - wn; y++) {
            int offset = y * widthE;
            int offsetn = (y + dy) * widthE;
            for (int x = wn; x < widthE - wn; x++) {
                int k = offset + x;
                int kn = offsetn + x + dx;
                int weight = computeWeight(S, widthE, x, y, weightPrecision);
                wMaxArr[k] = Math.max(weight, wMaxArr[k]);
                wSumArr[k] += weight;
                wMaxArr[kn] = Math.max(weight, wMaxArr[kn]);
                wSumArr[kn] += weight;
                for (int d = 0; d < dim; d++) {
                    int wk = weight * pixelsExpand[d][k];
                    int wkn = weight * pixelsExpand[d][kn];

                    targetArr[d][k] += wkn;
                    targetArr[d][kn] += wk;
                }
            }
        }
    }

    /**
     * Computes the Weight between the Pixel x,y and the Pixel that lies at x+
     * dx,y+dy. dx and dy are implicitly given because the Difference Image is
     * based on them.
     *
     * @param S Difference Image for a dx/dy pair
     * @param x x-Coordinate of the current Pixel
     * @param y y-Coordinate of the current Pixel
     * @param precision Precision of the Weight. Should be multiple of 10
     * @return
     */
    private int computeWeight(long[][] S, int widthE, int x, int y, int precision) {
        //Computes the Difference (distance) between the Surroundings of the Pixel x,y and the
        //Pixel that lies at x+dx,y+dy. dx and dy are implicitly given because
        //the Difference Image is based on them. Is used to compute the Weights.
        double distance = 0;
        for (int d = 0; d < dim; d++) {
            distance += S[d][(y + n) * widthE + (x + n)] + S[d][(y - n) * widthE + (x - n)]
                    - S[d][(y - n) * widthE + (x + n)] - S[d][(y + n) * widthE + (x - n)];
        }
        double exp = Math.max(distance - sigma2, 0.0);
//        exp /=h2;
//        double weight=Math.exp(-exp);
        double weight = h2 / (h2 + exp);
//        int iWeight=FastMathStuff.fastRound(weight*precision)+1;
//        if (iWeight ==0) iWeight=1;
        //Return fastRound(weight*precision);
        //Taken from http://www.java-gaming.org/index.php?topic=24194.0
        return (int) (weight * precision + BIG_ENOUGH_ROUND) - BIG_ENOUGH_INT;
    }

    /**
     * Computes the Difference Image for a given dx/dy Pair. As dx and dy can be
     * negative,image needs to be expanded to prevent out of bounds errors.
     *
     * @param image Expanded Version of Original Image
     * @param targetArr Target Array in which the Difference Image gets stored
     * into
     * @param dx
     * @param dy
     */
    private void computeDifferenceImage(int[][] image, long[][] targetArr,
            int dx, int dy, int widthE, int heightE) {
        int wn = w + n;
        long temp;

        // Compute very first Pixel of Image (x=0;y=0)
        for (int d = 0; d < dim; d++) {
            temp = image[d][wn * widthE + wn] - image[d][(wn + dy) * widthE + dx + wn];
            targetArr[d][wn * widthE + wn] = temp * temp;
        }
        // Compute first Row of Image (y=0)
        int offset = wn * widthE;
        int offsetdy = (wn + dy) * widthE;
        for (int x = wn + 1; x < widthE; x++) {
            for (int d = 0; d < dim; d++) {
                temp = image[d][offset + x] - image[d][offsetdy + x + dx];
                targetArr[d][offset + x] = targetArr[d][offset + x - 1] + temp * temp;
            }
        }
        // Compute first Column of Image (x=0)
        for (int y = wn + 1; y < heightE; y++) {
            int offsety = y * widthE;
            offsetdy = (y + dy) * widthE;
            for (int d = 0; d < dim; d++) {
                temp = image[d][offsety + wn] - image[d][offsetdy + wn + dx];
                targetArr[d][offsety + wn] = targetArr[d][offsety - widthE + wn] + temp * temp;
            }
        }
        // Compute rest of the Image
        for (int y = wn + 1; y < heightE; y++) {
            offset = y * widthE;
            int offset2 = (y + dy) * widthE;
            for (int x = wn + 1; x < widthE; x++) {
                for (int d = 0; d < dim; d++) {
                    targetArr[d][offset + x] = targetArr[d][offset + x - 1];
                    targetArr[d][offset + x] += targetArr[d][offset + x - widthE];
                    targetArr[d][offset + x] -= targetArr[d][offset + x - 1 - widthE];
                    temp = image[d][offset + x] - image[d][offset2 + x + dx];
                    double temp2 = temp * temp;
                    targetArr[d][offset + x] += temp2;
                }
            }
        }
    }

    /**
     * Expands the boundaries of an image in all four directions. The new
     * content of the Image gets filled with the adjacent parts of the Image. To
     * view a Preview of this Image,use display=true
     *
     * @param image Original Image
     * @param display Display Preview of generated Image
     * @return
     */
    private int[][] expandImage(int[][] image, int xstart, int ystart, int width,
            int height, int orgWidth, int orgHeight, boolean display) {
        int heightE = height + 2 * w + 2 * n;
        int widthE = width + 2 * w + 2 * n;
        int[][] result = new int[dim][widthE * heightE];

        for (int y = 0; y < heightE; y++) {
            int yr = y - w - n + ystart;
//            if (yr >=orgHeight) yr=(ystart-w-n)+yr-orgHeight;
            if (yr >= orgHeight) {
                yr = yr - orgHeight;
            }
            if (yr < 0) {
                yr = height + yr;
            }
            int offset = y * widthE;
            int offsetr = yr * orgWidth;
            for (int x = 0; x < widthE; x++) {
                int xr = x + (xstart - w - n);
//                if (xr >=orgWidth) xr=xstart+xr-orgWidth;
                if (xr >= orgWidth) {
                    xr = xr - orgWidth;
                }
                if (xr < 0) {
                    xr = width + xr;
                }
                for (int d = 0; d < dim; d++) {
                    result[d][offset + x] = image[d][offsetr + xr];
                }
            }
        }
        if (display) {
            int[] pixelsPicture = new int[result[0].length];
            for (int y = 0; y < heightE; y++) {
                int offset = y * widthE;
                for (int x = 0; x < widthE; x++) {
                    int p = offset + x;
                    int red = (int) result[0][p];
                    int green = (int) result[1][p];
                    int blue = (int) result[2][p];
                    int pixel = ((red & 0xff) << 16) + ((green & 0xff) << 8) + (blue & 0xff);
                    pixelsPicture[p] = pixel;
                }
            }
            BufferedImage bimg = convertToImage(widthE, heightE, pixelsPicture);
            ImagePlus imp2 = new ImagePlus("Expanded Image", bimg);
            imp2.show();
        }
        return result;
    }

    /**
     * Returns next dx/dy Pair dx and dy are needed to compute a specific
     * iteration of the Algorithm. This method provides the next unused dx/dy
     * Pair to be used in a denoising Thread.
     *
     * @return dx and dy as int array,in this respective order
     */
    private synchronized int[] getNextDV() {
        if (nextdy > 0) {
            return null;
        }
        int[] result = new int[]{nextdx, nextdy};
        if (nextdx == w) {
            nextdy++;
            nextdx = -w;
        } else {
            nextdx++;
        }
        return result;
    }

    public static BufferedImage convertToImage(int width, int height, int[] pixels) {
        int wh = width * height;
        int[] newPixels = new int[wh * 3];
        for (int i = 0; i < wh; i++) {
            int rgb = pixels[i];
            newPixels[i * 3] = (rgb >> 16) & 0xFF;
            newPixels[i * 3 + 1] = (rgb >> 8) & 0xFF;
            newPixels[i * 3 + 2] = rgb & 0xFF;
        }
        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        WritableRaster raster = (WritableRaster) image.getData();
        raster.setPixels(0, 0, width, height, newPixels);
        image.setData(raster);
        return image;
    }

   @Override
    public boolean checkInput() {
     return true;//   throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    @Override
    public Img<T> getResult() {
       return output;
    }

    @Override
    public void inputParametersGui() {
        GenericDialog gd = new GenericDialog("Non-local means");
        gd.addNumericField("Sigma", sigma, 0);
//        gd.addNumericField("Window Width",512,0);
//        gd.addNumericField("Window Height",512,0);
        gd.addNumericField("Smoothing_factor", 1, 0);
        gd.addCheckbox("Auto estimate sigma", false);
        gd.addHelp("http://fiji.sc/Non_Local_Means_Denoise");
        gd.showDialog();

        sigma = (int) gd.getNextNumber();
        smoothingFactor = (int) gd.getNextNumber();
        autoEstimate = gd.getNextBoolean();
       
  
    }

    class Worker extends Thread {

        private int[][] image;
        private long[][] u;
        private long[] wMaxArr;
        private long[] wSumArr;
        int width;
        int height;

        public Worker(int width, int height, int[][] image, long[][] u, long[] wMaxArr, long[] wSumArr) {
            this.width = width;
            this.height = height;
            this.image = image;
            this.u = u;
            this.wMaxArr = wMaxArr;
            this.wSumArr = wSumArr;
        }

        @Override
        public void run() {
            int[] vec;
            int dx, dy;
            int heightE = height + 2 * w + 2 * n;
            int widthE = width + 2 * w + 2 * n;
            long[] TwMaxArr = new long[widthE * heightE];
            long[] TwSumArr = new long[widthE * heightE];
            long[][] TimagePart = new long[dim][widthE * heightE];
            long[][] TS = new long[dim][widthE * heightE];

            vec = getNextDV();
            while (vec != null) {
                dx = vec[0];
                dy = vec[1];
                if ((2 * w + 1) * dy + dx >= 0) {
                    vec = getNextDV();
                    continue;
                }
                // compute Sdx
                computeDifferenceImage(image, TS, dx, dy, widthE, heightE);
                // denoise with Sdx
                denoise(TimagePart, image, TS, TwMaxArr, TwSumArr, widthE, heightE, dx, dy);
                // get next Vector
                vec = getNextDV();
            }
            // save to global variables
            deliverImagePart(TimagePart, u, widthE, heightE);
            deliverWMaxArr(TwMaxArr, wMaxArr, widthE, heightE);
            deliverWSumArr(TwSumArr, wSumArr, widthE, heightE);
        }
    }
}
