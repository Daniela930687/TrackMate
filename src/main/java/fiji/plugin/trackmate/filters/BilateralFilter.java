/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fiji.plugin.trackmate.filters;

import fiji.plugin.trackmate.detection.DetectionUtils;
import ij.gui.GenericDialog;
import net.imglib2.Cursor; //Valores de toda la imagen original
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.BenchmarkAlgorithm;
import net.imglib2.algorithm.OutputAlgorithm;
import net.imglib2.algorithm.neighborhood.Neighborhood;
import net.imglib2.algorithm.neighborhood.RectangleShape;
import net.imglib2.algorithm.neighborhood.RectangleShape.NeighborhoodsAccessible;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.util.Util;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

public class BilateralFilter< T extends RealType< T> & NativeType< T>> extends BenchmarkAlgorithm implements OutputAlgorithm< Img< T>>,Filter {

    private static final String BASE_ERROR_MSG = "[MedianFiler2D] ";

    private final RandomAccessibleInterval< T> source;

    private Img< T> output;

    private  int radius;

    int n;//Size of the spatial window
    int nMedios;//Size of the half of the spatial window
    int sd = 1;//Sigma espacial
    int sr = 10;// Sigma Rango

    public BilateralFilter(final RandomAccessibleInterval< T> source) {
         this.source=source;
        if (!DetectionUtils.filterApplied){
           DetectionUtils.filterApplied=true;
           inputParametersGui(); 
        }
 
    }

    @Override
    public boolean checkInput() {
        if (source.numDimensions() > 3) {
            errorMessage = BASE_ERROR_MSG + " Can only operate on 1D, 2D or 3D images. Got " + source.numDimensions() + "D.";
            return false;
        }
        if (radius < 1) {
            errorMessage = BASE_ERROR_MSG + "Radius cannot be smaller than 1. Got " + radius + ".";
            return false;
        }
        return true;
    }

    @Override
    public boolean process() {
        final long start = System.currentTimeMillis();

        final T type = source.randomAccess().get().createVariable();
        final ImgFactory< T> factory = Util.getArrayOrCellImgFactory(source, type);
        this.output = factory.create(source);

        if (source.numDimensions() > 2) {
            final long nz = source.dimension(2);
            for (long z = 0; z < nz; z++) {
                final IntervalView< T> slice = Views.hyperSlice(source, 2, z);
                final IntervalView< T> outputSlice = Views.hyperSlice(output, 2, z);
                processSlice(slice, outputSlice);
            }
        } else {
            processSlice(source, output);
        }

        this.processingTime = System.currentTimeMillis() - start;
        return true;
    }

    private void processSlice(final RandomAccessibleInterval< T> in, final IterableInterval< T> out) {
        final Cursor< T> cursor = out.localizingCursor();

        n = (int) (Math.sqrt(2 * sd * sd * 3.5065) * 2 + 1);//Side of the mask calculado para una atenuacion de 3sigma n=(sqrt-2*sd^2*ln(0.03))*2+1
        if (n % 2 == 0) {
            n++;//Evita que se presenten mascaras tamaño par.
        }
        nMedios = n / 2;

        final RectangleShape shape = new RectangleShape((n - 1) / 2, false);
        final NeighborhoodsAccessible< T> nracessible = shape.neighborhoodsRandomAccessible(Views.extendZero(in));
        final RandomAccess< Neighborhood< T>> nra = nracessible.randomAccess(in);

        final int size = (int) nra.get().size();
        final double[] values = new double[size];

        int x, y;
        double[] gaussian = new double[size];
        int cont = 0;
        for (y = -nMedios; y <= nMedios; y++) {
            for (x = -nMedios; x <= nMedios; x++) {
                double aux = (double) (Math.sqrt(x * x + y * y)) / sd;
                gaussian[cont] = (double) Math.exp(-0.5 * aux * aux);
                cont++;
            }
        }

        double cte = 0; // Suma los fcfs de cada punto
        double fcfs; // Valor Gaussiano * e-(aux/2)
        double aux; // Valor del pixel - el valor del centro / en SR
        double output = 0;
        double punto = 0; // Punto centro del vector de esa ventana para cualquier tamaño de ventana

        int tam = (n * n - 1) / 2;

        // Fill buffer with median values.
        while (cursor.hasNext()) {
            cursor.fwd();

            T puntoT = cursor.get();

            nra.setPosition(cursor); // Trae valores de la ventana en el punto que se esta trabajando

            cte = 0;
            cont = 0;
            output = 0;
            double vPix;

            int index = 0;
            for (final T pixel : nra.get()) {

                values[index++] = pixel.getRealDouble();

            }

            punto = values[(size - 1) / 2];

            for (int i = 0; i < size; i++) {

                vPix = values[i];
                aux = (Math.abs((vPix) - (punto)) / sr);
                fcfs = gaussian[cont] * Math.exp(-0.5 * aux * aux);
                cte += fcfs;
                output += vPix * fcfs;
                cont++;
            }

            puntoT.setReal(output / cte);

            //Arrays.sort( values, 0, index );
            //	cursor.get().setReal( values[ ( index - 1 ) / 2 ] );
        }
    }

    @Override
    public Img< T> getResult() {
        return output;
    }

    @Override
    public void inputParametersGui() {
      
        GenericDialog gd = new GenericDialog("2D Anisotropic Diffusion Tschumperle-Deriche v");
        gd.addNumericField("Radius", radius, 0);
        gd.showDialog();
        // the user presses the Cancel button
        if (gd.wasCanceled()) {
            
        }
        radius = (int) gd.getNextNumber();

    }

}
