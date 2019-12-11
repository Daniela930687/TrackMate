/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fiji.plugin.trackmate.filters;

import fiji.plugin.trackmate.detection.DetectionUtils;
import ij.IJ;
import ij.WindowManager;
import ij.gui.GenericDialog;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.BenchmarkAlgorithm;
import net.imglib2.algorithm.OutputAlgorithm;
import net.imglib2.img.Img;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;


public class BwGuiedFilter< T extends RealType< T > & NativeType< T >> extends BenchmarkAlgorithm implements OutputAlgorithm< Img< T >>,Filter {

    
       private double sd=0.2;//Sigma espacial
    private int sr=8;// Sigma Rango
    
    private Img< T> output;
        final RandomAccessibleInterval< T> input;
    
    public BwGuiedFilter(final RandomAccessibleInterval< T > source) {
           input=source;
        if (!DetectionUtils.filterApplied){
             DetectionUtils.filterApplied=true;
           inputParametersGui(); 
          
        }
    }

    
    
    @Override
    public boolean checkInput() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public boolean process() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

 @Override
    public Img< T> getResult() {
        return output;
    }

    @Override
    public void inputParametersGui() {
        
        int[] wList = WindowManager.getIDList();
        if (wList==null||wList.length<2)
        {
            IJ.showMessage("Filtro Guiado","This plugin requires two images of the same size.");
            return;
        }
        String[] titles=new String[wList.length];
        
        GenericDialog gd=new GenericDialog("Filtro Guiado");
        gd.addChoice("Ambient image:",titles,titles[0]);
        gd.addChoice("Flash image:",titles,titles[1]);
        gd.addNumericField("Space value (sigma d):", sd, 0);
        gd.addNumericField("Range value (sigma r):", sr, 0);
        gd.showDialog();
        
        
        
    }
    
}
