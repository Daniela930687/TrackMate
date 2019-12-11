/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fiji.plugin.trackmate.filters;

import fiji.plugin.trackmate.detection.DetectionUtils;
import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.plugin.filter.*;
import java.util.ArrayList;
import java.util.Random;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.BenchmarkAlgorithm;
import net.imglib2.algorithm.OutputAlgorithm;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;


public class K_SVD_Filter < T extends RealType< T> & NativeType< T>> extends BenchmarkAlgorithm implements OutputAlgorithm< Img< T>>, Filter 
{
    int height, width, chnls = 1;
    ImagePlus imp;
    double img_noisy[];
    float img[];
    float img_denoised[];
    double fSigma;
    boolean useAcceleration=false;
    boolean noGo = false;
    double mean=0;
    double variance;
    double         psnr=0;
    double      rmse=0;
    int N_iter =0;
    int col;//! Number of iterations

    private Img< T> output;

    final RandomAccessibleInterval< T> input;
    
    
    public K_SVD_Filter(final RandomAccessibleInterval< T> source) {
        input=source;
        if (!DetectionUtils.filterApplied){
             DetectionUtils.filterApplied=true;
           inputParametersGui(); 
          
        }
    }
    
    @Override
    public void inputParametersGui()
    {
        
        GenericDialog gd = new GenericDialog("K-SVD Algorithm");
        String[] items2={"No","Yes"};
        gd.addNumericField("Number Iterations:", N_iter, 1);
        gd.addNumericField("Standard desviation:", fSigma, 3);
        gd.addChoice("Noise",items2,items2[0]);
        gd.showDialog();
        col=(gd.getNextChoiceIndex())*255;
        N_iter = (int) gd.getNextNumber();
        fSigma = gd.getNextNumber();

        
        if (gd.wasCanceled()) {
            if (imp != null) {
                imp.unlock();
            }
            noGo = true;
        }
        variance = fSigma * fSigma;
        if(fSigma==0)
        {
        }
        

     
    }

    @Override
    public boolean process()
    { 
        ImagePlus imp = ImageJFunctions.wrap(input, "");
        ImageProcessor ip = imp.getChannelProcessor();
        long TInicio, TFin, tiempo; //Variables para determinar el tiempo de ejecuci�n
        TInicio=System.currentTimeMillis();
        width = ip.getWidth();
        height = ip.getHeight();
        img = (float[]) ip.convertToFloatProcessor().getPixels();
        img_noisy = new double[width*height];
        if (col==255)
        {
        Add_Gaussian_Noise(img,img_noisy); 
        }
        else
        {
        for(int i=0;i<width*height;i++)
            img_noisy[i]=(double)img[i];
        }
        ImagePlus imp1=NewImage.createByteImage("image noisy ",width,height,1,NewImage.FILL_BLACK);
        ImageProcessor ip1=imp1.getProcessor();
        byte[] pixels=(byte[])ip1.getPixels();
        for(int i=0;i<width*height;i++)
                pixels[i]=(byte)img_noisy[i];
     //   imp1.show();
        img_denoised=new float[width*height];
       //----------------------------------------------------------------------
        //! Denoising
        //! Initializations
        
        double sigma=fSigma/255.0;
        int N1 = (sigma * 255.0 <= 20 ? 5 : (sigma * 255.0 <= 60 ? 7 : 9)); //! Size of patches
        int N2 = 256;    //! size of the dictionary
        int N1_2 = N1 * N1;

        int T = (sigma * 255.0 > 40.0 ? 8 : (sigma * 255.0 > 10.0 ? 16 : 32));
        double gamma = 5.25;
        int num_patches = (width - N1 + 1) * (height - N1 + 1);
        //! C = sqrt(1/(chnls * N1 * N1)*chi2inv(0.93,chnls * N1 * N1));
        double C = (chnls == 1 ? (N1 == 5 ? 1.2017729876383829 : (N1 == 7 ? 1.1456550151825420 : 1.1139195378939404)) : (N1 == 5 ? 1.1182997771678573 : (N1 == 7 ? 1.0849724948297015 : 1.0662877194412401)));
        //! Assuming that img_noisy \isin [0, 255]
        for (int k = 0; k < width * height * chnls; k++)
            img_noisy[k]/=255.0;
        //! Declarations
        double [][] patches = new double[num_patches][N1_2 * chnls];
        double [][] dictionary = new double[N2][N1_2 * chnls];
        //! Decompose the image into patches
        im2patches(img_noisy, patches, N1);
        int h=0;
        //! Keep (1 / T) patches to learn the dictionary
     if (useAcceleration || (int) (num_patches / T) > N2)
      {
            //! Obtain less patches (divided by T)
            int num_sub_patches = (int) (num_patches / T);
            double [][]  sub_patches = new double[num_sub_patches][N1_2 * chnls];
            for (int k = 0; k < num_sub_patches; k++)
                sub_patches[k] = patches[k * T];
            obtain_dict(dictionary, sub_patches);//Obtain the initial dictionary
            int h0=0;
//! Learn dictionary
          ksvd_process( sub_patches, dictionary, sigma, N1, N2, N_iter, gamma, C, false);
            int h00=0;

            //! Denoising
//          
           ksvd_process( patches, dictionary, sigma, N1, N2, 1, gamma, C, true);
                  
/*ImagePlus imp3=NewImage.createByteImage("Dictionary ",N2,N1_2,1,NewImage.FILL_BLACK);
        ImageProcessor ip3=imp3.getProcessor();
        byte[] pixels2=(byte[])ip3.getPixels();
        for(int i=0;i<25;i++){
            for(int j=0;j<25;j++)  {  
                dictionary[i][j]*=255;
                pixels2[j*N1_2+i]=(byte)dictionary[i][j];
            }
        }
        imp3.show();*/
        
            int h1=0;
      }
     else
    {
        //! Obtain the initial dictionary
        obtain_dict(dictionary, patches);

        //! Denoising
        ksvd_process( patches, dictionary, sigma, N1, N2, 1, gamma, C, true);
    }

   
       //! Back to the [0, 255] for the image value
        for (int k = 0; k < width * height * chnls; k++)
     //       img_noisy[k]*=255.0;
         img_denoised[k]*=255.0;
       
      psnr_rmse(img, img_denoised,  width*height);
       
        ImagePlus imp2=NewImage.createByteImage("image denoised ",width,height,1,NewImage.FILL_BLACK);
        ImageProcessor ip2=imp2.getProcessor();
        byte[] pixels1=(byte[])ip2.getPixels();
        for(int i=0;i<width*height;i++)
                pixels1[i]=(byte)img_denoised[i];
        //imp2.show();
        
        TFin = System.currentTimeMillis(); //Tomamos la hora en que finaliz� el algoritmo y la almacenamos en la variable T
        tiempo = TFin - TInicio; //Calculamos los milisegundos de diferencia

        
                output = ImageJFunctions.wrap(imp2);

        return true;


 
    }
void psnr_rmse(float[]   img_1, float[]   img_2, int size)
{
    float tmp = 0.0f;
    for (int k = 0; k < size; k++)
        tmp += (img_1[k] - img_2[k]) * (img_1[k] - img_2[k]);

    rmse = Math.sqrt(tmp / (float) size);
    psnr = 20.0 *(float) Math.log10(255.0 / (rmse));
}  
void Add_Gaussian_Noise(float img[],double img_noise[])
{
       
        if (noGo) 
	{
            return;
        }
        int noise;
        double value;
        Random r = new Random();
                    double st_dev = Math.sqrt(variance);

        
           boolean inRange;
           for (int i = 0; i < width*height; i++)
            {
                
                    inRange = false;
                    do
                    {
                        noise = (int) Math.round(r.nextGaussian()*st_dev);
                        value = img[i]+noise;
                        inRange = value>=0 && value<=255;
                        if (inRange) img_noise[i] = (double)value;
                    } while (!inRange);
                
            }
    }

 
    void im2patches(double img[], double[][] patches, int N ) {
        //! Declarations
        int h_p = height - N + 1;
        int wh_p = (width - N + 1) * h_p;
        
        double[][] it_p =new double[wh_p][];
        double[] it=new double[N];
        for (int k = 0; k <wh_p ; k++)
        {
            it_p[k] = patches[k];

            int dk = k / h_p + (k % h_p) * width;
            for (int c = 0; c < chnls; c++)
            {
                int dc = c * width * height + dk;
                int q=0;
                int i=0;
                for (int p = 0; p < N; p++, dc++)
                {
                    it=it_p[k];
                    q=q;
                    i++;
                    int N0=i*N;
                    int j=0;
                    for (; q < N0; q++)
                    {
                        it[q] = (it_p[k][q]);
                        it[q] = img[dc + j * width];
                        j++;
                    }
                }
            }
        }
       
    }

    void obtain_dict(double dictionary[][], double patches[][])
    {
        
         ArrayList<ArrayList< Double>> patches1=new ArrayList<ArrayList< Double>>();
         
         for(int i=0;i<patches.length;i++){
            patches1.add(new ArrayList<Double>());
             for(int j=0;j<patches[0].length;j++){
                 patches1.get(i).add(patches[i][j]);}
                         }  
         double[][] patches0=new double[patches.length][patches[0].length];
         for(int i=0;i<patches.length;i++)
             for(int j=0;j<patches[0].length;j++)
                 patches0[i][j]=patches1.get(i).get(j);
                
        //! Declarations
        int[] perm = new int[(patches.length)];

        //! Obtain random indices
        randperm(perm);

        //! Getting the initial random dictionary from patches
        int[] it_p =new int[patches.length];

        double[][] it_d = dictionary;
        for (int i = 0; i < dictionary.length; i++)
        {
            it_p[i] = perm[i];
            it_d[i] = dictionary[i];
            it_d[i] = patches0[it_p[i]];
        }

        //! Normalize column
        double[][] it = new double[dictionary.length][];
        double norm;
        for (int i = 0; i < dictionary.length; i++)
        {
            norm = 0.0;
            it[i]=dictionary[i];
            double[] it_d1;
            it_d1=it[i];
            for (int i0 = 0; i0 < it[0].length; i0++)
            {
                it_d1[i0] = it[i][i0];
                norm += (it_d1[i0]) * (it_d1[i0]);
            }
            norm = 1 / Math.sqrt(norm);
            double[] it_d2;
            it_d2=it[i];
            for (int i1 = 0; i1 < it[0].length; i1++)
            {
                it_d2[i1]=it[i][i1];
                it_d2[i1] *= norm;
            }
        }
    }

    void randperm(int[] perm)
    {
        //! Initializations
        int N = perm.length;
        int[] tmp = new int[N+1];
        tmp[1] = 1;
        for (int i = 2; i < N+1; i++)
        {
            Random rand = new Random();
            int j = rand.nextInt(i)+1;
            tmp[i] = tmp[j];
            tmp[j] = i;
        }
        int[] it_t = new int[perm.length];
        int[] it_p = perm;
        for (int i0 = 0; i0 < perm.length; i0++)
        {
            it_t[i0] = tmp[i0+1]  ;
            (it_p[i0]) = (it_t[i0]) - 1;
            perm[i0]=it_p[i0];
        }
        int h=0;
    }

    void ksvd_process( double[][] patches, double[][] dictionary, double sigma, int N1, int N2, int N_iter, double gamma, double C, boolean doReconstruction) {
        //! Declarations

        int N1_2 = N1 * N1;
        double x = 1.0 + gamma;
        double corr = (Math.sqrt(x) - 1.0) / ((double) N1_2);
        double eps = ((double) (chnls * N1_2)) * C * C * sigma * sigma;
        int h_p = patches[0].length;
        int w_p = patches.length;

        //! Mat & Vec initializations
        double[][] dict_ormp = new double[N2][h_p];
        double[][] patches_ormp = new double[w_p][h_p];
        double[][] tmp = new double[h_p][N2];
        double[] normCol = new double[N2];
        double[][] Corr = new double[h_p][h_p];
        double[] U = new double[h_p];
        double[] V;
        double[][] E = new double[w_p][h_p];

       ArrayList<ArrayList<Double>> ormp_val= new ArrayList<ArrayList<Double>>();
       ArrayList<ArrayList<Integer>> ormp_ind= new ArrayList<ArrayList<Integer>>();
       
       for(int i=0;i<w_p;i++)
       {           
           ormp_val.add(new ArrayList<Double>());
           ormp_ind.add(new ArrayList<Integer>());
       }

        //! Vector for ORMP
       
        double[][] res_ormp = new double[N2][w_p];
        ArrayList<ArrayList<Integer>> omega_table= new ArrayList<ArrayList<Integer>>();
        int[] omega_size_table = new int[N2];
        ArrayList<ArrayList<Double>> alpha= new ArrayList<ArrayList<Double>>();
          for(int i=0;i<N2;i++)
       {           
           alpha.add(new ArrayList<Double>());
           omega_table.add(new ArrayList<Integer>());
       }

        

        //! To avoid reallocation of memory
        /*  for (int k = 0; k < w_p; k++)
    {
        ormp_val[k].reserve(N2);
        ormp_ind[k].reserve(N2);
    }*/
    //
        /*for (int i=0; i < omega_table.length; i++)
    {
        it[i]=omega_table[i];
        it[i]->reserve(w_p);
    }
         */
      //  V = new double[w_p];
        //! Correcting matrix
        for (int i = 0; i < h_p; i++) 
        {
            Corr[i][i] = 1.0;
        }
        for (int c = 0; c < chnls; c++) 
        {
            double[][] it_Corr = Corr;
            for (int i = 0; i < N1_2; i++) 
            {
                double[] it=it_Corr[i];
                for (int j = 0; j < N1_2; j++) {
                    it[j] += corr;
                }
            }
        }

       
        //#pragma omp parallel for
        for (int j = 0; j < w_p; j++) 
        {
            for (int c = 0; c < chnls; c++) 
            {
                double[] it_ormp =  patches_ormp[j];
                double[] it0 = patches[j];
                for (int i = 0; i < N1_2; i++) 
                {
                    it_ormp[i] = patches_ormp[j][i] + (double) c * N1_2;
                    it0[i] = patches[j][i] + (double) c * N1_2;
                    double val = 0.0;
                    double[] it_temp = patches[j];
                    for (int k = 0; k < N1_2; k++) 
                    {
                        it_temp[k] = patches[j][k] + (double) c * N1_2;
                        val += corr * (it_temp[k]);
                    }
                    it_ormp[i] = val + it0[i];
                }
            }
        }
         
        //! Big loop
        for (int iter = 0; iter <N_iter; iter++) 
        {
            
            //! Sparse coding
            

            for (int i = 0; i < h_p; i++) 
            {
                double[] it_tmp = tmp[i];
                for (int j = 0; j < N2; j++) 
                {
                    it_tmp[j] = tmp[i][j];
                    double val = 0.0;
                    double[] it_corr_i = Corr[i];
                    double[] it_dict_j = dictionary[j];
                    for (int k = 0; k < h_p; k++) 
                    {
                        it_corr_i[k] = Corr[i][k];
                        it_dict_j[k] = dictionary[j][k];
                        val += (it_corr_i[k]) * (it_dict_j[k]);
                    }
                    (it_tmp[j]) = val * val;
                    tmp[i][j]=it_tmp[j];
                }
            }
            double[] it_normCol = normCol;
            for (int j = 0; j < N2; j++) 
            {
                it_normCol[j] = normCol[j];
                double val = 0.0;
                for (int i = 0; i < h_p; i++) {
                    val += tmp[i][j];
                }
                it_normCol[j] = 1.0 / Math.sqrt(val);
            }
            for (int i = 0; i < h_p; i++) 
            {
                double[] it_normCol_j = normCol;
                
                for (int j = 0; j < N2; j++) {
                    double[] it_corr_i = Corr[i];
                    double[] it_dict_j =dictionary[j];
                    it_normCol_j[j] = normCol[j];
                    double val = 0.0;
                    for (int k = 0; k < h_p; k++)
                    {
                        it_corr_i[k] = Corr[i][k];
                        it_dict_j[k] = dictionary[j][k];
                        val += (it_corr_i[k]) * (it_dict_j[k]);
                    }
                    dict_ormp[j][i] = val * (it_normCol_j)[j];
                }
            }

            //! ORMP process
            
           ormp_process(patches_ormp, dict_ormp, ormp_ind, ormp_val, N2, eps );
            
            
            
            
            for (int i = 0; i < w_p; i++) 
            {
                ArrayList<Integer> it_ind = ormp_ind.get(i);
                ArrayList<Double> it_val =ormp_val.get(i);
                int size = ormp_val.get(i).size();
                for (int j = 0; j < size; j++) 
                {
                    it_ind.set(j, ormp_ind.get(i).get(j));
                    it_val.set(j, ormp_val.get(i).get(j));
                    it_val.set(j, it_val.get(j)*normCol[it_ind.get(j)]) ;
                    ormp_val.get(i).set(j, it_val.get(j));
                }
            }
int h90=0;
            //! Residus
            for (int i = 0; i < N2; i++)
                for(int j=0;j<w_p;j++)
                {
                    omega_size_table[i] = 0;
                    omega_table.get(i).clear();
                    alpha.get(i).clear();
                    res_ormp[i][j]=0.0;
                }
            
             
            for (int i = 0; i < w_p; i++) 
            {
                ArrayList<Integer> it_ind = ormp_ind.get(i);
                ArrayList<Double> it_val =ormp_val.get(i);

                for (int j = 0; j < ormp_val.get(i).size() ; j++) 
                {
                    it_ind.set(j, ormp_ind.get(i).get(j));
                    it_val.set(j, ormp_val.get(i).get(j));
                    omega_table.get(it_ind.get(j)).add(i);
                    omega_size_table[it_ind.get(j)]++;
                    alpha.get(it_ind.get(j)).add(it_val.get(j));
                    res_ormp[it_ind.get(j)][i]=( it_val.get(j));
                    
                }
            }
           
           
            //! Dictionary update
            for (int l = 0; l < N2; l++) 
            {
                //! Initializations

                int omega_size = omega_size_table[l];
                double[] it_dict_l = dictionary[l];
                ArrayList<Double> it_alpha_l = alpha.get(l);
                ArrayList<Integer> it_omega_l = omega_table.get(l);

                for (int i = 0; i < U.length; i++) 
                {
                    U[i] = 0.0;
                }
               //
                if (omega_size > 0) {
                    ArrayList<Double> it_a =  it_alpha_l;
                    ArrayList<Integer> it_o = it_omega_l ;
                    for (int j = 0; j < omega_size; j++) {
                        it_a.set(j, it_alpha_l.get(j));
                        it_o.set(j, it_omega_l.get(j));
                        double[] it_d = it_dict_l;
                        double[] it_e = E[j];
                        double[] it_p = patches[it_o.get(j)];
                        for (int i = 0; i < h_p; i++) {
                            it_d[i] = it_dict_l[i];
                            it_e[i] = E[j][i];
                            it_p[i] = patches[it_o.get(j)][i];
                            (it_e[i]) = it_p[i] + (it_d[i]) * (it_a.get(j));
                        }
                    }
                    double[][] it_res = res_ormp;
                    for (int k = 0; k < N2; k++) {
                        it_res[k] = res_ormp[k];
                        ArrayList<Integer> it_o0 = it_omega_l;
                        double[] it_dict_k = dictionary[k];
                        for (int j = 0; j < omega_size; j++) {
                            it_o0.set(j, it_omega_l.get(j));
                            double val = (it_res[k])[(it_o0.get(j))];
                            if (Math.abs(val) > 0.0) {
                                double[] it_d = it_dict_k;
                                double[] it_e = E[j];
                                for (int i = 0; i < h_p; i++) {
                                    it_d[i] = it_dict_k[i];
                                    it_e[i] = E[j][i];
                                    (it_e[i]) -= (it_d[i]) * val;
                                }
                            }
                        }
                    }

                    //! SVD truncated
              
             
            
             //W(k,I) = S*V';
                    V = new double[omega_size];
                    double S = svd_trunc(E, U, V);
                    double[] um=new double[h_p] ;
                   int h00=0;
                    for(int i=0;i<h_p;i++){
                        um[i]= U[i];
                        dictionary[l][i] = um[i];
                    }
                    

                   
                    double[] it_v = V;
                    it_a = it_alpha_l;
                    it_o = it_omega_l;
                    for (int j = 0; j < omega_size; j++) {
                        it_a.set(j, it_alpha_l.get(j));
                        it_v[j] = V[j];
                        it_o.set(j, it_omega_l.get(j));
                        it_a.set(j,it_v[j]* S);
                        res_ormp[l][it_o.get(j)] =  it_a.get(j);
                        

                    }
                }
            }
            int h1=0;
        }
        int h1=0;


        if (doReconstruction) {
            //! Final patches estimation
            for (int i = 0; i < patches.length; i++) 
            {
                for (int j = 0; j < patches[0].length; j++) 
                {
                    patches[i][j] = 0.00;
                }
            }
//		#pragma omp parallel for
            for (int l = 0; l < N2; l++) {
                ArrayList<Double> it_a = alpha.get(l);
                ArrayList<Integer> it_o = omega_table.get(l);
                int omega_size = omega_size_table[l];
                for (int j = 0; j < omega_size; j++) {
                    it_a.set(j, alpha.get(l).get(j));
                    it_o.set(j, omega_table.get(l).get(j));
                    double val = (it_a.get(j));
                    double[] it_d = dictionary[l];
                    double[] it_p = patches[it_o.get(j)];
                    for (int k = 0; k < h_p; k++) {
                        it_d[k] = dictionary[l][k];
                        it_p[k] = patches[it_o.get(j)][k];
                        (it_p[k]) += (it_d[k]) * val;
                        
                    }
                }
            }

            //! First, obtention of the denoised image without weighting with lambda
            patches2im(patches, img_denoised, img_noisy, width, height, chnls, 0.0, N1);

            //! Second, obtention of lambda from the norm between img_denoised and img_noisy
            double d = 0.0;
            for (int k = 0; k < width * height * chnls; k++) 
            {
                d += (img_denoised[k] - img_noisy[k]) * (img_denoised[k] - img_noisy[k]);
            }
            d /= (height * width * chnls * sigma * sigma);
            double lambda;
            if(fSigma==0)
            {
               lambda=0.05;
            }
            else
            {
            lambda = Math.abs(Math.sqrt(d) - 1.0);
            }
            //! Finally, obtention of the denoised image with lambda
            patches2im( patches, img_denoised, img_noisy, width, height, chnls, lambda, N1);
        }
    }

    void patches2im(  double[][] patches, float[] img, double[] img_ref, int width, int height, int chnls, double lambda, int N) {
        //! Declarations
        int size=height * width;
        int h_p = height - N + 1;
        int wh_p = (width - N + 1) * h_p;

        for (int c = 0; c < chnls; c++) 
        {
            double[] denominator = new double[size];
            double[] numerator = new double[size];
            double[][] it_p = patches;


            //! Aggregation
            for (int k = 0; k < wh_p; k++)
            {
                it_p[k] = patches[k];
                int ind = (k % h_p) * width + k / h_p;
                double[] it = it_p[k];
                int i=0;
                int p=0;
                for (int q = 0; q < N; q++, ind++) 
                {
                    it=it_p[k];
                    p=p;
                    i++;
                    int j=0;
                    int N0=i*N;
                    for (; p < N0; p++) 
                    {
                        it[p] = it_p[k][p] + c * N * N;
                        numerator[j * width + ind] += it[p];
                        denominator[j * width + ind]++;
                        j++;
                    }
                }
            }

            //! Weighting
            int dc_i = 0;
            double[] it_d = denominator;
            double[] it_n = numerator;
            
            for (int k = 0; k < height * width; k++, dc_i++) 
            {
                it_d[k] = denominator[k];
                it_n[k] = numerator[k];
                img[k] = (float)( ((it_n[k]) + lambda * (img_ref[dc_i])) / ((it_d[k]) + lambda));
            
            
            }
            

        }
       
       
    }
    

    void ormp_process(double[][] X, double[][] D, ArrayList<ArrayList<Integer>> ind_v, ArrayList<ArrayList<Double>> val_v, int L, double eps) {
        //! Declarations
        int n = X[0].length;
        int Np = X.length;
        int k = D.length;

        //! Initializations
        if (L <= 0) {
            return;
        }

        L = Math.min(n, Math.min(L, k));
        //   #pragma omp parallel shared(ind_v, val_v, X, D)
        {
            //! Declarations
            double[] norm = new double[k], scores = new double[k], x_T = new double[k];
            double[][] A = new double[L][L], D_ELj = new double[k][L], D_D = new double[k][k], D_DLj = new double[k][L];

            //! Compute the scalar products between the atoms of the dictionary
            for (int j = 0; j < k; j++) {
                double[] it_D_D = D_D[j];
                for (int i = 0; i < k; i++) {
                    it_D_D[i] = D_D[j][i];
                    double[] D_i = D[i];
                    double[] D_j = D[j];
                    double val = 0;
                    for (int s = 0; s < n; s++) {
                        D_i[s] = D[i][s];
                        D_j[s] = D[j][s];
                        val += (double) (D_i[s]) * (double) (D_j[s]);
                    }
                    (it_D_D[i]) = (double) val;
                    
                }
            }


            //       #pragma omp for schedule(dynamic) nowait
            for (int i = 0; i < Np; ++i) 
            {
                //! Initialization
                ind_v.get(i).clear();
                val_v.get(i).clear();
                double[] X_i = X[i];

                //! Compute the norm of X[i]
                double normX = 0.0;
                double[] it_x = X_i;

                for (int j = 0; j < X[i].length; j++) 
                {
                    it_x[j] = X[i][j];
                    normX += (it_x[j]) * (it_x[j]);
                }

                //! Compute the scalar products between X[i] and the elements
                //! of the dictionary
                for (int j = 0; j < k; j++) 
                {
                    double[] it_d = D[j];
                    double[] it_x0 = X_i;
                    double val = 0.0;
                    for (int s = 0; s < n; s++) {
                        it_d[s] = D[j][s];
                        it_x0[s] = X_i[s];
                        val += (it_d[s]) * (it_x0[s]);
                    }
                    x_T[j] = (double) val;
                }

                coreORMP(D, D_D, scores, norm, A, D_ELj, D_DLj, x_T, ind_v.get(i), val_v.get(i), eps, (double) normX);
                
            } //! End of I-lop

        } //! End of parallel section
    }

    void coreORMP(double[][] D, double[][] D_D, double[] scores, double[] norm, double[][] A, double[][] D_ELj, double[][] D_DLj, double[] x_T, ArrayList<Integer> ind, ArrayList<Double> coord, double eps, double normr) {
        //! Declarations
        int L = A.length;
        int p = D.length;

        if (normr <= eps || L == 0) {
            return;
        }

        //! Initializations
        for(int i=0;i<x_T.length;i++)
        scores[i] = x_T[i];
        for (int i = 0; i < norm.length; i++) 
            norm[i]=1.0;
        
        for(int i=0;i<L;i++)
            for(int j=0;j<L;j++)
                A[i][j]=0;
        
        //! Declarations
        double[] A_j, it_A_j, it_D_ELj_lj, it_gs, x_T_i, norm_i, scores_i, it_A_j_tmp=new double[L];
        double[] D_lj, it_D_DLj;
        double[][] it_D_ELj, it_A=new double[L][L];
        ArrayList<Double> x_el=new ArrayList<Double>();
        
        
        
        //! Loop over j
        int j;
        for (j = 0; j < L; j++) 
        {
            //! Initialization
            int lj = ind_fmax(scores);
            A_j = A[j];

            //! Stop if we cannot inverse norm[lj]
            if (norm[lj] < 1e-6) 
            {
                break;
            }

            double invNorm = 1.0 / (double) Math.sqrt(norm[lj]);
            double x_elj = (double) x_T[lj] * invNorm;
            double delta = x_elj * x_elj;
            
            
           
            x_el.add(x_elj);
            coord.add(x_elj); //! The coordinate of x on the last chosen vector
            normr = (double) ((double) normr - delta);//! Update of the residual
            ind.add(lj); //! Memorize the chosen index
            

            //! Gram-Schmidt Algorithm, Update of Aj
            it_A_j = A_j;
            it_D_ELj_lj = D_ELj[lj];
            for (int i = 0; i < j; i++) 
            {
                it_A_j[i] = A_j[i];
                it_D_ELj_lj[i] = D_ELj[lj][i];
                (it_A_j[i]) = (it_D_ELj_lj[i]);
                
            }
           int  g=0;
            for (int i = 0; i < j; i++,g++) 
            {
                it_A_j[g] = A_j[i];
                double sum = 0.0;
                int g0=0;

                for (int s = 0; s < j - i; s++,g0++) {
                    it_A[s] = A[s+i];
                    it_A_j_tmp[g0] = it_A_j[g+g0];
                    sum -= (double) ((it_A[s])[i]) * (double) (it_A_j_tmp[g0]);
                }
                (it_A_j[g]) = (double) (sum * invNorm);
            }
            (it_A_j[g]) = (double) invNorm;

            if (j == L - 1 || (normr <= eps)) {
                j++;
                break;
            }

            //! Update of D_
            it_D_DLj = D_D[lj];
            for (int i = 0; i < p; i++) {
                it_D_DLj[i] = D_D[lj][i];
                D_DLj[i][j] = (it_D_DLj[i]);
            }
            
            //! Compute the scalar product D[j]_D[lj] and memorize it now
            double val = 0.0;
            D_lj = D[lj];
            double[] D_j = D[j];
            for (int i = 0; i < D[j].length; i++) {
                D_j[i] = D[j][i];
                D_lj[i] = D[lj][i];
                val += (double) (D_j[i]) * (double) (D_lj[i]);
            }
            D_DLj[j][j] = (double) val;

            //! Update of D_ELj, x_T, norm, and scores.
            x_T_i = x_T;
            norm_i = norm;
            scores_i = scores;
            it_D_ELj = D_ELj;
            for (int i = 0; i < p; i++) {
                x_T_i[i] = x_T[i];
                norm_i[i] = norm[i];
                scores_i[i] = scores[i];
                it_D_ELj[i] = D_ELj[i];
                double val0 = 0.0;
                double [] it_A_j0=A_j;
                it_gs = D_DLj[i];
                for (int s = 0; s < j + 1; s++) {
                   it_A_j0[s] = A_j[s]+1;
                  //it_A_j0[s] = A_j[s];
                    it_gs[s] = D_DLj[i][s];
                    val0 += (double) (it_A_j0[s]) * (double) (it_gs[s]);
                }
                it_D_ELj[i][j] = val0;
                D_ELj[i][j]=it_D_ELj[i][j];
                (x_T_i[i]) =  ((x_T_i[i]) - x_elj * val0);
                x_T[i]=x_T_i[i];
                (norm_i[i]) =  ((double) (norm_i[i]) - val0 * val0); 
                norm[i]=norm_i[i];
                (scores_i[i]) =  ((double) (x_T_i[i]) * (double) (x_T_i[i]) / (double) (norm_i)[i]);
                scores[i]=scores_i[i];
            }

        }
        
        //! Compute the final coordinates of x on the chosen atoms of the dictionary
        ArrayList<Double> it_coord = coord;
        for (int i = 0; i < j; i++) {
            it_coord.set(i, coord.get(i));
            double sum = 0;
            double[][] it_a = new double[L][L];
            ArrayList<Double> it_x = x_el;
            for (int s = 0; s < j - i; s++) {
                it_a[s] = A[s+i] ;
                it_x.set(s, x_el.get(s+i));
                sum += (double) ((it_a[s])[i]) * (double) (it_x.get(s));
            }
            it_coord.set(i, sum);
        }
    }

    int ind_fmax(double[] V) {
        double val = 0.0f;
        int ind = 0;
        double[] it_v = V;
        int j = 0;

        //! Find the index of the maximum of V
        for (int i = 0; i < V.length; i++, j++) 
        {
            it_v[i] = V[i];
            if (val < Math.abs(it_v[i])) 
            {
                val = Math.abs(it_v[i]);
                ind = j;
            }
        }

        return ind;
    }

    double svd_trunc(double[][] tX, double[] U, double[] V) {
        //! Declarations
        double epsilon = 10e-6;
        int max_iter = 100;
        int iter = 0;
        boolean go_on = true;
        double S_old = 0;
        double S = 0;
        int m = U.length;
        int n = V.length;

        double norm = 0.0;
        double[] it_v = V;
        for (int j = 0; j < n; j++) {
            it_v[j] = V[j];
            double val = 0.0;
            double[] it_x = tX[j];
            for (int i = 0; i < m; i++) {
                it_x[i] = tX[j][i];
                val += (double) Math.abs(it_x[i]);
            }
            (it_v[j]) = (double) val;
            norm += val * val;
        }

        double s_inv = -1.0 / Math.sqrt(norm);

        for (int j = 0; j < n; j++) {
            (it_v[j]) *= (double) s_inv;
        }

        while (iter < max_iter && go_on) {
            S_old = S;
            double[] it_u = U;
            norm = 0.0;
            for (int i = 0; i < m; i++) {
                it_u[i] = U[i];
                double value = 0.0;
                it_v = V;
                for (int j = 0; j < n; j++) {
                    it_v[j] = V[j];
                    value += (double) tX[j][i] * (double) (it_v[j]);
                }
                (it_u[i]) = (double) value;
                norm += value * value;
            }
            s_inv = 1.0 / Math.sqrt(norm);

            for (int i = 0; i < m; i++) {
                it_u[i] = U[i];
                (it_u[i]) *= s_inv;
            }
            for (int i = 0; i < n; i++) {
                it_v[i] = V[i];
                (it_v[i]) = 0.0;
            }
            for (int j = 0; j < n; j++) {
                it_v[j] = V[j];
                double[] it_x = tX[j];
                for (int i = 0; i < m; i++) {
                    it_u[i] = U[i];
                    it_x[i] = tX[j][i];
                    (it_v[j]) += (it_x[i]) * (it_u[i]);
                }
            }

            norm = 0.0;
            for (int i = 0; i < n; i++) {
                it_v[i] = V[i];
                norm += (it_v[i]) * (it_v[i]);
            }
            S = Math.sqrt(norm);
            s_inv = 1.0 / S;

            for (int i = 0; i < n; i++) {
                it_v[i] = V[i];
                (it_v[i]) *= s_inv;
            }
            iter++;
            go_on = Math.abs(S - S_old) > epsilon * S;
        }

        return S;
    }

    @Override
    public boolean checkInput() {
     return true;//   throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public Img<T> getResult() {
       return output;
    }


}
