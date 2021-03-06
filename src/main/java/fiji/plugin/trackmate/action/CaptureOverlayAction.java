package fiji.plugin.trackmate.action;

import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.io.File;

import javax.swing.ImageIcon;
import javax.swing.JOptionPane;

import org.scijava.plugin.Plugin;

import fiji.plugin.trackmate.LoadTrackMatePlugIn_;
import fiji.plugin.trackmate.TrackMate;
import fiji.plugin.trackmate.gui.TrackMateGUIController;
import fiji.plugin.trackmate.gui.TrackMateWizard;
import fiji.plugin.trackmate.visualization.TrackMateModelView;
import fiji.plugin.trackmate.visualization.trackscheme.TrackSchemeFrame;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.Calibration;
import ij.process.ColorProcessor;

public class CaptureOverlayAction extends AbstractTMAction
{

	public static final ImageIcon ICON = new ImageIcon( TrackSchemeFrame.class.getResource( "resources/camera_go.png" ) );

	public static final String NAME = "Capture overlay";

	public static final String KEY = "CAPTURE_OVERLAY";

	public static final String INFO_TEXT = "<html>" +
			"If the current displayer is the HyperstackDisplayer, this action <br>" +
			"will capture the TrackMate overlay with current display settings. <br>" +
			"That is: a new RGB stack will be created (careful with large data) where <br>" +
			"each frame contains a RGB snapshot of the TrackMate display. " +
			"<p>" +
			"It can take long since we pause between each frame to ensure the whole <br>" +
			"overlay is redrawn. The current zoom is taken into account. <br>" +
			"Also, make sure nothing is moved over the image while capturing. " +
			"</html>";

	private final TrackMateWizard gui;

	private static int firstFrame = -1;

	private static int lastFrame = -1;

	private ImagePlus capture;

	public CaptureOverlayAction( final TrackMateWizard gui )
	{
		this.gui = gui;
	}

	@Override
	public void execute( final TrackMate trackmate )
	{
		final ImagePlus imp = trackmate.getSettings().imp;

		if ( firstFrame < 0 )
			firstFrame = 1;
		firstFrame = Math.max( firstFrame, 1 );
		if ( lastFrame < 0 )
			lastFrame = imp.getNFrames();
		lastFrame = Math.min( lastFrame, imp.getNFrames() );

		final CaptureOverlayPanel panel = new CaptureOverlayPanel( firstFrame, lastFrame );
		final int userInput = JOptionPane.showConfirmDialog( gui, panel, "Capture TrackMate overlay", JOptionPane.OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, TrackMateWizard.TRACKMATE_ICON );
		if ( userInput != JOptionPane.OK_OPTION )
			return;

		final int first = panel.getFirstFrame();
		final int last = panel.getLastFrame();
		firstFrame = Math.min( last, first );
		lastFrame = Math.max( last, first );
		firstFrame = Math.max( 1, firstFrame );
		lastFrame = Math.min( imp.getNFrames(), lastFrame );

		logger.log( "Capturing TrackMate overlay from frame " + firstFrame + " to " + lastFrame + ".\n" );
		final Rectangle bounds = imp.getCanvas().getBounds();
		final int width = bounds.width;
		final int height = bounds.height;
		final ImageStack stack = new ImageStack( width, height );
		logger.log( " done.\n" );

		final int nCaptures = lastFrame - firstFrame + 1;
		logger.log( "Performing capture..." );
		final int channel = imp.getChannel();
		final int slice = imp.getSlice();
		imp.getCanvas().hideZoomIndicator( true );
		for ( int frame = firstFrame; frame <= lastFrame; frame++ )
		{
			logger.setProgress( ( float ) ( frame - firstFrame ) / nCaptures );
			imp.setPositionWithoutUpdate( channel, slice, frame + 1 );
			final BufferedImage bi = new BufferedImage( width, height, BufferedImage.TYPE_INT_ARGB );
			imp.getCanvas().paint( bi.getGraphics() );
			final ColorProcessor cp = new ColorProcessor( bi );
			final int index = imp.getStackIndex( channel, slice, frame );
			stack.addSlice( imp.getImageStack().getSliceLabel( index ), cp );
		}
		imp.getCanvas().hideZoomIndicator( false );
		this.capture = new ImagePlus( "TrackMate capture of " + imp.getShortTitle(), stack );
		transferCalibration( imp, capture );
		capture.show();

		logger.log( " done.\n" );
		logger.setProgress( 0. );
	}

	/**
	 * Transfers the calibration of an {@link ImagePlus} to another one,
	 * generated from a capture of the first one. Pixels sizes are adapter
	 * depending on the zoom level during capture.
	 * 
	 * @param from
	 *            the imp to copy from.
	 * @param to
	 *            the imp to copy to.
	 */
	private void transferCalibration( final ImagePlus from, final ImagePlus to )
	{
		final Calibration fc = from.getCalibration();
		final Calibration tc = to.getCalibration();

		tc.setUnit( fc.getUnit() );
		tc.setTimeUnit( fc.getTimeUnit() );
		tc.frameInterval = fc.frameInterval;

		final double mag = from.getCanvas().getMagnification();
		tc.pixelWidth = fc.pixelWidth / mag;
		tc.pixelHeight = fc.pixelHeight / mag;
		tc.pixelDepth = fc.pixelDepth;
	}

	/**
	 * Returns a the {@link ImagePlus} resulting from the last
	 * {@link #execute(TrackMate)} call.
	 * 
	 * @return a RGB {@link ImagePlus}, or <code>null</code> if the
	 *         {@link #execute(TrackMate)} has not been called yet.
	 */
	public ImagePlus getCapture()
	{
		return capture;
	}

	@Plugin( type = TrackMateActionFactory.class )
	public static class Factory implements TrackMateActionFactory
	{

		@Override
		public String getInfoText()
		{
			return INFO_TEXT;
		}

		@Override
		public String getKey()
		{
			return KEY;
		}

		@Override
		public TrackMateAction create( final TrackMateGUIController controller )
		{
			return new CaptureOverlayAction( controller.getGUI() );
		}

		@Override
		public ImageIcon getIcon()
		{
			return ICON;
		}

		@Override
		public String getName()
		{
			return NAME;
		}

	}

	public static void main( final String[] args )
	{
		ImageJ.main( args );
		final File file = new File( "/Users/tinevez/Google Drive/Projects/Contacts/raw data/2015-09-17/Trackmate files/SiC + SAg2_1_20_BCells.xml" );
		final LoadTrackMatePlugIn_ loader = new LoadTrackMatePlugIn_();
		loader.run( file.getAbsolutePath() );

		loader.getSettings().imp.setDisplayMode( IJ.GRAYSCALE );
		loader.getSettings().imp.getCanvas().zoomIn( 50, 50 );
		loader.getSettings().imp.getCanvas().zoomIn( 50, 50 );
		loader.getSettings().imp.getCanvas().zoomIn( 50, 50 );
		loader.getSettings().imp.getCanvas().zoomIn( 50, 50 );
		loader.getSettings().imp.getCanvas().zoomIn( 50, 50 );

		for ( final TrackMateModelView view : loader.getController().getGuimodel().getViews() )
		{
			view.setDisplaySettings( TrackMateModelView.KEY_TRACK_DISPLAY_DEPTH, 100 );
			view.setDisplaySettings( TrackMateModelView.KEY_TRACK_DISPLAY_MODE, TrackMateModelView.TRACK_DISPLAY_MODE_LOCAL_BACKWARD );
		}

		final TrackMate trackmate = loader.getController().getPlugin();
		new CaptureOverlayAction( loader.getController().getGUI() ).execute( trackmate );
	}
}
