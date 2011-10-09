package fiji.plugin.trackmate.tracking;

import fiji.plugin.trackmate.gui.SimpleLAPTrackerSettingsPanel;
import fiji.plugin.trackmate.gui.TrackerSettingsPanel;

public class SimpleLAPTracker extends LAPTracker {

	@Override
	public String toString() {
		return "Simple LAP tracker";
	}

	@Override
	public String getInfoText() {
		return "<html>" +
				"This tracker is identical to the LAP tracker present in this plugin, except that it <br>" +
				"proposes fewer tuning options. Namely, only gap closing is allowed, based solely on <br>" +
				"a distance and time condition. Track splitting and merging are not allowed, resulting <br>" +
				"in having non-branching tracks." +
				" </html>";
	}
	
	@Override
	public TrackerSettings createDefaultSettings() {
		TrackerSettings ts = new TrackerSettings() {
			public TrackerSettingsPanel createConfigurationPanel() {
				return new SimpleLAPTrackerSettingsPanel();
			}
		};
		ts.allowMerging = false;
		ts.allowSplitting = false;
		return ts;
	}

}