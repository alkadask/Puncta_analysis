# Puncta_analysis

SAMPLE PREP:
Grow worms at 20 degrees
Mount young adults on groove agar pads in M9
Use 5 mM Levamisole in M9 for immobilization.

DATA COLLECTION:
Microscope system: Keyence
Objective: 40x
Digital zoom: 1x
Exposure time: 2.5 sec
Mode: High resolution (No binning)
Only image the neuron closest to the cover slip.
Image stitching in uncompressed mode

FILE NAMING:
Group folder name: yyyymmdd_Strain_xx
File prefix: yyyymmdd_Strain
Stitched file name: yyyymmdd_Strain_xx_Neuron(s)

IMAGE PRE-PROCESSING:
Open image in Fiji.
Stitched files are rotated to orient anterior end towards left unless there are multiple worms in the frame.
Shortcut to rotate 90 degrees right: L
To determine L/R: Look for alae position relative to TRN. If alae is above TRN then ALML, otherwise ALMR. In case of ALMR, often the AVM cell body is also visible.
Trace neuron using segmented line tool. Use spline fit. Always trace starting from cell body moving towards distal end.
Save ROIs.
Press spacebar to switch to hand mode for easy scrolling through zoomed image.
Before straightening make sure to reset brightness and contrast levels.
Straighten shortcut: q
Line width for straightening: 20 px
Save straightened image.
Straightened file name: yyyymmdd_Strain_xx_Neuron-x
Every straightened image should be easily traceable to its raw image file.

Note: The ImageJ macro script, Neuron_tracing.ijm can partially automate the repetitive steps. Manual verification and adjustments are still needed but the workload is considerably reduced.

QUALITY CONTROL:
Do not straighten images with the following issues (trace neuron and save ROI):
-- Puncta out of focus
-- Movement of worm during acquisition
-- Improper stitching

PUNCTA ANALYSIS:
Open script Peakfinder25.py
Change the file path to the folder containing straightened images
Run script
Algorithm should output (display only, not saved) traces of each neuron, with puncta identified as green dots.
Verify visually if peaks are appropriately identified
Algorithm will save an excel file with all raw data, peaks data, and overall analysis data in the output folder. The algorithm will save to a new file with timestamp every time the program is run (Files will not be overwritten)
For all further analysis and plotting, read this excel file in pandas and extract data as needed.

Changes in Peakfinder26.py compared to Peakfinder25.py:
The most major change is that in this version we estimate the baseline fluorescence between puncta by fitting a cubic spline through points sampled by a argrelmin function with high order. We then subtract the baseline from the total neurite fluorescence. This operation brings the base of most major peaks down to zero. The find_peaks function is applied on this array.
Minor differences - No smoothing is applied and also not using the controversial 'signal to noise ratio' calculation used in the last version for normalization. This affords direct comparison of fluorescence intensity levels between strains.
