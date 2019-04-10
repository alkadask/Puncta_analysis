# Puncta_analysis
Code for analyzing puncta distribution of labeled proteins in C. elegans touch receptor neurons

SAMPLE PREP:
Grow worms at 20 degrees
Mount young adults on groove agar pads in 150 mOSm Buffer
Use 5 mM Levamisole in M9 for immobilization.

DATA COLLECTION:
Microscope system: Keyence BZ-X800
Objective: 40x oil
Digital zoom: 1x
Exposure time: 2.5 sec
Mode: High resolution
Only image the neuron closest to the cover slip.
Stitching in uncompressed mode

FILE NAMING:
Group folder name: yyyymmdd_Strain_xx
File prefix: yyyymmdd_Strain
Stitched file name: yyyymmdd_Strain_xx_Neuron(s)
Stitched files are rotated to orient anterior end towards left unless there are multiple worms in the frame.
To determine L/R: Look for alae position relative to TRN. If alae is above TRN then ALML, otherwise ALMR.

IMAGE PRE-PROCESSING:
Open image in Fiji.
Trace neuron using segmented line tool. Use spline fit. Always trace starting from cell body moving towards distal end.
Save ROIs.
Press spacebar to switch to hand mode for easy scrolling through zoomed image while tracing.
Before straightening make sure to reset brightness and contrast levels to original.
For images with puncta out of focus, movement during image aquisition, improper stitching or other quality issues, trace neuron and save ROI but do not straighten.
Line width for straightening: 20 px
Save straightened image.
Straightened file name: yyyymmdd_Strain_xx_Neuron-x
Every straightened image should be easily traceable to its raw image file.
