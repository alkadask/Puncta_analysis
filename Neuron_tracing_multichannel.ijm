/* Before using this macro, make sure of the following points:
 * Each image has only one neuron of interest. If image has more then crop out each neuron containing regions and save as separate files.
 * In each image the anterior end of the worm faces RIGHT. NOTE: This is opposite of general C. elegans image display convention but this is better for image processing.
 */
keepgoing = "yes"
Dialog.create("Which channel has neuron label?");
Dialog.addChoice("Channel selection", newArray("red", "green"));
Dialog.show()
neuronlabel = Dialog.getChoice();
if (neuronlabel=="red"){
	punctalabel="green";
	}
else{
	punctalabel="red";
}
while (keepgoing == "yes"){
	run("Close All");
	run("Clear Results");

	f= File.openDialog("Select file");
	fpath = File.getParent(f);
	name = File.getName(f);
	parentdir = substring(fpath, 0,lengthOf(fpath)-10);
	open(name);
	selectWindow(name);
	imheight=getHeight();		//get image height, this will be needed later to convert y-cordinate values from plot to array indexes
	run("Split Channels");
	selectWindow(name+" (blue)");
	close();
	selectWindow(name+" ("+neuronlabel+")");
	run("Duplicate...", " ");	//Duplicate image
	run("Subtract Background...", "rolling=50"); //subtract background using default settings
	run("Auto Local Threshold", "method=Bernsen radius=15 parameter_1=0 parameter_2=0 white"); //auto local threshold using default settings
	run("Invert");

	//Find selection area of neuron and clear all other thresholded points
	run("Create Selection");
	roiManager("Add");
	roiManager("Show All with labels");
	roiManager("Select", 0);
	run("Enlarge...", "enlarge=20");	//enlarge selected areas by a certain number of pixel radius. this operation effectively joins the puncta in the neuron to get a single area selection spanning the entire neuron
	roiManager("Add");
	roiManager("Select", 1);
	if (selectionType()==9){
		roiManager("Split");
		run("Select None");
		roiManager("Deselect");
		
		Dialog.create("Check the ROIs to see if the entire neuron is covered in one selection area");
		Dialog.addChoice("ROI selection", newArray("automatic", "manual"));
		Dialog.show()
		choice = Dialog.getChoice();
	
		//Automatic ROI selection
		if (choice=="automatic"){
			run("Set Measurements...", "perimeter bounding shape feret's redirect=None decimal=3");
			roiManager("Measure");
			//find the index of the roi with highest perimeter
			maxPerimindex=2;
			for (j = 3; j < roiManager("count"); j++) {
				if(getResult("Perim.",j)>getResult("Perim.",maxPerimindex)){
					maxPerimindex=j;
				}
			}
			roiManager("Select", maxPerimindex);
		}

		//Manual ROI selection
		else{
			waitForUser("Select ROIs representing the neuron. Click OK when done");
			roiManager("Combine");
			roiManager("Deselect");
			roiManager("Deselect");
			roiManager("Delete");
			roiManager("Add");
			roiManager("Select", 0);
		}
	}
	
	//clear all threshold regions outside of the selected roi
	setBackgroundColor(255, 255, 255);
	run("Clear Outside");
	run("Select None");
	mask = replace(name, ".tif", "_mask.tif");
	saveAs("Tiff", parentdir+"/masks and rois/"+mask);
	
	roiManager("Delete");
	
	//Fit a line along the neuron and get its coordinates
	run("Analyze Line Graph");
	selectWindow("Line Graph");
	Plot.getValues(xpoints, ypoints);
	for (k=1; k<lengthOf(ypoints); k++){
		ypoints[k]=imheight-ypoints[k];
	}
	//Downsample points by a factor of 20
	size=xpoints.length;
	xp=newArray(size/20-1);
	yp=newArray(size/20-1);
	for (k=1;k<xp.length;k++){
		xp[k]=xpoints[k*20];
		yp[k]=ypoints[k*20];
	}
	close();
	close();
	
	//add the spline selection to the neuron label channel
	selectWindow(name+" ("+neuronlabel+")");
	Roi.setPolylineSplineAnchors(xp, yp)
	roiManager("Add");
	roiManager("Select", 0);
	waitForUser("Visually cross check to make sure the line correctly traces the entire neuron. Trim the proximal end to remove cell body and extend the distal end if necessary. Click OK when done");
	roiManager("update");
	roiManager("save selected", parentdir+"/masks and rois/"+name+".roi");
	
	//Add selection line to puncta channel and straighten line
	selectWindow(name+" ("+punctalabel+")");
	roiManager("Select", 0);
	waitForUser("Visually cross check to make sure the line correctly traces the entire neuron. Trim the proximal end to remove cell body and extend the distal end if necessary. Click OK when done");
	run("Straighten...", "line=20");
	straightened = replace(name, ".tif", "_straightened.tif");
	saveAs("Tiff", parentdir+"/straightened/"+straightened);
	
	roiManager("Delete");
	run("Close All");
	run("Clear Results");

	Dialog.create("Do you want to process another image?");
	Dialog.addChoice("Keep going?", newArray("yes", "no"));
	Dialog.show()
	keepgoing = Dialog.getChoice();
}