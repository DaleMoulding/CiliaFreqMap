// Cilia beat analysis via ImageJ array.fourier 
// Takes a time series image
// extracts movement by making a rolling static projection over time
// subtracting this from the original time series
// as a 32 bit result to give oscillations around the mean
// as a sine wave with mean value near to zero
// Reslices to place the z axis along the y axis
// allows fast data extraction to arrays
// Uses Fourier array analysis to get the beating frequency of each pixel

// Macro under development by Dale Moulding.
// MIT License

//Copyright (c) 2025 DaleMoulding

//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:

//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.

//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.

//***********changelog*******************
//v001 for single files. Add a folder loop
//v002 change from inoput image durastion to fps /Hz to allow adjustment of duration when trimming movies to 128 512 1024 etc frames.
//v003 try and find way reducing image duration is giving double the expected beatFreq
// the problem seems to be repeatable in the test fourier array. use 1029 then 1028 array length. ? reset the fourier array with a small data set?
// solution: cropping the images to 1-256 is actually 257 frames!
//v004a why is a different peak strength giving a different result and not just more zeros?
// 004a fixed, needed to reset arrays within x loop
// v005 filter the outputs to include only frequencies within a pre-defined range
// ?? i will do this on the output results, removing pixels with values outside the defined range.
// could also do this on the maxima of the fourier array? So if the first peak is outside the range, use the next one, until on is on range?
// ?? or filter the fourier array before finding maxima? use array.slice to remove first x values and last x values.
// or better still, just work throught the maxima. If 1st maxima not in range, go to next etc etc. Until find one in range OR run out.
// v006 add measurments. Frequency. % image moving. Freq distrubution.
// v101 Remove adaptive?
// v101 measure the StDev of the fourier array excluing the top ? 5%, 10% user defined. Use this StdDev multiplied by a user defined noise tolerance 
// to set the maxima tolerance to find the primary frequency.
// v101 or use the mean
// v102 set the images to 2^n frames. 256 not 255. 
// v103 add print out start. Input frequency, smoothing options, fourier results noise tolerance
// v103 when doing exactly 2^n pixels, the peak frequency is 2x the frequency of beating! Therefore divide by 2.
// v103 save a list in excel of the histogram plot for frequency distribution
// v104 had to exit batch mode to do the histogram. ? do via getimage stats?
// v104 tidy the results table
// v105 add a time stamp to results file and log
// v106 fixed bug in timestring2 if time was before 10am.
// v107 change histograms to 400 bins.
// v108 fix crash when there are no results. Line 269 clear outside crashes as there is no selection.
// v109n increase input speed up to 500 fps & Average frames up to 100
// v110 reslice the output movie back from a kymograph
// v111 try some filters on the processed moving image before the fourier analysis
// v112 fix bug line 283 that an image with a few moving pixels that then get median filtered to zero don't crash when there is nothing to clear as there is no selection
// v113 reduce min image length to 128 Frames
// v113 make macro simply skip files that are too short. & record to the log
// v113 fresh start or close all to help clear memory?
// v113 rename result table columns
// v113 tidy up log file
// v114 swap collect garbage & close all lines to close before clearing RAM
// v114 further log file clean up & timings added

// v201 make a power spectrum. Save every fourier array as a new slice of an image. Then make an average z project of the image and plot that.
// v202 power spectrum works, make a filtered one, just from pixels that pass the background/noise test
// v202 use just the filtered power spectrum
// v202 seems not to work? Power spectrum as plotted doesn't fit well. Abandon this.
// v202 instead make a nicer plot from the peak frequency distribution.
// v203 troubleshoot v202
// v203 solution: the height of the resliced image is 1 more than the line number (y) of the last line. i.e. a 512 line image has line number 0-511! So use 'height-1' to get the profile!
// v204 incorporate the solution, height-1 for the profile for the fourier analysis, and then no longer dived by 2 to get the frequency
// v204 the power spectrum should now be fine also.
// v204 also fix nomalise image line ~201
// v205 remove background freq % selection, just set to 80%
// v205 add an option to normalise the image between frames (to remove low freq oscialltions due to image capture artefacts / fluctuating ination)
// instead of measuring the image variance for the normalistion step if loop.
// v206 Array sizes for BeatFreq & Confidence increased from width to width*slices so they arte opened at teh correct size = total number of pixels
// v206 change image name to name without extension
// v207 make a grpah from the power spectrum, converted to Hz
// v207 make the histogram bins equal the number of fourier measurments
// v207b AVI files through bioformats open as 3 channels, squish them back to one
// v207c on my PC only clear outside images were giving a value of 1 not zero, add:
// setForegroundColor(255, 255, 255);
// setBackgroundColor(0, 0, 0);
// v208 set min colour balaqnce to zero rather than the lowest freq cut off measured, so low freqs are visible on the beatfreq image.
// v209 save the results table each image? And overright it each time?
// v209 !! Limit to threshold gives the area correctly without having to measure twice!
//v209 change output from BeatFreq to FreqMap
// v210 change macro name, make it installable, include MIT license
// v211 option to save or not the static movie
// v211. add min and max Hz that can be measured from the image, and the Hz resolution
// v211 set the scale of the power spectrum tif to the freq resolution
// v211 histogram output bins now equal the Hz resolution
// v future... add an image of Beat Freq, Confidence, both with scale bars and % image moving and Means Hz. As a jpg.

macro "CiliaFreqMap"{

run("Options...", "iterations=1 count=1 black"); // black backgrounds
setForegroundColor(255, 255, 255);
setBackgroundColor(0, 0, 0); // belt and braces black backgrounds v208B for "Clear Outside"
run("Set Measurements...", "mean standard modal min median area_fraction limit display redirect=None decimal=3"); //v209 limit to threshold
run("Clear Results"); //v104
if (isOpen("Log")){ //v105
	selectWindow("Log");
	run("Close");
}

windowType="Hann"; //None, Hamming, Hann or Flattop for fourier array

#@File(label = "Input directory", style = "directory") input
#@File(label = "Output directory", style = "directory") output

// get image parameters for pre-processing
#@ Float InputFPS (label = "Capture speed of movies (FPS / Hz)", style="slider,format:0.00", value=87.07, min=10, max=500, stepSize=0.01, persist=true)
#@ Integer 	LowerHzCutOff	(label = "Lowest frequency (Hz) to measure:", value=3, min=0, max=100, persist=true, style="slider")
#@ Integer 	UpperHzCutOff	(label = "Highest frequency (Hz) to measure:", value=40, min=0, max=100, persist=true, style="slider")
#@ String (visibility=MESSAGE, value="It is recommended to down sample by a factor of 2 (half size) to remove image noise...", required=false) msg
#@ Integer 	Downsample (label = "Downsample? 1 = Full size, 2 = half size, 3 = 1/3 rd etc...", value=2, min=1, max=25, persist=true)
#@ String NormaliseIllumination	(choices={"Yes","No"}, label = "Normalise illumination between frames? Recommended unless using a computer generated test file", value="Yes", persist=false, style="radioButtonHorizontal")
#@ String (visibility=MESSAGE, value="Static background works best with the number of frames that capture a 1-2 oscillations...", required=false) msg2
#@ Integer 	BackgroundCreationFrames	(label = "How many frames to average in time to create static background?", value=20, min=0, max=100, persist=true, style="slider")
#@ String SaveStatic	(choices={"Yes","No"}, label = "Save static movie? Only needed for parameter optimisation", value="No", persist=false, style="radioButtonHorizontal")
#@ String (visibility=MESSAGE, value="Smooth noise between time frames? 0 = no smoothing, 1 or 2 = average of current frame with 1 or 2 frames before and after", required=false) msg3
#@ Integer 	SmoothingFrames	(label = "How many frames to average in time?", value=1, min=0, max=2, persist=true, style="slider")
#@ String (visibility=MESSAGE, value="What amplitude should the prominent frequency be vs background noise? Lower values will detect more background signals", required=false) msg4
#@ Float 	noise (label = "Primary Frequency prominence (2 to 100):", style="format:0.0", value=20, min=1.9, max=100, stepSize=0.1, persist=true)
//#@ Float 	outliers (label = "Proportion of low strength (background) frequencies used to set primary frequency cut off", style="slider,format:0.00", value=0.8, min=0.2, max=0.95, stepSize=0.05, persist=true)

outliers = 0.8;//v205 this takes the lowest 80% of fourier frequencies (by amplitude /RMS) to set the mean background frequency amplitude.

print("***Image Processing parameters***");
	MonthNames = newArray("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
     DayNames = newArray("Sun", "Mon","Tue","Wed","Thu","Fri","Sat");
     getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
     TimeString ="Date: "+DayNames[dayOfWeek]+" ";
     if (dayOfMonth<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+dayOfMonth+"-"+MonthNames[month]+"-"+year+"\nTime: ";
     if (hour<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+hour+":";
     if (minute<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+minute+":";
     if (second<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+second;
	 print(TimeString);
print("Capture speed of movies (FPS / Hz) = "+InputFPS);
print("Frequency cut off (Hz):\nLow: "+LowerHzCutOff+"\nHigh: "+UpperHzCutOff);
print("Image down-sampled in x & y by a factor of "+Downsample);
print("Static background created by averaging over "+BackgroundCreationFrames+" timepoints (frames)");
print("Image smoothed over time by smoothing current frame +/- "+SmoothingFrames+" frames");
print("Dominant frequency amplitude set at "+noise+" times the mean background frequency amplitude");
//print("Mean background frequency taken from the lowest "+outliers*100+"% of frequency signal amplitudes"); // v206 deleted

setBatchMode(true);

listdir = getFileList(input); 
for (i = 0; i < listdir.length; i++) { 
	//V201 ADD BACK THE CLEAR MEMORY
		run("Close All"); // v113 run this to succesfully clear memory. v114 moved up one line
		run("Collect Garbage"); // try and clear the memory
	
       path = input + File.separator + listdir[i];
      //if (File.isDirectory(path)){
       	run("Bio-Formats Importer", "open=[path] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
       	time1=getTime();
 //v003 move this to start  	
       	FF = newArray(); // reset arrays
		profile = newArray();
		maxima = newArray();
		RMS = newArray();
		
ImageName = File.nameWithoutExtension; //v206
//ImageName = substring(getTitle(), 0, lengthOf(getTitle())-4);

print("\n******\n"+ImageName+" analysis..."); //v113

// make the image 2^n frames. 256 512 1024 etc
getDimensions(width, height, channels, slices, frames);
//v207b squish 3 ch images back to one channel
if (channels >1) {run("Arrange Channels...", "new=2");
run("Grays");
}

if (slices > frames) { // if the movie progresses through slices rather than frames, take the number of frames as the timepoints.
	frames = slices;
}
if (frames < 128) {print("Image duration too short, must be at least 128 frames / timepoints. \n"+ImageName+" not analysed.");}

else{ // v113 only run the file if it is at least 128 frames
if (frames == 128) {
print("\nFile of "+frames+ " (2^6) frames is suitable for analysis.");	
}
else if (frames < 256) {run("Slice Keeper", "first=1 last=128 increment=1"); 
print("\nFile of "+frames+ " frames has been reduced to 128 frames for fourier analysis.");
} //v113 these lines inserted
else if (frames == 256) {
print("\n'File of "+frames+ " (2^7) frames is suitable for analysis.");	
}
else if (frames < 512) {run("Slice Keeper", "first=1 last=256 increment=1"); // set 1-256 (= 256 frames) v102
print("\nFile of "+frames+ " frames has been reduced to 256 frames for fourier analysis.");
}
else if (frames == 512) {
print("\nFile of "+frames+ " (2^8) frames is suitable for analysis.");	
}
else if (frames < 1024) {run("Slice Keeper", "first=1 last=512 increment=1");
print("\nFile of "+frames+ " frames has been reduced to 512 frames for fourier analysis.");
}
else if (frames == 1024) {
print("\nFile of "+frames+ " (2^9) frames is suitable for analysis.");	
}
else if (frames < 2048) {run("Slice Keeper", "first=1 last=1024 increment=1");
print("\nFile of "+frames+ " frames has been reduced to 1024 frames for fourier analysis.");
}
else if (frames == 2048) {
print("\nFile of "+frames+ " (2^10) frames is suitable for analysis.");	
}
else if (frames > 2048) {run("Slice Keeper", "first=1 last=2048 increment=1");
print("\nFile of "+frames+ " frames has been reduced to 2048 frames for fourier analysis.");
}

getDimensions(width, height, channels, slices, frames); //v002
duration = frames * slices / InputFPS;

print("This gives "+duration+" seconds of '"+ImageName+"' for analysis.");//v102, v113
MaxMeasurableFreq = minOf((frames*slices)/2, InputFPS/2); //v211
print("Maximum measurable frequency is "+MaxMeasurableFreq+" Hz."); //v211
HzResolution = 2 * MaxMeasurableFreq / (frames * slices); // v211
print("Frequency resolution is "+HzResolution+" Hz (we can measure beating in "+HzResolution+" Hz increments)."); //v211
rename("InputImage");

if (Downsample >= 2) {
	print("\nDownsampling by a factor of "+Downsample+" to speed up processing...");
	run("Scale...", "x="+1/Downsample+" y="+1/Downsample+" interpolation=Bilinear average process create"); // need to close orig image and rename the smaller one
	close("InputImage");
	selectWindow("InputImage-1");
	rename("InputImage");
	time2 = getTime();
	elapsed = (time2-time1)/1000;
	print("Completed in "+elapsed+" seconds.");
}

print("\nExtracting moving structures by removing a time smoothed ("+BackgroundCreationFrames+" frames) background \n& reslicing to speed up subsequent steps...");//v114


time2 = getTime();
run("Duplicate...", "title=meanZbackground duplicate");
run("Mean 3D...", "x=0 y=0 z="+BackgroundCreationFrames);
imageCalculator("Subtract create 32-bit stack", "InputImage","meanZbackground");

//v111 normalise the moving image
// v203 this will make every slice = zero if each image is a single pixel value i.e a test image with a perfect sigma wave
	if(NormaliseIllumination=="Yes"){		//v205
	for (j = 1; j <= nSlices; j++) {
	    setSlice(j);
	   getStatistics(area, mean, min, max, std, histogram);
// v203
		   if(mean>0) run("Subtract...", "value="+mean+" slice");
		   if(mean<0) run("Add...", "value="+Math.abs(mean)+" slice");
	}
	}
	
run("Reslice [/]...", "output=1.000 start=Top avoid"); 
time3 = getTime();
elapsed = (time3-time2)/1000;
print("Moving structures isolated in "+elapsed+" seconds.");

if (SmoothingFrames >= 1) {																			// ?? quicker without reslice first?
	print("\nSmoothing over time with "+SmoothingFrames+" frame(s) before and after each time point.");
	print("See progress of 3D filtering in ImageJ status bar...");
	run("Mean 3D...", "x=0 y="+SmoothingFrames+" z=0");
	time4 = getTime(); //v114
	elapsed = (time4-time3)/1000;
	print("Time smoothed in "+elapsed+" seconds.");//v114
	}
else time4 = getTime();

rename(ImageName+"-Processed&Resliced");
selectWindow("meanZbackground");
rename(ImageName+"-static");
selectWindow("Result of InputImage");
rename(ImageName+"-moving");
selectWindow("InputImage");
if (Downsample >= 2) rename(ImageName+"-downsampled");
else rename(ImageName);
      
print("\nPerforming Fourier analysis, hang tight your PC is working hard! See progress bar in Fiji window...");

selectWindow(ImageName+"-Processed&Resliced"); // this is resliced!!

getDimensions(width, height, channels, slices, frames);

Score = newArray(width*slices); // this array stores the RMS of the maximal frequency //v206 increased size to width*slices
BeatFreq = newArray(width*slices);
c = 0; // start a counter for number of pixels analysed

// v202 filtered power spectrum only
newImage("Power Spectrum filtered", "32-bit black", height/2, 1, 1); //v202

//time6 = getTime();

// Fourier analysis loops. Sequentially over every pixel...
	for (z = 0; z < slices; z++) {
	//v201 reselect the correct image!
	selectWindow(ImageName+"-Processed&Resliced"); // this is resliced!!	
	
	Stack.setSlice(z);
	showProgress(z, slices); // will show progress each time it does one z slice

		for (x = 0; x < width; x++) {
			//showProgress(x/width);	// not really working less than 30ms per cycle
			selectWindow(ImageName+"-Processed&Resliced");
			makeLine(x, 0, x, height-1, 1); // use the height of the image for the line to measure, MINUS 1 as the line cound (y) starts at ZERO not ONE v204
			profile = getProfile();
		
			FF = Array.fourier(profile, windowType); // v003 change to "Hann" windowing

Array.getStatistics(FF, min, max, mean, stdDev); // v004 get the SD
/*
 * Do a filter on the fourier data, ignore maxima above & below preset thresholds. Then take the best maxima that is remaining
 * Array.findMaxima Returns an array holding the peak positions (sorted with descending strength).
 */

	keep = (lengthOf(FF) * outliers); // proportion of Fourier array RMS values to keep i.e discard the top 10% and take the mean of the rest = background frequencies
	
	FFsorted = Array.copy(FF);// copy the fourier array
  	Array.sort(FFsorted);//sort the copy low to high
  	FFshorter = Array.trim(FFsorted, keep); // remove the top x % of values
  	Array.getStatistics(FFshorter, min, max, mean, stdDev); //get the mean
  	
					maxima = Array.findMaxima(FF, noise*mean, 1); //v101 use the noise * mean RMS of background as a cut off
					
					Array.getStatistics(maxima, min, max, mean, stdDev); // min = Infinity if array is empty.
//v202 add the filtered power spectrum in this if loop?					
					if(min != "Infinity"){
								PeakF = maxima[0]; // this is the postion of the maxima of the array FF
								BeatF = PeakF/duration;							//v103 divide by 2 when using 2^n frame input files // v204 no need to divide by 2 FF array now correct length
								Confidence = FF[PeakF]; // this is the value of the array at postion PeakF
								
								selectWindow("Power Spectrum filtered");//v202
									for (ps = 0; ps < height/2; ps++) {
									setPixel(ps, 0, FF[ps]);
   								}	
   								run("Add Slice"); // there will be one too many slices at the end.

							}
						else{ BeatF = 0; Confidence = 0;}
						
						Score[c] = Confidence;
						BeatFreq[c] = BeatF;	
						c=c+1; // add 1 to the counter
				
		FF = newArray(); // reset arrays added as v004a 
		profile = newArray();
		maxima = newArray();
		RMS = newArray();
		FFsorted =newArray();
		FFshorter = newArray();
		} //x loop
		
	} // z slices loop
		
	
time5 = getTime();

FFelapsed = (time5-time4)/1000;

print(c+" pixels with "+height+" timepoints analysed by Fourier transform in "+FFelapsed+" seconds."); //v115 height in resliced image is the number of time points analysed

//v202 save the power spectrum

selectWindow("Power Spectrum filtered");
if(nSlices >1) {run("Delete Slice"); // delete the last empty slice v202  v204 add the if command
run("Z Project...", "projection=[Average Intensity]");
}//v102


//v207 make a graph
run("Select All");
powerspectrum = getProfile();

axis = newArray(lengthOf(powerspectrum));
  for (p=0; p<lengthOf(powerspectrum); p++)
    axis[p] = p/duration; // divide by duration to get the value in Hz
  Plot.create("Power Spectrum: ", "frequency (Hz)", "amplitude (RMS)", axis, powerspectrum);
  Plot.show();
  saveAs("Jpeg",  output+File.separator+ImageName+"-PowerSpectrum");

selectWindow("AVG_Power Spectrum filtered"); //v211
Stack.setXUnit("Hz");//v211
run("Properties...", "channels=1 slices=1 frames=1 pixel_width="+HzResolution+" pixel_height=1 voxel_depth=1");
saveAs("Tiff",  output+File.separator+ImageName+"-Power_Spectrum"); //v211

// make images of the results

// BeatFreq
c=0; // reset the counter
newImage(ImageName+"FreqMap", "32-bit black", width, slices, 1); // slices (z) in resliced image is the original image height (y)

Array.getStatistics(BeatFreq, min, max, mean, stdDev); // v108 make sure there are some results greater than the lowerHz cut off
if (max > LowerHzCutOff){
	for (y = 0; y < slices; y++) {
		for (x = 0; x < width; x++) {
		setPixel(x, y, BeatFreq[c]);
		c = c+1;
		}
		//run("Enhance Contrast", "saturated=0.35");
		//run("Fire");
	}		
	run("Median...", "radius=0.5"); //v006
	setThreshold(0, LowerHzCutOff);//v112 select anything below lower Hz and delete
	run("Create Selection");
	resetThreshold;
	if (selectionType()>=0) run("Clear"); //v112 selection type -1 = no selection
	//run("Select None");
	setThreshold(UpperHzCutOff, 1000000); //v112 select anything above upper Hz and delete
	run("Create Selection");
	resetThreshold;
	if (selectionType()>=0) run("Clear"); //v112
}			// v108 skip these lines in the 'if' section if no results
	setThreshold(LowerHzCutOff, UpperHzCutOff);//v112 select the beating pixels within range to measure just these
	run("Create Selection");//v112
	setMinAndMax(0, UpperHzCutOff); // moved v102 & V208 make the lower display range value = zero
	run("Fire");
	//run("Histogram", "bins="+height/2+" x_min=0 x_max=&UpperHzCutOff y_max=Auto"); // v107 set bins to 400 to give better graphs //v207 set bins = height/2 = number of fourier measurements
	//v211
	run("Histogram", "bins="+Math.round(UpperHzCutOff/HzResolution)+" x_min=0 x_max=&UpperHzCutOff y_max=Auto");
	rename("Histogram");
	
	selectWindow(ImageName+"FreqMap");//v006
	rename(ImageName); //v210 to just have the name in the results table
	run("Select None"); //v209
	setThreshold(LowerHzCutOff, UpperHzCutOff); //v209
	run("Measure"); // measure the frequencies and area in one step v209 by setting a threshold and using 'Limit to threshold' in measurements
	rename(ImageName+"FreqMap");
	
//v209
	if (dayOfMonth<10) {TimeString2 ="0";
     TimeString2 = TimeString2+dayOfMonth+"-"+MonthNames[month]+"-"+year+"-";
     }
     else{TimeString2 = ""+dayOfMonth+"-"+MonthNames[month]+"-"+year+"@";}
     if (hour<10) {TimeString2 = TimeString2+"0";} //v016 second Timestring changed to Timestring2
     TimeString2 = TimeString2+hour;
     if (minute<10) {TimeString2 = TimeString2+"0";}
     TimeString2 = TimeString2+minute;
		
saveAs("results",  output+File.separator+"Results_Table-"+TimeString2+".csv"); 
// v209
	
// Confidence (Fourier RMS)
c=0; // reset the counter
newImage("Confidence", "32-bit black", width, slices, 1); // slices in resliced is the original image height (y)

	for (y = 0; y < slices; y++) {
		for (x = 0; x < width; x++) {
		setPixel(x, y, Score[c]);
		c = c+1;
		}
		run("Enhance Contrast", "saturated=0.35");
		run("Fire");
	}
	run("Median...", "radius=0.5"); //v006
	
// save the images. Static movie. Moving. Beat Freq. Confidence. Processed & Resliced. Histogram. 
if(SaveStatic=="Yes"){ // v211 option to save
	selectWindow(ImageName+"-static");
	saveAs("Tiff",  output+File.separator+ImageName+"-Static");
}
selectWindow(ImageName+"FreqMap");								//v006
saveAs("Tiff",  output+File.separator+ImageName+"-FreqMap");
selectWindow("Confidence");
saveAs("Tiff",  output+File.separator+ImageName+"-Confidence");
selectWindow(ImageName+"-Processed&Resliced");	//v003
run("Select All");
run("Reslice [/]...", "start=Top avoid"); 							//v110
saveAs("Tiff",  output+File.separator+ImageName+"-moving&processed");
selectWindow("Histogram");//v006
setBatchMode("show");
	Table.showHistogramTable;	//v103 save histogram list CRASHING HERE? //v103 have to jump out of batch mode
	Table.rename("Histogram", "HistogramPlot"); 
	Table.deleteColumn("index");
	Table.renameColumn("bin start", "Frequency");
	saveAs("results",  output+File.separator+ImageName+"FreqPlot.csv");
	close(ImageName+"FreqPlot.csv");
selectWindow("Histogram");//v006
saveAs("Jpeg",  output+File.separator+ImageName+"-FrequencyPlot");
close("Histo*"); // close the histogram window
setBatchMode(true);


time6 = getTime();
elapsed = (time6-time1)/1000;
print("\nTotal processing time for '"+ImageName+"' was "+elapsed+" seconds.\nMemory Status:");

run("Collect Garbage"); // try and clear the memory
IJ.freeMemory(); // report the memory usage

} //v113 internal loop only if file >=128 frames
}

selectWindow("Results");  
Table.renameColumn("Mean", "Mean Beating Frequency (Hz)");
Table.renameColumn("%Area", "Beating cilia area (%)"); //v211
Table.renameColumn("MinThr", "Lower cut off (Hz)");
Table.renameColumn("MaxThr", "Upper cut off (Hz)");
updateResults();		
saveAs("results",  output+File.separator+"Results_Table-"+TimeString2+".csv");

// add a time stamp v105 correct at end v209	
 getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec); //v209 correct time at end
if (dayOfMonth<10) {TimeString2 ="0";
     TimeString2 = TimeString2+dayOfMonth+"-"+MonthNames[month]+"-"+year+"@";
     }
     else{TimeString2 = ""+dayOfMonth+"-"+MonthNames[month]+"-"+year+"@";}
     if (hour<10) {TimeString2 = TimeString2+"0";} //v016 second Timestring changed to Timestring2
     TimeString2 = TimeString2+hour;
     if (minute<10) {TimeString2 = TimeString2+"0";}
     TimeString2 = TimeString2+minute;
 


print("\nAnalysis finished "+TimeString2+".");
selectWindow("Log");
saveAs(output+File.separator+"Log-"+TimeString2+".txt");

exit("All done!\n"+i+" images processed");

}

