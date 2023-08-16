//Lab 1: Intro to Remote Sensing and EE


// 1 Searching (and finding) Landsat imagery


// Note that we need to cast the result of first() to Image.
var image = ee.Image(landsat

    // Filter to get only images in the specified range.
    .filterDate('2014-01-01', '2014-12-31')

    // Filter to get only images at the location of the point.
    .filterBounds(point)

    // Sort the collection by a metadata property.
    .sort('CLOUD_COVER')

    // Get the first image out of this collection.
    .first());

// Print the image to the console.
print('A Landsat scene:', image);

//2 Visualizing Landsat imagery

// Define visualization parameters in a JavaScript dictionary.
var trueColor = {
  bands: ['B4', 'B3', 'B2'],
  min: 4000,
  max: 12000
};

// Add the image to the map, using the visualization parameters.
Map.addLayer(image, trueColor, 'true-color image');

// Define false-color visualization parameters.
var falseColor = {
  bands: ['B5', 'B4', 'B3'],
  min: 4000,
  max: 13000
};

// Add the image to the map, using the visualization parameters.
Map.addLayer(image, falseColor, 'false-color composite');

//3 Plot at-sensor radiance at several locations

// Use these bands. 
var bands = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10', 'B11'];

// Get an image that contains only the bands of interest.
var dnImage = image.select(bands);

// Apply the transformation.
var radiance = ee.Algorithms.Landsat.calibratedRadiance(dnImage);

// Display the result.
var radParams = {bands: ['B4', 'B3', 'B2'], min: 0, max: 100};
Map.addLayer(radiance, radParams, 'radiance');

//4 Plot top-of-atmosphere (TOA) reflectance at several locations

var toaImage = ee.Image('LANDSAT/LC8_L1T_TOA/LC80440342014077LGN00');

Map.addLayer(toaImage, {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.3}, 'toa');

// Hardcode a point in Golden Gate Park.
var ggPark = ee.Geometry.Point([-122.4860, 37.7692]);

// Define reflective bands as bands B1-B7.  See the docs for slice().
var reflectiveBands = bands.slice(0, 7);

// See http://landsat.usgs.gov/band_designations_landsat_satellites.php
var wavelengths = [0.44, 0.48, 0.56, 0.65, 0.86, 1.61, 2.2];

// Select only the reflectance bands of interest.
var reflectanceImage = toaImage.select(reflectiveBands);

// Define an object of customization parameters for the chart.
var options = {
  title: 'Landsat 8 TOA spectrum in Golden Gate Park',
  hAxis: {title: 'Wavelength (micrometers)'},
  vAxis: {title: 'Reflectance'},
  lineWidth: 1,
  pointSize: 4
};

// Make the chart, using a 30 meter pixel.
var chart = Chart.image.regions(
  reflectanceImage, ggPark, null, 30, null, wavelengths)
    .setOptions(options);

// Display the chart.
print(chart);

//End






//Lab 2: Characteristics of remotely sensed data


/*The following code summarizes the Lab 2: Characteristics of remotely sensed data of the Google Earth Engine
Educational Resources. Detailed description of the Lab can be found at https://goo.gl/hAiyeJ. The purpose of this lab 
is to demonstrate concepts of spatial, spectral, temporal and radiometric resolution, introducing image data from 
several sensors aboard various platforms.

Datasets
var myd09 = ee.ImageCollection("MODIS/MYD09GA"),
    mss = ee.ImageCollection("LANDSAT/LM5_L1T"),
    tm = ee.ImageCollection("LANDSAT/LT5_L1T_TOA"),
    naip = ee.ImageCollection("USDA/NAIP/DOQQ");
 
 1 Spatial resolution */
 
 //A.MODIS
 
// Define a region of interest as a point at SFO airport.
var sfoPoint = ee.Geometry.Point(-122.3774, 37.6194);

// Center the map at that point.
Map.centerObject(sfoPoint, 16);

// Get a surface reflectance image from the MODIS MYD09GA collection.
var modisImage = ee.Image(myd09.filterDate('2011-08-28').first());

// Use these MODIS bands for red, green, blue, respectively.
var modisBands = ['sur_refl_b01', 'sur_refl_b04', 'sur_refl_b03'];

// Define visualization parameters for MODIS.
var modisVis = {bands: modisBands, min: 0, max: 3000};

// Add the MODIS image to the map.
Map.addLayer(modisImage, modisVis, 'MODIS');

// Get the scale of the data from the first band's projection:
var modisScale = modisImage.select('sur_refl_b01')
    .projection().nominalScale();

print('MODIS scale:', modisScale);

 //B.MSS
 
// Filter MSS imagery by location, date and cloudiness.
var mssImage = ee.Image(mss
    .filterBounds(Map.getCenter())
    .filterDate('2012-05-01', '2012-10-01')
    .sort('CLOUD_COVER')
    // Get the least cloudy image.
    .first());

// Display the MSS image as a color-IR composite.
Map.addLayer(mssImage, {bands: ['B3', 'B2', 'B1'], min: 0, max: 200}, 'MSS')

// Get the scale of the MSS data from its projection:
var mssScale = mssImage.select('B1')
    .projection().nominalScale();

print('MSS scale:', mssScale);

 //C.TM
 
// Filter TM imagery by location, date and cloudiness.
var tmImage = ee.Image(tm
    .filterBounds(Map.getCenter())
    .filterDate('2011-05-01', '2011-10-01')
    .sort('CLOUD_COVER')
    .first());
    
// Display the TM image as a color-IR composite.
Map.addLayer(tmImage, {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.4}, 'TM');

// Get the scale of the TM data from its projection:
var tmScale = tmImage.select('B1')
    .projection().nominalScale();

print('TM scale:', tmScale);

 //D.NAIP
 
// Get NAIP images for the study period and region of interest.
var naipImages = naip.filterDate('2012-01-01', '2012-12-31')
    .filterBounds(Map.getCenter());

// Mosaic adjacent images into a single image.
var naipImage = naipImages.mosaic();

// Display the NAIP mosaic as a color-IR composite.
Map.addLayer(naipImage, {bands: ['N', 'R', 'G']}, 'NAIP');

// Get the NAIP resolution from the first image in the mosaic.
var naipScale = ee.Image(naipImages.first()).projection().nominalScale();

print('NAIP scale:', naipScale);

// 2 Spatial resolution 

// Get the MODIS band names as a List
var modisBands = modisImage.bandNames();

// Print the list.
print('MODIS bands:', modisBands);

// Print the length of the list.
print('Length of the bands list:', modisBands.length());

// 3 Temporal resolution

 //A.MODIS
 
// Filter the MODIS mosaics to one year.
var modisSeries = myd09.filterDate('2011-01-01', '2011-12-31');

// Print the filtered MODIS ImageCollection.
print('MODIS series:', modisSeries);

 //B.LANDSAT

// Filter to get a year's worth of TM scenes.
var tmSeries = tm
    .filterBounds(Map.getCenter())
    .filterDate('2011-01-01', '2011-12-31');

// Print the filtered TM ImageCollection. 
print('TM series:', tmSeries);

var getDate = function(image) {
  // Note that you need to cast the argument
  var time = ee.Image(image).get('system:time_start');
  // Return the time (in milliseconds since Jan 1, 1970) as a Date
  return ee.Date(time);
};

var dates = tmSeries.toList(100).map(getDate);

print(dates);

// 4 Radiometric resolution

// 5 Orbits and sensor motion

var aquaImage = ee.Image(myd09.filterDate('2013-09-01').first());

// Zoom to global level.
Map.setCenter(81.04, 0, 3);

// Display the sensor-zenith angle of the Aqua imagery.
var szParams = {bands: 'SensorZenith', min: 0, max: 70*100};
Map.addLayer(aquaImage, szParams, 'Aqua sensor-zenith angle')

var aquaOrbit = ee.FeatureCollection('ft:1ESvPygQ76WvVflKMN2nc14sS2wwtzqv3j2ueTqg');

Map.addLayer(aquaOrbit, {color: 'FF0000'}, 'Aqua Positions');

// Load Landsat ETM+ data directly, filter to one day.
var landsat7 = ee.ImageCollection('LANDSAT/LE7')
.filterDate('2013-09-01', '2013-09-02');

// Display the images by specifying one band and a single color palette.
Map.addLayer(landsat7, {bands: 'B1', palette: 'blue'}, 'Landsat 7 scenes');

//End






//Lab 3: Digital imagery and image processing

//1 What is a digital image?

// Get a single NAIP image over the area of interest.
var image = ee.Image(naip
    .filterBounds(point)
    .sort('system:time_start', false)
    .first());
    
// Print the image to the console.
print('Inspect the image object:', image);

// Display the image with the default visualization.
Map.centerObject(point, 18);
Map.addLayer(image, {}, 'Original image');

// Display the projection of band 0.
print('Inspect the projection of band 0:', image.select(0).projection());

//2 Digital image visualization

// Display gamma stretches of the input image.
Map.addLayer(image.visualize({gamma: 0.5}), {}, 'gamma = 0.5');
Map.addLayer(image.visualize({gamma: 1.5}), {}, 'gamma = 1.5');


// Define a RasterSymbolizer element with '_enhance_' for a placeholder.
var histogram_sld =
  '<RasterSymbolizer>' +
    '<ContrastEnhancement><Histogram/></ContrastEnhancement>' +
    '<ChannelSelection>' +
      '<RedChannel>' +
        '<SourceChannelName>R</SourceChannelName>' +
      '</RedChannel>' +
      '<GreenChannel>' +
        '<SourceChannelName>G</SourceChannelName>' +
      '</GreenChannel>' +
      '<BlueChannel>' +
        '<SourceChannelName>B</SourceChannelName>' +
      '</BlueChannel>' +
    '</ChannelSelection>' +
  '</RasterSymbolizer>';

// Display the image with a histogram equalization stretch.
Map.addLayer(image.sldStyle(histogram_sld), {}, 'Equalized');



//3 Linear filtering

// Print a uniform kernel to see its weights.
print('A uniform kernel:', ee.Kernel.square(2));

// Define a square, uniform kernel.
var uniformKernel = ee.Kernel.square({
  radius: 2,
  units: 'meters',
});

// Filter the image by convolving with the smoothing filter.
var smoothed = image.convolve(uniformKernel);
Map.addLayer(smoothed, {min: 0, max: 255}, 'smoothed image');

// Print a Gaussian kernel to see its weights.
print('A Gaussian kernel:', ee.Kernel.gaussian(2));

// Define a square Gaussian kernel:
var gaussianKernel = ee.Kernel.gaussian({
  radius: 2,
  units: 'meters',
});

// Filter the image by convolving with the Gaussian filter.
var gaussian = image.convolve(gaussianKernel);
Map.addLayer(gaussian, {min: 0, max: 255}, 'Gaussian smoothed image');

// Define a Laplacian filter.
var laplacianKernel = ee.Kernel.laplacian8();

// Print the kernel to see its weights.
print(laplacianKernel);

// Filter the image by convolving with the Laplacian filter.
var edges = image.convolve(laplacianKernel)
		     .reproject('EPSG:26910', null, 1);
Map.addLayer(edges, {min: 0, max: 255}, 'Laplacian filtered image');

// Compute the image gradient in the X and Y directions.
var xyGrad = image.select('N').gradient();

// Compute the magnitude of the gradient.
var gradient = xyGrad.select('x').pow(2)
          .add(xyGrad.select('y').pow(2)).sqrt()
    .reproject('EPSG:26910', null, 1);

// Compute the direction of the gradient.
var direction = xyGrad.select('y').atan2(xyGrad.select('x'))
.reproject('EPSG:26910', null, 1);

// Display the results.
Map.setCenter(-122.054, 37.7295, 10);
Map.addLayer(direction, {min: -3, max: 3, format: 'png'}, 'direction');
Map.addLayer(gradient, {min: -10, max: 50, format: 'png'}, 'gradient');

// Define a "fat" Gaussian kernel.
var fat = ee.Kernel.gaussian({
  radius: 3,
  sigma: 3,
  magnitude: -1,
  units: 'meters'
});

// Define a "skinny" Gaussian kernel.
var skinny = ee.Kernel.gaussian({
  radius: 3,
  sigma: 0.5,
  units: 'meters'
});

// Compute a difference-of-Gaussians (DOG) kernel.
var dog = fat.add(skinny);

// Add the DoG filtered image to the original image.
var sharpened = image.add(image.convolve(dog));
Map.addLayer(sharpened, {min: 0, max: 255}, 'Edges enhanced');


//4 Non-linear filtering

var median = image.reduceNeighborhood({
  reducer: ee.Reducer.median(), 
  kernel: uniformKernel
});

Map.addLayer(median, {min: 0, max: 255}, 'Median');

// Create and display a simple two-class image.
var veg = image.select('N').gt(200);

// Display the two-class (binary) result.
var binaryVis = {min: 0, max: 1, palette: ['black', 'green']};
Map.addLayer(veg, binaryVis, 'veg');


// Compute the mode in each 5x5 neighborhood and display the result.
var mode = veg.reduceNeighborhood({
  reducer: ee.Reducer.mode(), 
  kernel: uniformKernel
});

Map.addLayer(mode, binaryVis, 'mode');


// Dilate by takaing the max in each 5x5 neighborhood.
var max = veg.reduceNeighborhood({
  reducer: ee.Reducer.max(), 
  kernel: uniformKernel
});

Map.addLayer(max, binaryVis, 'max');

// Erode by takaing the min in each 5x5 neighborhood.
var min = veg.reduceNeighborhood({
  reducer: ee.Reducer.min(), 
  kernel: uniformKernel
});

Map.addLayer(min, binaryVis, 'min');

// Perform an opening by dilating the eroded image.
var opened = min.reduceNeighborhood({
  reducer: ee.Reducer.max(), 
  kernel: uniformKernel
});

Map.addLayer(opened, binaryVis, 'opened');

// Perform a closing by eroding the dilated image.
var closed = max.reduceNeighborhood({
  reducer: ee.Reducer.min(), 
  kernel: uniformKernel
});

Map.addLayer(closed, binaryVis, 'closed');

//5 Texture

// Define a big neighborhood with a 7-meter radius kernel.
var bigKernel = ee.Kernel.square({
  radius: 7,
  units: 'meters'
});

// Compute SD in a neighborhood.
var sd = image.reduceNeighborhood({
  reducer: ee.Reducer.stdDev(), 
  kernel: bigKernel
});

Map.addLayer(sd, {min: 0, max: 70}, 'SD');


// Compute entropy in a neighborhood.
var entropy = image.entropy(bigKernel);

Map.addLayer(entropy, {min: 1, max: 5}, 'entropy');

// Use the GLCM to compute a large number of texture measures.
var glcmTexture = image.glcmTexture(7);

// Display the 'contrast' results for the red, green and blue bands.
var contrastVis = {
  bands: ['R_contrast', 'G_contrast', 'B_contrast'], 
  min: 40,
  max: 2000
};

Map.addLayer(glcmTexture, contrastVis, 'contrast');

// Create a list of weights for a 9x9 kernel.
var list = [1, 1, 1, 1, 1, 1, 1, 1, 1];
// The center of the kernel is zero.
var centerList = [1, 1, 1, 1, 0, 1, 1, 1, 1];
// Assemble a list of lists: the 9x9 kernel weights as a 2-D matrix.
var lists = [list, list, list, list, centerList, list, list, list, list];
// Create the kernel from the weights.
// Non-zero weights represent the spatial neighborhood.
var kernel = ee.Kernel.fixed(9, 9, lists, -4, -4, false);

// Use the max among bands as the input.
var maxBands = image.reduce(ee.Reducer.max());

// Convert the neighborhood into multiple bands.
var neighs = maxBands.neighborhoodToBands(kernel);

// Compute local Geary's C, a measure of spatial association.
var gearys = maxBands.subtract(neighs).pow(2).reduce(ee.Reducer.sum())
             .divide(Math.pow(9, 2));

Map.addLayer(gearys, {min: 20, max: 2500}, "Geary's C");


//6 Resampling and Reprojection

// Zoom all the way in.
Map.centerObject(point, 21);

// Display edges computed on a reprojected image.
Map.addLayer(image.convolve(laplacianKernel), {min: 0, max: 255}, 
    'Edges with little screen pixels');

// Display edges computed on the image at native resolution.
Map.addLayer(edges, {min: 0, max: 255}, 
    'Edges with 1 meter pixels');

// Resample the image with bilinear instead of nearest neighbor.
var bilinearResampled = image.resample('bilinear');
Map.addLayer(bilinearResampled, {}, 'input image, bilinear resampling');

// Resample the image with bicubic instead of nearest neighbor.
var bicubicResampled = image.resample('bicubic');
Map.addLayer(bicubicResampled, {}, 'input image, bicubic resampling');

//End

// Lab 4: Spectral indices and transformations


/*The following code summarizes the Lab 4: Spectral indices and
transformations. Detailed description of the Lab can be found at https://goo.gl/hAiyeJ. 
The purpose of this lab is introducing to methods for creating vegetation, water
snowm bare soil and burned area indices.*/

// 1 Spectral indices

// A.NDVI


var image = ee.Image(landsat8
    .filterBounds(point)
    .filterDate('2015-06-01', '2015-09-01')
    .sort('CLOUD_COVER')
    .first())

var trueColor = {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.3};
Map.addLayer(image, trueColor, 'image');

var ndvi = image.normalizedDifference(['B5', 'B4']);

var vegPalette = ['red', 'blue', 'yellow', 'green'];
Map.addLayer(ndvi, {min: -1, max: 1, palette: vegPalette}, 'NDVI');

// B.EVI

var evi = image.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
      'NIR': image.select('B5'),
      'RED': image.select('B4'),
      'BLUE': image.select('B2')
});

Map.addLayer(evi, {min: -1, max: 1, palette: vegPalette}, 'EVI');

// C.NDWI

var ndwi = image.normalizedDifference(['B5', 'B6']);

var waterPalette = ['red', 'yellow', 'green', 'blue'];
Map.addLayer(ndwi, {min: -0.5, max: 1, palette: waterPalette}, 'NDWI');

// D.NDWBI

var ndwi = image.normalizedDifference(['B3', 'B5']);
Map.addLayer(ndwi, {min: -1, max: 0.5, palette: waterPalette}, 'NDWBI');

// E.NDBI

var ndbi = image.normalizedDifference(['B6', 'B5']);
var barePalette = waterPalette.slice().reverse();
Map.addLayer(ndbi, {min: -1, max: 0.5, palette: barePalette}, 'NDBI');

// F.BAI

var burnImage = ee.Image(landsat8
    .filterBounds(ee.Geometry.Point(-120.083, 37.850))
    .filterDate('2013-08-17', '2013-09-27')
    .sort('CLOUD_COVER')
    .first());
Map.addLayer(burnImage, trueColor, 'burn image');

var bai = burnImage.expression(
    '1.0 / ((0.1 - RED)**2 + (0.06 - NIR)**2)', {
      'NIR': burnImage.select('B5'),
      'RED': burnImage.select('B4'),
});

var burnPalette = ['green', 'blue', 'yellow', 'red'];
Map.addLayer(bai, {min: 0, max: 400, palette: burnPalette}, 'BAI');

// G.NBRT

var nbrt = burnImage.expression(
  '(NIR - 0.0001 * SWIR * Temp) / (NIR + 0.0001 * SWIR * Temp)', {
    'NIR': burnImage.select('B5'),
    'SWIR': burnImage.select('B7'),
    'Temp': burnImage.select('B11')
});

Map.addLayer(nbrt, {min: 1, max: 0.9, palette: burnPalette}, 'NBRT');

// H.NDSI

var snowImage = ee.Image(landsat8
    .filterBounds(ee.Geometry.Point(-120.0421, 39.1002))
    .filterDate('2013-11-01', '2014-05-01')
    .sort('CLOUD_COVER')
    .first());

Map.addLayer(snowImage, trueColor, 'snow image');

var ndsi = snowImage.normalizedDifference(['B3', 'B6']);
var snowPalette = ['red', 'green', 'blue', 'white'];
Map.addLayer(ndsi, {min: -0.5, max: 0.5, palette: snowPalette}, 'NDSI');

// 2 Linar transformations

// A.Tasseled cap (TC)

var coefficients = ee.Array([
  [0.3037, 0.2793, 0.4743, 0.5585, 0.5082, 0.1863],
  [-0.2848, -0.2435, -0.5436, 0.7243, 0.0840, -0.1800],
  [0.1509, 0.1973, 0.3279, 0.3406, -0.7112, -0.4572],
  [-0.8242, 0.0849, 0.4392, -0.0580, 0.2012, -0.2768],
  [-0.3280, 0.0549, 0.1075, 0.1855, -0.4357, 0.8085],
  [0.1084, -0.9022, 0.4120, 0.0573, -0.0251, 0.0238]
]);

var tcImage = ee.Image(landsat5
    .filterBounds(point)
    .filterDate('2008-06-01', '2008-09-01')
    .sort('CLOUD_COVER')
    .first());

var bands = ['B1', 'B2', 'B3', 'B4', 'B5', 'B7'];
// Make an Array Image, with a 1-D Array per pixel.
var arrayImage1D = tcImage.select(bands).toArray();
// Make an Array Image with a 2-D Array per pixel, 6x1.
var arrayImage2D = arrayImage1D.toArray(1);

var componentsImage = ee.Image(coefficients)
  .matrixMultiply(arrayImage2D)
  // Get rid of the extra dimensions.
  .arrayProject([0])
  // Get a multi-band image with TC-named bands.
  .arrayFlatten(
    [['brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth']]);

var vizParams = {
  bands: ['brightness', 'greenness', 'wetness'],
  min: -0.1, max: [0.5, 0.1, 0.1]
};
Map.addLayer(componentsImage, vizParams, 'TC components');

// B.Principal Component Analysis (PCA)

var bands = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10', 'B11'];
var arrayImage = image.select(bands).toArray();

var covar = arrayImage.reduceRegion({
  reducer: ee.Reducer.covariance(),
  maxPixels: 1e9
});
var covarArray = ee.Array(covar.get('array'));

var eigens = covarArray.eigen();

var eigenVectors = eigens.slice(1, 1);

var principalComponents = ee.Image(eigenVectors)
.matrixMultiply(arrayImage.toArray(1));

var pcImage = principalComponents
  // Throw out an an unneeded dimension, [[]] -> [].
  .arrayProject([0])
  // Make the one band array image a multi-band image, [] -> image.
  .arrayFlatten([['pc1', 'pc2', 'pc3', 'pc4', 'pc5', 'pc6', 'pc7', 'pc8']]);
Map.addLayer(pcImage.select('pc3'), {}, 'PC');

// C.Spectral unmixing

var unmixImage = image.select(['B2', 'B3', 'B4', 'B5', 'B6', 'B7']);

Map.addLayer(image, {bands: ['B5', 'B4', 'B3'], max: 0.4}, 'false color');

print(Chart.image.regions(unmixImage, ee.FeatureCollection([
    ee.Feature(bare, {label: 'bare'}), 
    ee.Feature(water, {label: 'water'}),
    ee.Feature(veg, {label: 'vegetation'})]), 
  ee.Reducer.mean(), 30, 'label', [0.48, 0.56, 0.65, 0.86, 1.61, 2.2]));

var bareMean = unmixImage
.reduceRegion(ee.Reducer.mean(), bare, 30).values();
var waterMean = unmixImage
.reduceRegion(ee.Reducer.mean(), water, 30).values();
var vegMean = unmixImage
.reduceRegion(ee.Reducer.mean(), veg, 30).values();

var endmembers = ee.Array.cat([bareMean, vegMean, waterMean], 1);

var arrayImage = unmixImage.toArray().toArray(1);

var unmixed = ee.Image(endmembers).matrixSolve(arrayImage);

var unmixedImage = unmixed
.arrayProject([0])
.arrayFlatten([['bare', 'veg', 'water']]);

Map.addLayer(unmixedImage, {}, 'Unmixed');

// 3 The HSV transform

// Convert Landsat RGB bands to HSV
var hsv = image.select(['B4', 'B3', 'B2']).rgbToHsv();
// Convert back to RGB, swapping the image panchromatic band for the value.
var rgb = ee.Image.cat([
  hsv.select('hue'), 
  hsv.select('saturation'),
  image.select(['B8'])
]).hsvToRgb();
Map.addLayer(rgb, {max: 0.4}, 'Pan-sharpened');

// Lab 5: Supervised classification and regression


/*The following code summarizes the Lab 5: Supervised classification and regression. 
Detailed description of the Lab can be found at https://goo.gl/hAiyeJ. 
supervised classification and regression: prediction of nominal or numeric values of a 
geographic variable from other geographic variables.*/

//Datasets 

var mod44b = ee.ImageCollection("MODIS/051/MOD44B"),
    l5raw = ee.ImageCollection("LANDSAT/LT5_L1T"),
    bare = /* color: #d63000 */ee.FeatureCollection(
        [ee.Feature(
            ee.Geometry.Polygon(
                [[[-121.453857421875, 37.45959832290546],
                  [-121.4208984375, 37.46177847961746],
                  [-121.4263916015625, 37.49011473195046]]]),
            {
              "class": 0,
              "system:index": "0"
            })]),
    vegetation = /* color: #98ff00 */ee.FeatureCollection(
        [ee.Feature(
            ee.Geometry.Polygon(
                [[[-121.8878173828125, 37.34832607355296],
                  [-121.85211181640625, 37.34832607355296],
                  [-121.87957763671875, 37.37233994582318]]]),
            {
              "class": 1,
              "system:index": "0"
            })]),
    water = /* color: #0b4a8b */ee.FeatureCollection(
        [ee.Feature(
            ee.Geometry.Polygon(
                [[[-122.56072998046875, 37.39416407012379],
                  [-122.486572265625, 37.374522644077246],
                  [-122.51678466796875, 37.42252593456306]]]),
            {
              "class": 2,
              "system:index": "0"
            })]);

// 1. Introduction to classification and regression

// 2. Regression

// A.Ordinary Least Squares (OLS)

var tree = ee.Image(mod44b.sort('system:time_start', false).first());

var percentTree = tree.select('Percent_Tree_Cover')
    .where(tree.select('Percent_Tree_Cover').eq(200), 0);
Map.addLayer(percentTree, {max: 100}, 'percent tree cover');

var l5filtered = l5raw.filterDate('2010-01-01', '2010-12-31')
                      .filterMetadata('WRS_PATH', 'equals', 44)
                      .filterMetadata('WRS_ROW', 'equals', 34);

var landsat = ee.Algorithms.Landsat.simpleComposite({
  collection: l5filtered,
  asFloat: true
});

Map.addLayer(landsat, {bands: ['B4', 'B3', 'B2'], max: 0.3}, 'composite');

var predictionBands = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7'];

var trainingImage = ee.Image(1)
    .addBands(landsat.select(predictionBands))
    .addBands(percentTree);

var training = trainingImage.sample({
  region: l5filtered.first().geometry(), 
  scale: 30, 
  numPixels: 1000
});

var trainingList = ee.List(predictionBands)
    .insert(0, 'constant')
    .add('Percent_Tree_Cover');

var regression = training.reduceColumns({
  reducer: ee.Reducer.linearRegression(8), 
  selectors: trainingList
});

var coefficients = ee.Array(regression.get('coefficients'))
    .project([0])
    .toList();
var predictedTreeCover = 
ee.Image(1).addBands(landsat.select(predictionBands))
            .multiply(ee.Image.constant(coefficients))
            .reduce(ee.Reducer.sum())
            .rename('predictedTreeCover');
Map.addLayer(predictedTreeCover, {min: 0, max: 100}, 'prediction');

// B.Non-linear regression functions

var cartRegression = ee.Classifier.cart()
    .setOutputMode('REGRESSION')
    .train({
      features: training, 
      classProperty: 'Percent_Tree_Cover', 
      inputProperties: predictionBands
    });

var cartRegressionImage = landsat.select(predictionBands)
    .classify(cartRegression, 'cartRegression');

Map.addLayer(cartRegressionImage, {min: 0, max: 100}, 'CART regression');

// 3. Classification

var trainingFeatures = bare.merge(vegetation).merge(water);

var classifierTraining = landsat.select(predictionBands)
    .sampleRegions({
      collection: trainingFeatures, 
      properties: ['class'], 
      scale: 30
    });

var classifier = ee.Classifier.cart().train({
  features: classifierTraining, 
  classProperty: 'class', 
  inputProperties: predictionBands
});

var classified = landsat.select(predictionBands).classify(classifier);
Map.addLayer(classified, {min: 0, max: 2, palette: ['red', 'green', 'blue']}, 'classified');

// 4. Accuracy Assessment

var trainingTesting = classifierTraining.randomColumn();
var trainingSet = trainingTesting
.filter(ee.Filter.lessThan('random', 0.6));
var testingSet = trainingTesting
.filter(ee.Filter.greaterThanOrEquals('random', 0.6));

var trained = ee.Classifier.cart().train({
  features: trainingSet, 
  classProperty: 'class', 
  inputProperties: predictionBands
});

var confusionMatrix = ee.ConfusionMatrix(testingSet.classify(trained)
    .errorMatrix({
      actual: 'class', 
      predicted: 'classification'
    }));

print('Confusion matrix:', confusionMatrix);
print('Overall Accuracy:', confusionMatrix.accuracy());
print('Producers Accuracy:', confusionMatrix.producersAccuracy());
print('Consumers Accuracy:', confusionMatrix.consumersAccuracy());

//END

// Lab 6: Time Series Analysis

/*The following code summarizes the Lab 5: Time Series Analysis. 
Detailed description of the Lab can be found at https://goo.gl/hAiyeJ. 
The purpose of this lab is to establish a foundation for time series analysis 
of remotely sensed data, usually in the form of a temporally ordered stack of images.  
You will be introduced to concepts of smoothing, interpolation, linear modeling and phenology.  
At the completion of the lab, you will be able to perform analysis of multi-temporal data for 
determining trend and seasonality on a per-pixel basis.*/

// 1. Multi-temporal data in Earth Engine

// 2. Data preparation and preprocessing

// A.Load a time series of Landsat data

// B. Filtering, masking and preparing bands of interest

// This field contains UNIX time in milliseconds.
var timeField = 'system:time_start';
// Use this function to mask clouds in Landsat 8 imagery.
var maskClouds = function(image) {
  var quality = image.select('BQA');
  var cloud01 = quality.eq(61440);
  var cloud02 = quality.eq(53248);
  var cloud03 = quality.eq(28672);
  var mask = cloud01.or(cloud02).or(cloud03).not();
  return image.updateMask(mask);
};
// Use this function to add variables for NDVI, time and a constant
// to Landsat 8 imagery.
var addVariables = function(image) {
  // Compute time in fractional years since the epoch.
  var date = ee.Date(image.get(timeField));
  var years = date.difference(ee.Date('1970-01-01'), 'year');
  // Return the image with the added bands.
  return image
    // Add an NDVI band.
    .addBands(image.normalizedDifference(['B5', 'B4']).rename('NDVI'))
.float()
    // Add a time band.
    .addBands(ee.Image(years).rename('t').float())
    // Add a constant band.
    .addBands(ee.Image.constant(1));
};
// Remove clouds, add variables and filter to the area of interest.
var filteredLandsat = l8toa
  .filterBounds(roi)
  .map(maskClouds)
  .map(addVariables);

// C. Plot the time series at the location of interest

// Plot a time series of NDVI at a single location.
var l8Chart = ui.Chart.image.series(filteredLandsat.select('NDVI'), roi)
    .setChartType('ScatterChart')
    .setOptions({
      title: 'Landsat 8 NDVI time series at ROI',
      trendlines: {0: {
        color: 'CC0000'
      }},
      lineWidth: 1,
      pointSize: 3,
    });
print(l8Chart);

// 3. Linear modeling of time

// A.Estimate linear trend over time

// List of the independent variable names
var independents = ee.List(['constant', 't']);
// Name of the dependent variable.
var dependent = ee.String('NDVI');
// Compute a linear trend.  This will have two bands: 'residuals' and 
// a 2x1 band called coefficients (columns are for dependent variables).
var trend = filteredLandsat.select(independents.add(dependent))
    .reduce(ee.Reducer.linearRegression(independents.length(), 1));
// Map.addLayer(trend, {}, 'trend array image');
// Flatten the coefficients into a 2-band image
var coefficients = trend.select('coefficients')
  .arrayProject([0])
  .arrayFlatten([independents]);

// Compute a de-trended series.
var detrended = filteredLandsat.map(function(image) {
  return image.select(dependent).subtract(
          image.select(independents).multiply(coefficients).reduce('sum'))
          .rename(dependent)
          .copyProperties(image, [timeField]);
});

// Plot the detrended results.
var detrendedChart = ui.Chart.image.series(detrended, roi, null, 30)
    .setOptions({
      title: 'Detrended Landsat time series at ROI',
      lineWidth: 1,
      pointSize: 3,
    });
print(detrendedChart);

// B.Estimate seasonality with a harmonic model

// Use these independent variables in the harmonic regression.
var harmonicIndependents = ee.List(['constant', 't', 'cos', 'sin']);
// Add harmonic terms as new image bands.
var harmonicLandsat = filteredLandsat.map(function(image) {
  var timeRadians = image.select('t').multiply(2 * Math.PI);
  return image
    .addBands(timeRadians.cos().rename('cos'))
    .addBands(timeRadians.sin().rename('sin'));
});

// The output of the regression reduction is a 4x1 array image.
var harmonicTrend = harmonicLandsat
  .select(harmonicIndependents.add(dependent))
  .reduce(ee.Reducer.linearRegression(harmonicIndependents.length(), 1));

// Turn the array image into a multi-band image of coefficients.
var harmonicTrendCoefficients = harmonicTrend.select('coefficients')
  .arrayProject([0])
  .arrayFlatten([harmonicIndependents]);
// Compute fitted values.
var fittedHarmonic = harmonicLandsat.map(function(image) {
  return image.addBands(
    image.select(harmonicIndependents)
      .multiply(harmonicTrendCoefficients)
      .reduce('sum')
      .rename('fitted'));
});
// Plot the fitted model and the original data at the ROI.
print(ui.Chart.image.series(
fittedHarmonic.select(['fitted','NDVI']), roi, ee.Reducer.mean(), 30)
    .setSeriesNames(['NDVI', 'fitted'])
    .setOptions({
      title: 'Harmonic model: original and fitted values',
      lineWidth: 1,
      pointSize: 3,
}));

// Compute phase and amplitude.
var phase = harmonicTrendCoefficients.select('cos').atan2(
            harmonicTrendCoefficients.select('sin'));
            
var amplitude = harmonicTrendCoefficients.select('cos').hypot(
                harmonicTrendCoefficients.select('sin'));
// Use the HSV to RGB transform to display phase and amplitude
var rgb = phase.unitScale(-Math.PI, Math.PI).addBands(
          amplitude.multiply(2.5)).addBands(
          ee.Image(1)).hsvToRgb();
Map.addLayer(rgb, {}, 'phase (hue), amplitude (saturation)');

// C.More on harmonic models

// 4. Autocovariance and autocorrelation

// A.Create a lagged ImageCollection

var lag = function(leftCollection, rightCollection, lagDays) {
  var filter = ee.Filter.and(
    ee.Filter.maxDifference({
      difference: 1000 * 60 * 60 * 24 * lagDays,
      leftField: timeField, 
      rightField: timeField
    }), 
    ee.Filter.greaterThan({
      leftField: timeField, 
      rightField: timeField
  }));
  
  return ee.Join.saveAll({
    matchesKey: 'images',
    measureKey: 'delta_t',
    ordering: timeField, 
    ascending: false, // Sort reverse chronologically
  }).apply({
    primary: leftCollection, 
    secondary: rightCollection, 
    condition: filter
  });
};

var lagged17 = lag(detrended, detrended, 17);

// B.Compute autocovariance and autocorrelation

var merge = function(image) {
  // Function to be passed to iterate.
  var merger = function(current, previous) {
    return ee.Image(previous).addBands(current);
  };
  return ee.ImageCollection.fromImages(
image.get('images')).iterate(merger, image);
};

var merged17 = ee.ImageCollection(lagged17.map(merge));

var covariance = function(mergedCollection, band, lagBand) {
  return mergedCollection.select([band, lagBand]).map(function(image) {
    return image.toArray();
  }).reduce(ee.Reducer.covariance(), 8);
};
var lagBand = dependent.cat('_1');
var covariance17 = ee.Image(covariance(merged17, dependent, lagBand));

Map.addLayer(covariance17.arrayGet([0, 1]), {}, 'covariance (lag=34 days)');

var correlation = function(vcArrayImage) {
  var covariance = ee.Image(vcArrayImage).arrayGet([0, 1]);
  var sd0 = ee.Image(vcArrayImage).arrayGet([0, 0]).sqrt();
  var sd1 = ee.Image(vcArrayImage).arrayGet([1, 1]).sqrt();
  return covariance.divide(sd0).divide(sd1).rename('correlation');
};

var correlation17 = correlation(covariance17);
Map.addLayer(correlation17, {min: -1, max: 1}, 
'correlation (lag = 34 days)');

// 5. Cross-covariance and Cross-correlation

// Precipitation (covariate)
var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/PENTAD');
// Join the t-l (l=1 pentad) precipitation images to the Landsat.
var lag1PrecipNDVI = lag(filteredLandsat, chirps, 5);
// Add the precipitation images as bands.
var merged1PrecipNDVI = ee.ImageCollection(lag1PrecipNDVI.map(merge));
// Compute and display cross-covariance.
var cov1PrecipNDVI = covariance(merged1PrecipNDVI, 'NDVI', 'precipitation');
Map.addLayer(cov1PrecipNDVI.arrayGet([0, 1]), {}, 'NDVI - PRECIP cov (lag = 5)');
// Compute and display cross-correlation.
var corr1PrecipNDVI = correlation(cov1PrecipNDVI);
Map.addLayer(corr1PrecipNDVI, {min: -0.5, max: 0.5}, 'NDVI - PRECIP corr (lag = 5)');

// Join the precipitation images from the previous month
var lag30PrecipNDVI = lag(filteredLandsat, chirps, 30);
print(lag30PrecipNDVI);
var sum30PrecipNDVI = ee.ImageCollection(lag30PrecipNDVI.map(function(image) {
  var laggedImages = ee.ImageCollection.fromImages(image.get('images'));
  return ee.Image(image).addBands(laggedImages.sum().rename('sum'));
}));
// Compute covariance.
var cov30PrecipNDVI = covariance(sum30PrecipNDVI, 'NDVI', 'sum');
Map.addLayer(cov1PrecipNDVI.arrayGet([0, 1]), {}, 'NDVI - sum cov (lag = 30)');
// Correlation.
var corr30PrecipNDVI = correlation(cov30PrecipNDVI);
Map.addLayer(corr30PrecipNDVI, {min: -0.5, max: 0.5}, 'NDVI - sum corr (lag = 30)');

// 6. Auto-regressive models

var lagged34 = ee.ImageCollection(lag(filteredLandsat, filteredLandsat, 34));

var merged34 = lagged34.map(merge).map(function(image) {
  return image.set('n', ee.List(image.get('images')).length());
}).filter(ee.Filter.gt('n', 1));

var arIndependents = ee.List(['constant', 'NDVI_1', 'NDVI_2']);
var ar2 = merged34
  .select(arIndependents.add(dependent))
  .reduce(ee.Reducer.linearRegression(arIndependents.length(), 1));
// Turn the array image into a multi-band image of coefficients.
var arCoefficients = ar2.select('coefficients')
  .arrayProject([0])
  .arrayFlatten([arIndependents]);

// Compute fitted values.
var fittedAR = merged34.map(function(image) {
  return image.addBands(
    image.expression('beta0 + beta1 * p1 + beta2 * p2', {
      p1: image.select('NDVI_1'),
      p2: image.select('NDVI_2'),
      beta0: arCoefficients.select('constant'),
      beta1: arCoefficients.select('NDVI_1'),
      beta2: arCoefficients.select('NDVI_2')
    }).rename('fitted'));
});

print(ui.Chart.image.series(
fittedAR.select(['fitted', 'NDVI']), roi, ee.Reducer.mean(), 30)
    .setSeriesNames(['NDVI', 'fitted'])
    .setOptions({
      title: 'AR(2) model: original and fitted values',
      lineWidth: 1,
      pointSize: 3,
}));


