
//Challenge: Navigate to any region and present 2 Landsat images for the rainy season of
//2013 and one for 2016, using the median value of the composite

var amazon = /* color: #d63000 */ee.Geometry.Point([-64.40322875976562, 6.685061972567749]);

var landsat8TOA = ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA')
                  .filterBounds(amazon);  
print(landsat8TOA);

var landsat8TOA2013 = landsat8TOA
                      .filterDate('2013-06-01', '2013-09-30');
print(landsat8TOA2013);

var landsat8TOA2013AMedian = landsat8TOA2013.median();

Map.addLayer(landsat8TOA2013AMedian, {'bands': ['B4', 'B5', 'B3']});

print(landsat8TOA2013AMedian);

var image = ee.Algorithms.Landsat.simpleComposite({
  collection: landsat8TOA2013AMedian,
  asFloat: true,
  cloudScoreRange: 80
});

print(image);

//Just repeat the above code with the other period of interest to get the other inmages
//End