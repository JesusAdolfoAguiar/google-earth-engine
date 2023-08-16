
var startdate = ee.Date('2015-01-01');
var enddate = ee.Date('2016-12-31');


var countryNames = [
  'Myanmar (Burma)', 'Thailand', 'Laos', 'Vietnam', 'Cambodia'
]; 

var countries = ee.FeatureCollection(
  'ft:1tdSwUL7MVpOauSgRzqVTOwdfy17KDbw-1d9omPw');

var mekongCountries = countries.filter(
  ee.Filter.inList('Country', countryNames));

var mekongRegion = mekongCountries.geometry();


var l8images = L8.filterDate(startdate, enddate)
    .filterBounds(mekongRegion);


var cloudThresh = 20;

var cloudfunction = function(image){
  var scored = ee.Algorithms.Landsat.simpleCloudScore(image);
  var score = scored.select('cloud');
  var cloudy = quality.gt(cloudThresh);
  var cloudmask = cloudy.not();
  return image.updateMask(cloudmask);
};


var l8CloudMasked = l8images.map(cloudFunction);


Map.centerObject(mekongRegion, 5);
Map.addLayer(l8CloudMasked.median().clip(mekongRegion), {
    min: 0, 
    max: 0.5,
    bands: ['B4', 'B3', 'B2']
  }, 'Landsat 8 True color');


function addNdwi(img) {
  var ndwi = img.normalizedDifference(['green', 'nir']).rename('NDWI');
  return img.addBands(ndwi);
}


var bands = ['green', 'nir'];

var l8ndwi = l8CloudMasked
    .select(['green', 'nir'], bands)
    .map(addNdwi);


print(l8ndwi.first());


var ndwiViz = {bands: 'NDWI', min: 0.0, max: 0.3, palette: '0000FF'};

var ndwimean = l8ndwi.select('NDWI').mean();

Map.addLayer(ndwimean.mask(ndwimean), ndwiViz, 'NDWI');


var l5ndwi = L5
  .filterBounds(mekongRegion)
  .filterDate('1984-01-01', '1998-12-31')
  .map(cloudFunction)
  .select(['B2', 'B4'], bands)
  .map(addNdwi);

var l7ndwi = L7
  .filterBounds(mekongRegion)
  .filterDate('1999-01-01', '2012-12-31')
  .select(['B2', 'B4'], bands)
  .map(addNdwi);

var collection = ee.ImageCollection(
  l5ndwi.merge(l7ndwi).merge(l8ndwi));


var year = 1990;

var ndwi1990 = collection.filterDate({
  start: ee.Date.fromYMD(year, 1, 1),
  end: ee.Date.fromYMD(year, 12, 31)
}).max();

Map.addLayer(ndwi1990.mask(ndwi1990), ndwiViz, '1990 NDWI');


var waterBinary = collection.select('NDWI').map(function(image) {
  return image.gt(0);
});

var frequency = waterBinary.sum().divide(waterBinary.count());
Map.addLayer(frequency.mask(frequency), 
  {palette: ['white', 'magenta', 'blue']}, 'inundation frequency');

Export.image.toDrive(frequency);


