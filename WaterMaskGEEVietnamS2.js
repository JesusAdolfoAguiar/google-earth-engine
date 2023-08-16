
var countryNames = [
  'Myanmar (Burma)', 'Thailand', 'Laos', 'Vietnam', 'Cambodia'
]; 

var countries = ee.FeatureCollection(
  'ft:1tdSwUL7MVpOauSgRzqVTOwdfy17KDbw-1d9omPw');

var mekongCountries = countries.filter(
  ee.Filter.inList('Country', countryNames));

var mekongRegion = mekongCountries.geometry();

var computeQAbits = function(image, start, end, newName) {
    var pattern = 0;

    for (var i=start; i<=end; i++) {
        pattern += Math.pow(2, i);
    }

    return image.select([0], [newName]).bitwiseAnd(pattern).rightShift(start);
};

var sentinel2 = function(image) {

  var cloud_mask = image.select("QA60");
  var opaque = computeQAbits(cloud_mask, 10, 10, "opaque");
  var cirrus = computeQAbits(cloud_mask, 11, 11, "cirrus");
  var mask = opaque.or(cirrus);

 return image.updateMask(mask.not());
};

var s2mask = require('users/fitoprincipe/geetools:cloud_masks').sentinel2;

var AOI = ee.Geometry.Rectangle(-85.14404296875, 46.9502622421856, 
  -71.87255859375, 41.60722821271717);

var collection = ee.ImageCollection("COPERNICUS/S2")
  .filterdate('2015-01-01', '2016-12-31')
  .filterBounds(mekongRegion)
  .filterMetadata('CLOUDY_PIXEL_PERCENTAGE', 'less_than', 70)
  .map(s2mask);

var image2 = ee.Image(collection.mean());
Map.addLayer(image, {bands:['B8', 'B12', 'B4'], min:0, max:3000});
Map.centerObject(mekongRegion);
