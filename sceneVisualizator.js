Map.setCenter(-63.5000000, 4, 6);

// Setup parameters

//************ PARAMETERS FOR IMAGE VISUALIZATION *****************
var vpars = {
  'year'    : 2016, 
  'path'    : 1,
  'row'     : 59, 
  'clInit'  : 0,
  'clStop'  : 30
}
//*****************************************************************


// import GEE Landsat collection depending of year
var coll_id;
if(vpars.year <= 2010) {
  coll_id = 'LT5_L1T_TOA'
} else if(vpars.year < 2013){
  coll_id = 'LE7_L1T_TOA'
} else {
  coll_id = 'LC8_L1T_TOA'
}
print(
  "LANDSAT/" + coll_id + " collection selected", 
  vpars.year, 'cloud ' + vpars.clInit + '-' + vpars.clStop,
  'scene ' + vpars.path + '-' + vpars.row
)
var coll = ee.ImageCollection('LANDSAT/' + coll_id);

// attempt to remove clouds.
var maskClouds = function(raster) {
  var scored  = ee.Algorithms.Landsat.simpleCloudScore(raster);
  var rasters = scored.select(['cloud']); 
  return raster.updateMask(rasters.lt(vpars.clStop));
};

//filter by date and cloud cover
var fmyr = ee.Number(vpars.year).format()
var coll_fil = coll.filterDate(fmyr.cat("-01-01"), fmyr.cat("-12-31"))
  .filter(ee.Filter.gt('CLOUD_COVER', vpars.clInit))
  .filter(ee.Filter.lt('CLOUD_COVER', vpars.clStop))
  .filter(ee.Filter.eq('WRS_PATH', vpars.path))
  .filter(ee.Filter.eq('WRS_ROW', vpars.row))
  .map(maskClouds)


var size = coll_fil.size().getInfo();
print(size + ' images found and added to map');
//print('Image in ' + vpars.path + '-' + vpars.row + ' ' + ' added to map')

if(size !== 0){
  var listImages = coll_fil.toList(size);
  // Define the visualization parameters.
  var vizParams = {
    bands: ['B4', 'B3', 'B2'],
    min: 0, max: 0.5,
    gamma: [0.95, 1.1, 1]
  };
  
  //Display layers
  for (var i = 0; i < size; i++){
    var data = ee.Image(listImages.get(i));
    var id  = data.id().getInfo();
    Map.addLayer(data, vizParams, id);
  }
}
