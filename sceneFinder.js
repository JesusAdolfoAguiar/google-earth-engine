Map.setCenter(-63.5000000, 4, 6);

// Setup parameters

//************ PARAMETERS FOR IMAGE VERIFICATION ******************
var params = {
  'year'    : 2017,
  'path'    : [1, 2, 3, 4, 232, 233],
  'row'     : [53, 54, 55, 56, 57, 58, 59, 60], 
  'clInit'  : 90,
  'clStop'  : 100
}
//*****************************************************************

// import GEE Landsat collection depending of year
var coll_id;
if(params.year <= 2010) {
  coll_id = 'LT5_L1T_TOA'
} else if(params.year < 2013){
  coll_id = 'LE7_L1T_TOA'
} else {
  coll_id = 'LC8_L1T_TOA'
}
print(
  "LANDSAT/" + coll_id + " collection selected", 
  params.year, 'cloud ' + params.clInit + '-' + params.clStop
)
var coll = ee.ImageCollection('LANDSAT/' + coll_id);

// attempt to remove clouds.
var maskClouds = function(raster) {
  var scored  = ee.Algorithms.Landsat.simpleCloudScore(raster);
  var rasters = scored.select(['cloud']); 
  return raster.updateMask(rasters.lt(params.clStop));
};

//filter by date and cloud cover
var fmyr = ee.Number(params.year).format()
var coll_fil = coll.filterDate(fmyr.cat("-01-01"), fmyr.cat("-12-31"))
  .filter(ee.Filter.gt('CLOUD_COVER', params.clInit))
  .filter(ee.Filter.lt('CLOUD_COVER', params.clStop))
  .map(maskClouds)


//************ FUNCTION FOR IMAGE VERIFICATION ******************


var vizbyPathRow = function(collection, path, row) {
  for(var k = 0; k < path.length; k++) {
    for(var l = 0; l < row.length; l++) {
      if(
        // this condition is for display only the scenes
        // of venezuelan amazonia
        path[k] == 1 && row[l] < 60 || path[k] == 2 || 
        path[k] == 3 && row[l] > 53 && row[l] < 60 || 
        path[k] == 4 && row[l] > 54 && row[l] < 59 ||
        path[k] == 232 && row[l] > 53 && row[l] < 58 ||
        path[k] == 233 && row[l] > 52 && row[l] < 60
      ){
        //print(params.path[k] + '-' + params.row[l])
        coll_filt = collection
          .filter(ee.Filter.eq('WRS_PATH', path[k]))
          .filter(ee.Filter.eq('WRS_ROW', row[l]))
  
        var dates = ee.List(coll_filt.get('date_range'));
        var dateRange = ee.DateRange(dates.get(0), dates.get(1));
  
        var size = coll_filt.size().getInfo();
        var output = {
          'path-row' : path[k] + '-' + row[l],
          'avl_date' : dateRange,
          'avl_data' : size
        };
        print(path[k] + ',' + row[l] + ',' + output.avl_data);
      }
    }
  }
};

// Run verification
var coll_filt = vizbyPathRow(coll_fil, params.path, params.row) 


//**************************************************************