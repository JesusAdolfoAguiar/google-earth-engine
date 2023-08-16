
var startyear = 2000; 
var endyear = 2016; 

var startmonth = 1;
var endmonth = 1;

var startdate = ee.Date.fromYMD(startyear, startmonth, 1);
var enddate = ee.Date.fromYMD(endyear + 1, endmonth, 1);

var years = ee.List.sequence(startyear, endyear);


var redriver = ee.FeatureCollection(
    'ft:1aA9E2E6SXMWI5PwWtO6esGYCoA2ZfUJ151U8g-eO','geometry');


var annualPrecip = chirps.filterDate(startdate, enddate)
  .sort('system:time_start', false)
  .filterBounds(redriver);


var pViz = {
  min: 0, 
  max: 2400, 
  palette: '000000, 0000FF, FDFF92, FF2700, FF00E7'
};

var title = {
  title: 'Annual precipitation',
  hAxis: {title: 'Time'},
  vAxis: {title: 'Precipitation (mm)'},
};


var chart = ui.Chart.image.seriesByRegion({
  imageCollection: annualPrecip, 
  regions: redriver,
  reducer: ee.Reducer.mean(),
  band: 'precipitation',
  scale: 2500,
  xProperty: 'system:time_start',
  seriesProperty: 'SITE'
}).setOptions(title)
  .setChartType('ColumnChart');


print(chart);

var annualMean = annualPrecip.mean().clip(redriver);
Map.centerObject(redriver, 7);
Map.addLayer(annualMean, pViz, 'mean yearly P');


var months = ee.List.sequence(1, 12);

var monthlyPrecip =  ee.ImageCollection.fromImages(
  years.map(function (y) {
    return months.map(function(m) {
      var w = annualPrecip.filter(ee.Filter.calendarRange(y, y, 'year'))
                    .filter(ee.Filter.calendarRange(m, m, 'month'))
                    .sum();
      return w.set('year', y)
              .set('month', m)
              .set('system:time_start', ee.Date.fromYMD(y, m, 1));
                        
    });
  }).flatten()
);

var monthlyMean = monthlyPrecip.mean().clip(redriver);
pViz.max = 300;
Map.addLayer(monthlyMean, pViz, 'mean monthly P');


var meanMonthlyP =  ee.ImageCollection.fromImages(
  months.map(function (m) {
    var w = monthlyPrecip.filter(ee.Filter.eq('month', m)).mean();
    return w.set('month', m)
            .set('system:time_start',ee.Date.fromYMD(1, m, 1)); 
  }).flatten()
);

title.title= 'Monthly precipitation';
var chartMonthly = ui.Chart.image.seriesByRegion({
  imageCollection: meanMonthlyP, 
  regions: redriver,
  reducer: ee.Reducer.mean(),
  band: 'precipitation',
  scale: 2500,
  xProperty: 'system:time_start',
  seriesProperty: 'SITE'
}).setOptions(title)
  .setChartType('ColumnChart');

print(chartMonthly);

var et = ee.ImageCollection('users/atepoortinga/etdata');



print(et);

var etViz = {
  min: 0, 
  max: 150, 
  palette: '000000,0000FF,FDFF92,FF2700,FF00E7'
};
Map.addLayer(et.mean(), etViz, 'Mean evapotranspiration');


var meanMonthlyET = months.map(function(m) {
    return et.filter(ee.Filter.calendarRange(m, m, 'month'))
        .mean()
        .rename('et')
        .set('month', m)
        .set('system:time_start', ee.Date.fromYMD(1, m, 1)); 
});

title.title = 'Monthly mean evapotranspiration';
title.vAxis.title = 'Evapotranspiration (mm)';
var monthlyETChart = ui.Chart.image.seriesByRegion({
  imageCollection: meanMonthlyET, 
  regions: redriver,
  reducer: ee.Reducer.mean(),
  band: 'et',
  scale: 2500,
  xProperty: 'system:time_start',
  seriesProperty: 'SITE'
}).setOptions(title)
  .setChartType('ColumnChart');

print(monthlyETChart);


var joined = ee.Join.saveAll({
  matchesKey: 'etImages'
}).apply({
  primary: monthlyPrecip, 
  secondary: et.map(function(image) {
    return image.set('month', image.date().get('month'))
                .set('year', image.date().get('year'));
  }), 
  condition: ee.Filter([
    ee.Filter.equals({
      leftField: 'month', 
      rightField: 'month'
    }),
    ee.Filter.equals({
      leftField: 'year', 
      rightField: 'year'
    })
  ])
});


var effectiveRainfall = joined.map(function(image) {
  var et = ee.ImageCollection.fromImages(image.get('etImages'))
      .sum();
  return ee.Image(image).subtract(et)
      .rename('effective_rainfall')
      .copyProperties(image, ['system:time_start']);
});


var meanER = ee.ImageCollection(effectiveRainfall).mean();
Map.addLayer(meanER, {}, 'mean effective rainfall');

title['title'] = 'Effective Rainfall';
title['vAxis']['title'] = 'Effective Rainfall (mm)';
var erChart = ui.Chart.image.seriesByRegion({
  imageCollection: effectiveRainfall, 
  regions: redriver,
  reducer: ee.Reducer.mean(),
  band: 'effective_rainfall',
  scale: 2500,
  xProperty: 'system:time_start',
  seriesProperty: 'SITE'
}).setOptions(title)
  .setChartType('ColumnChart');
print(erChart);





