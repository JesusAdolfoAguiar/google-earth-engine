/*
Charts (Exercise)
Analysis of annual trends in the SOL (sum of lights)
(Image collection: “DMSP-OLS Nighttime Lights Time Series Version 4”)

Please analyze annual trends in the SOL in 2 cities: San Diego and Tijuana, in the
period between 2000 and 2016.

First, create a point for each city. Then, create a buffer of 10Km around each point.

Create a graph that presents the annual trends in the SOL within these buffers.
*/

var DMSP = ee.ImageCollection('NOAA/DMSP-OLS/NIGHTTIME_LIGHTS')
              .filterDate('2000-01-01', '2016-12-31');

var SanFran = ee.Geometry.Point(-117.14489936828613, 32.72801412622118).buffer(10000);
    
var Tijuana = ee.Geometry.Point(-116.99443817138672, 32.56070352232516).buffer(10000);

var chartSanFran = {
title: 'SOL in San Francisco 2002-2016',
hAxis: {title: 'SOL'},
vAxis: {title: 'Years'},
lineWidth: 1,
pointSize: 4,
series: {
0: {color: 'd63000'}, // built-up
1: {color: '98ff00'}, // green
2: {color: '8b823b'}, // bare land
  }
};

var chartTijuana = {
title: 'SOL in San Tijuana 2002-2016',
hAxis: {title: 'SOL'},
vAxis: {title: 'Years'},
lineWidth: 1,
pointSize: 4,
series: {
0: {color: 'd63000'}, // built-up
1: {color: '98ff00'}, // green
2: {color: '8b823b'}, // bare land
  }
};
print(ui.Chart.image.series(DMSP, SanFran, ee.Reducer.mean()).setOptions(chartSanFran));

print(ui.Chart.image.series(DMSP, Tijuana, ee.Reducer.mean()).setOptions(chartTijuana));

Map.addLayer(SanFran);
Map.addLayer(Tijuana);
