var catchments = ee.FeatureCollection('projects/ee-ramira24/assets/yachats');
Map.addLayer(catchments, {color: 'FF0000'}, 'Sample Site Catchments');

// Visualization parameters for NDVI
var ndviVis = {
  min: -1,
  max: 1,
  palette: ['blue', 'white', 'green']
};

for (var i = 2013; i < 2023; i++) {
  var dataset = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterDate(i + '-01-01', (i + 1) + '-01-01')
    .filterBounds(catchments)
    .map(function(image) {
      var nir = image.select('SR_B5').multiply(0.0000275).add(-0.2); // NIR
      var red = image.select('SR_B4').multiply(0.0000275).add(-0.2); // Red
      var ndvi = nir.subtract(red).divide(nir.add(red)).rename('NDVI');
      return ndvi.copyProperties(image, ["system:time_start"]);
    });

  // Median composite of NDVI for the year
  var annualNDVI = dataset.median().clip(catchments).toFloat();

  Map.addLayer(annualNDVI, ndviVis, 'NDVI ' + i);

  Export.image.toDrive({
    image: annualNDVI,
    description: 'Landsat8_NDVI_' + i,
    scale: 30,
    region: catchments.geometry().bounds(),
    maxPixels: 1e13
  });
}

Map.centerObject(catchments, 10);
