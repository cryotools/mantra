/****************************************************

A tool to evaluate glacier surface classification and Transient Snowline Altitude (TSLA) calculation results from the 
"MountAiN glacier Transient snowline Retrieval Algorithm" (MANTRA).

Please note that results provided by this tool may be flawed, e.g. by cloud cover. 
It is advised to check the results by comparing the classification with the Landsat image by toggling layers before.

FEATURES
- Classify surface materials based on Landsat band-ratios (snow, ice, debris/rock, clouds)
- Create a catalogue of relevant Landsat scenes from multiple missions (Landsat 4 to 8)
- Merge scenes taken at the same date that are covered by the glacier's outline
- Mask pixels with insufficient illumination
- Visualize the classification
- Calculate and display statistics


ACKNOWLEDGEMENTS:
Routine to merge multiple Landsat sensors into one collection intensively inspired by Matthias Baumann's (HU Berlin) scripts.
Inge Grünberg (AWI Potsdam) developed the two ratios additional to the NDSI and helped finding adequate thresholds.


2022, David Loibl

*****************************************************/

/* CONFIGURATION AND IMPORTS*/

var rgi_id = 'RGI60-13.53885';                  // RGI id of glacier to analyze 
var start = ee.Date('2019-09-09');              // First date to include in the search for Landsat scenes
var finish = ee.Date('2019-09-10');             // Last date to include in the search for Landsat scenes

var tsl_median_binsize = 2;                     // Percentile of altitude range to obtain TSL from
var cf_threshold = 5;                           // Minimum cloudfree portion (%)of glacier surface to conduct TSL analysis

var rgi_path = 'users/___your_GEE_account___/rgi60_Asia_combined'; // Path to RGI shapefile containing the glacier to investigate
var core = require('users/___your_GEE_account___/tsl-retrieval:core');

var dem = 'JAXA/ALOS/AW3D30_V1_1';


/* LOAD DATASETS */

var rgi6_fc = ee.FeatureCollection(rgi_path);
var glaciers = rgi6_fc.filter(ee.Filter.eq('RGIId', rgi_id));

var selected_glacier = ee.Feature(glaciers.first());
var glacier_area = selected_glacier.area(0.1).divide(1000 * 1000);
Map.centerObject(selected_glacier);   

var alos_DEM = ee.Image(dem);
var alos_DEM = alos_DEM.select('AVE');

//Build a joined Landsat 4-8 collection ()
var bands         = ['B1', 'B2', 'B3', 'B4',  'B5',    'B7',    'B6',         'BQA'];
var band_names    = ['B',  'G',  'R',  'NIR', 'SWIR1', 'SWIR2', 'T',          'BQA'];
var l7bands       = ['B1', 'B2', 'B3', 'B4',  'B5',    'B7',    'B6_VCID_2',  'BQA'];
var l7band_names  = ['B',  'G',  'R',  'NIR', 'SWIR1', 'SWIR2', 'T',          'BQA'];
var l8bands       = ['B2', 'B3', 'B4', 'B5',  'B6',    'B7',    'B10',        'BQA'];
var l8band_names  = ['B',  'G',  'R',  'NIR', 'SWIR1', 'SWIR2', 'T',          'BQA'];

var ls4 = ee.ImageCollection('LANDSAT/LT04/C01/T1_TOA').map(function(image){return image.select(bands).rename(band_names)});
var ls5 = ee.ImageCollection('LANDSAT/LT05/C01/T1_TOA').map(function(image){return image.select(bands).rename(band_names)});
var ls7 = ee.ImageCollection('LANDSAT/LE07/C01/T1_TOA').map(function(image){return image.select(l7bands).rename(l7band_names)}).filterMetadata('GAIN_BAND_5', 'equals', 'L');  // Remove scenes where TOA was processed with high gain
var ls8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA').map(function(image){return image.select(l8bands).rename(l8band_names)});

var landsatCollection = ee.ImageCollection(ls4.merge(ls5).merge(ls7).merge(ls8))
                  
                  .filterDate(start, finish)  
                  .filterBounds(selected_glacier.geometry());
                  

// print('Glacier', selected_glacier);               
print ('landsatCollection', landsatCollection);




/* MOSAIC SATELLITE IMAGES BY DATE */

// Difference in days between start and finish
var diff = finish.difference(start, 'day');
var temporalResolution = 1;  // days

// Make a list of all dates
var range = ee.List.sequence(0, diff.subtract(1), temporalResolution).map(function(day){return start.advance(day,'day')});


// Funtion for iteraton over the range of dates
var day_list = function(date, newlist) {
  // Cast
  date = ee.Date(date);
  newlist = ee.List(newlist);

  // Filter collection between date and the next day
  var filtered = landsatCollection.filterDate(date, date.advance(temporalResolution, 'day'));
  var LS_DATE = filtered.first().get('DATE_ACQUIRED'); //date.format('YYYY-MM-dd');

  return ee.List(ee.Algorithms.If(filtered.size(), newlist.add(LS_DATE), newlist));
};

var sceneList = ee.List(range.iterate(day_list, ee.List([])));

// Check if multiple images exist for individual days. If so, create mosaics.
var get_day_images = function(day) {
  var filtered = landsatCollection.filterDate(ee.Date(day), ee.Date(day).advance(1, 'day'));
  var ls_id = filtered.first().get('LANDSAT_PRODUCT_ID');
  var sensor = filtered.first().get('SPACECRAFT_ID');
  var sun_azi = filtered.first().get('SUN_AZIMUTH');
  var sun_ele = filtered.first().get('SUN_ELEVATION');
  var image = ee.Algorithms.If(filtered.size().gt(1), 
                  ee.Image(filtered.mosaic().clip(selected_glacier))
                          .set('LS_DATE', day)
                          .set('LS_ID', ls_id)
                          .set('SENSOR', sensor)
                          .set('SUN_AZIMUTH', sun_azi)
                          .set('SUN_ELEVATION', sun_ele),
                  filtered.first().clip(selected_glacier)
                          .set('LS_DATE', day)
                          .set('LS_ID', ls_id)
                          .set('SENSOR', sensor)
                          .set('SUN_AZIMUTH', sun_azi)
                          .set('SUN_ELEVATION', sun_ele)
                  );
  return image;
};



/* FUNCTIONS */


var obtain_DEM_stats = function(vector_outline) {
  var alos_clipped = alos_DEM.clip(vector_outline);
  
  // Transient snowline estimation (-> median of lowest x% elevation of snow-covered region)
  var lowest_percentiles = alos_clipped.reduceRegion({reducer: ee.Reducer.percentile([tsl_median_binsize]),
                                geometry: vector_outline.geometry(),
                                bestEffort: true}).get('AVE');
  
  var masked_DEM = alos_clipped.updateMask(alos_clipped.lte(ee.Number(lowest_percentiles)));
  
  var reducer1 = ee.Reducer.mean();
  var reducers = reducer1.combine({reducer2: ee.Reducer.median(), sharedInputs: true})
                         .combine({reducer2: ee.Reducer.max(), sharedInputs: true})
                         .combine({reducer2: ee.Reducer.min(), sharedInputs: true})
                         .combine({reducer2: ee.Reducer.stdDev(), sharedInputs: true});

  var low_perc_stats = masked_DEM.reduceRegion({reducer: reducers,
                                  geometry: vector_outline.geometry(),
                                  bestEffort: true});

  // Put key information into a dictionary and create a feature from it 
  var dict = {
    RGIId: vector_outline.get('RGIId'), 
    LS_ID: vector_outline.get('LS_ID'),
    LS_DATE: vector_outline.get('LS_DATE'),
    // CLOUDFREE_AREA: lsScene.get('CLOUDFREE_AREA'),
    GLACIER_AREA: vector_outline.get('GLACIER_AREA'),
    STATS_MEAN: low_perc_stats.get('AVE_mean'),
    STATS_MEDIAN: low_perc_stats.get('AVE_median'),
    STATS_MIN: low_perc_stats.get('AVE_min'),
    STATS_MAX: low_perc_stats.get('AVE_max'),
    STATS_STDEV: low_perc_stats.get('AVE_stdDev')
    
  };
  var DEM_stats_result = ee.Feature(vector_outline.geometry(), dict); //dict);
  
  return DEM_stats_result;

};


var create_nan_dict = function(vector_outline) {
  var dict = {
    RGIId: vector_outline.get('RGIId'), 
    LS_ID: vector_outline.get('LS_ID'),
    LS_DATE: vector_outline.get('LS_DATE'),
    CLOUDFREE_AREA: "Less "+ cf_threshold +"%",
    GLACIER_AREA: vector_outline.get('GLACIER_AREA'),
    STATS_MEAN: "NaN",
    STATS_MEDIAN: "NaN",
    STATS_MIN: "NaN",
    STATS_MAX: "NaN",
    STATS_STDEV: "NaN"
  };
  var nan_result = ee.Feature(vector_outline.geometry(), dict);
  return nan_result;
};



/* USER INTERFACE */

var sceneList = ee.List(range.iterate(day_list, ee.List([])));

// Create a panel with vertical flow layout.
var panel = ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: {width: '300px'}
});

var removeLayer = function(name) {
  var layers = Map.layers();
  // list of layers names
  var names = [];
  layers.forEach(function(lay) {
    var lay_name = lay.getName();
    names.push(lay_name);
  });
  // get index
  var index = names.indexOf(name);
  if (index > -1) {
    // if name in names
    var layer = layers.get(index);
    Map.remove(layer);
  }
};


var sat_scene_dict = sceneList.map(function(scene) {
  var feature = selected_glacier.set('system:index', scene)
                                .set('value', scene)
                                .set('label', scene);
  return feature;
});


var features = sat_scene_dict.getInfo();
var select_items = [];
for (var i = 0; i < features.length; i++) {
  select_items.push({
    label: features[i]['properties']['label'],
    value: features[i]['properties']['value']
  });
}

var cleanup_map = function(scene_id){
  removeLayer('LS ' + scene_id);
  removeLayer('Snow cover raster ' + scene_id);
  removeLayer('Ice cover raster ' + scene_id);
  removeLayer('Cloud cover raster ' + scene_id);
  removeLayer('Debris cover raster ' + scene_id);
  removeLayer('CC TSL Range ' + scene_id);
  removeLayer('Masked DEM ' + scene_id);
};          

var last_scene_displayed = null;

var imageSelect = ui.Select({
  items: select_items,
  onChange: function(value) {
    
    // Conduct snow cover and TSL analysis for selected date

    var selected_scene = ee.Image(get_day_images(value));
    
    // Write some image parameters to vars
    var scene_date  = ee.Number(selected_scene.get('LS_DATE'));
    var scene_id    = ee.Number(selected_scene.get('LS_ID'));    
    

    ee.Algorithms.If(ee.Algorithms.IsEqual(last_scene_displayed, null),
      null,
      cleanup_map(last_scene_displayed)
    );
    
    last_scene_displayed = value;
    
    var scene_outline = selected_scene.geometry();
    
    // Mask shaded areas using DEM-based illumination calculation
    selected_scene = core.calculate_illumination(selected_scene);
    //selected_scene = selected_scene.updateMask(selected_scene.select('T').gt(0)); 
    
    // Apply band ratio to classify different surface cover types
    var snow_cover = ee.Image(core.snow_identification(selected_scene))
      .set('RGIId', rgi_id)
      .set('GLACIER_AREA', glacier_area);
    var SC_masked = snow_cover.updateMask(snow_cover.eq(1));    
    
    var ice_cover = ee.Image(core.ice_identification(selected_scene))
      .set('RGIId', rgi_id)
      .set('GLACIER_AREA', glacier_area);    
    var IC_masked = ice_cover.updateMask(ice_cover.eq(1));
    
    var debris_cover = ee.Image(core.debris_identification(selected_scene))
      .set('RGIId', rgi_id)
      .set('GLACIER_AREA', glacier_area);
    var DC_masked = debris_cover.updateMask(debris_cover.eq(1));
        
    var debris_plus_ice = debris_cover.add(ice_cover);


    // Make cloud cover / unknown dataset as all pixels not belonging to the classes identfied before 
    var detected_areas = debris_plus_ice.add(snow_cover);
    var area_false = detected_areas.gt(0);
    var area_true = detected_areas.eq(0);
    var cloud_cover = ee.Image.cat([area_true, area_false]).select(
    ['SWIR1', 'SWIR1_1'], ['covered', 'free']).select('covered');
    //boolraster_zones = boolraster_zones.updateMask(boolraster_zones.neq(0));
    
    
    var total_class_coverage = detected_areas.add(cloud_cover);
    total_class_coverage = total_class_coverage.updateMask(total_class_coverage.neq(0));
    
    cloud_cover = cloud_cover.updateMask(cloud_cover.neq(0));
    
    var CC_masked = cloud_cover.updateMask(cloud_cover.gte(1));
    

	// Create vector outlines from classified rasters
    var snow_cover_outline        = core.boolraster_to_outline(snow_cover);
    var ice_cover_outline         = core.boolraster_to_outline(ice_cover);
    var debris_cover_outline      = core.boolraster_to_outline(debris_cover);    
    var cloud_cover_outline       = core.boolraster_to_outline(cloud_cover);
    var debris_ice_cover_outline  = core.boolraster_to_outline(debris_plus_ice);
    var total_class_cov_outline   = core.boolraster_to_outline(total_class_coverage);
    
    
    // Calculate spatial extents of different surface classes
    var IC_area         = ee.Number(ice_cover_outline.area(0.1).divide(1000 * 1000));
    var DC_area         = ee.Number(debris_cover_outline.area(0.1).divide(1000 * 1000));    
    var CC_area         = ee.Number(cloud_cover_outline.area(0.1).divide(1000 * 1000));
    var SC_area         = ee.Number(snow_cover_outline.area(0.1).divide(1000 * 1000));
    var DIC_area        = ee.Number(debris_ice_cover_outline.area(0.1).divide(1000 * 1000));
    var TCC_area        = ee.Number(total_class_cov_outline.area(0.1).divide(1000 * 1000));
    var TCC_area_ctrl   = ee.Number(IC_area.add(DC_area).add(SC_area).add(CC_area));
    
    
    // For LS7 SLC off scenes: calculate coverage of total glacier (%)
    var class_coverage  = ee.Number(TCC_area.divide(glacier_area).multiply(100));
    var class_cov_ctrl  = ee.Number(TCC_area.divide(glacier_area).multiply(100));
    
    // Calculate the portion of classfied surfaces areas that is cloud-covered (%)
    var CC_total_port   = ee.Number(CC_area.divide(TCC_area).multiply(100));
    
    
    // Obtain DEM stats for snow-covered and ice-plus-debris-covered regions
    snow_cover_outline  = ee.Feature(ee.Algorithms.If(SC_area.gt(0),
                            obtain_DEM_stats(snow_cover_outline),
                            create_nan_dict(snow_cover_outline)
                        ));
    
    /*
    debris_ice_cover_outline = ee.Feature(ee.Algorithms.If(DIC_area.gt(0),
                            obtain_DEM_max(debris_ice_cover_outline),
                            create_nan_dict(debris_ice_cover_outline)
                        ));    
    
    
    var alos_DIC        = alos_DEM.clip(debris_ice_cover_outline);
    var DIC_min_max     = core.get_DEM_min_max(alos_DIC);
    //var glacier_DEM_min = alos_min_max.get('AVE_min');
    var DIC_DEM_max     = DIC_min_max.get('AVE_max');
    
    */
    
    
    
    var obtain_DIC_max = function(debris_ice_cover_outline) {
      var alos_DIC        = alos_DEM.clip(debris_ice_cover_outline);
      var DIC_min_max     = core.get_DEM_min_max(alos_DIC);
      //var glacier_DEM_min = alos_min_max.get('AVE_min');
      var DIC_max         = DIC_min_max.get('AVE_max');
      return DIC_max;
    };
    
    var DIC_max = ee.Number(ee.Algorithms.If(DIC_area.gt(0),
                            obtain_DIC_max(debris_ice_cover_outline),
                            'NaN'
                        )); 
    
    var SC_stdev      = ee.Number(snow_cover_outline.get('STATS_STDEV'));
    var SC_mean       = ee.Number(snow_cover_outline.get('STATS_MEAN'));
    var SC_median     = ee.Number(snow_cover_outline.get('STATS_MEDIAN'));
    
    
    // Identify cloud-covered regions close to TSL (mean +/- 1sigma)
    var alos_glacier  = alos_DEM.clip(selected_glacier);
    var alos_min_max  = core.get_DEM_min_max(alos_glacier);
    var glacier_DEM_min = alos_min_max.get('AVE_min');
    var glacier_DEM_max = alos_min_max.get('AVE_max');
    

    var get_DEM_TSLrange_outline = function(snow_cover_outline) {
      var TSL_range_min = ee.Number(snow_cover_outline.get('STATS_MEAN')).subtract(snow_cover_outline.get('STATS_STDEV'));
      var TSL_range_max = ee.Number(snow_cover_outline.get('STATS_MEAN')).add(snow_cover_outline.get('STATS_STDEV'));

      var DEM_TSL_range = alos_glacier.updateMask(alos_glacier.lte(TSL_range_max).and(alos_glacier.gt(TSL_range_min)));           
      
      var DEM_TSLrange_outline = DEM_TSL_range.reduceToVectors({
        reducer: ee.Reducer.countEvery(), 
        geometryType: 'polygon',
        labelProperty: 'title',
        scale: 30,
        bestEffort: true
      });
      
      var dict = {
        RGIId: rgi_id, 
        LS_ID: scene_id,
        LS_DATE: scene_date
      };
    
      DEM_TSLrange_outline = ee.Feature(DEM_TSLrange_outline.geometry(), dict);
    
      return DEM_TSLrange_outline;
    };
    
    //var masked_DEM = ee.Image(ee.Algorithms.If(ee.Algorithms.IsEqual('NaN', SC_median),
    
    
    
    var DEM_TSLrange_outline = ee.Feature(ee.Algorithms.If(ee.Algorithms.IsEqual('NaN', SC_mean),
      ee.Feature(null, {status: 0}),
      ee.Algorithms.If(SC_stdev.gt(0).and(SC_mean.gt(0)),
        get_DEM_TSLrange_outline(snow_cover_outline).set('status', 1),
        ee.Feature(null, {status: 0})
      )
    ));
    
    
    var DEM_TSLrange_status = ee.Number(DEM_TSLrange_outline.get('status'));
    var DEM_TSLrange_area = ee.Number(ee.Algorithms.If(DEM_TSLrange_status.eq(1),
      DEM_TSLrange_outline.area(0.1).divide(1000 * 1000),
      NaN
    ));
    
    // Calculate portion of cloud cover in TSL elevation range    
    var cc_TSL_range = ee.Feature(ee.Algorithms.If(DEM_TSLrange_status.eq(1),
      DEM_TSLrange_outline.intersection(cloud_cover_outline, 0.1),
      ee.Feature(null, {status: 0})
    ));
        
    var cc_TSLrange_area = ee.Number(ee.Algorithms.If(DEM_TSLrange_status.eq(1),
      ee.Number(cc_TSL_range.area(0.1).divide(1000 * 1000)),
      NaN
    ));
        
    var cc_TSLrange_percent = ee.Number(ee.Algorithms.If(DEM_TSLrange_status.eq(1),
      cc_TSLrange_percent = cc_TSLrange_area.divide(DEM_TSLrange_area).multiply(100),
      100
    ));
    
    // Calculate portion of total class coverage in TSL elevation range
    var tcc_TSL_range = ee.Feature(ee.Algorithms.If(DEM_TSLrange_status.eq(1),
      DEM_TSLrange_outline.intersection(total_class_cov_outline, 0.1),
      ee.Feature(null, {status: 0})
    ));
        
    var tcc_TSLrange_area = ee.Number(ee.Algorithms.If(DEM_TSLrange_status.eq(1),
      ee.Number(tcc_TSL_range.area(0.1).divide(1000 * 1000)),
      NaN
    ));
        
    var tcc_TSLrange_percent = ee.Number(ee.Algorithms.If(DEM_TSLrange_status.eq(1),
      tcc_TSLrange_percent = tcc_TSLrange_area.divide(DEM_TSLrange_area).multiply(100),
      100
    ));
        
    // Resemble the TSL determination from DEM for display
    var create_percentile_masked_dem = function(snow_cover_outline) {
      var alos_clipped = alos_DEM.clip(snow_cover_outline);
      var lowest_percentiles = alos_clipped.reduceRegion({reducer: ee.Reducer.percentile([tsl_median_binsize]),
                                      geometry: snow_cover_outline.geometry(),
                                      bestEffort: true}).get('AVE');
        
      var masked_DEM = alos_clipped.updateMask(alos_clipped.lte(ee.Number(lowest_percentiles)));
      return masked_DEM;
    };

    var masked_DEM = ee.Image(ee.Algorithms.If(ee.Algorithms.IsEqual('NaN', SC_median),
      ee.Image().byte(),
      create_percentile_masked_dem(snow_cover_outline)
    ));
    
    
    /* // Debugging messages 
    print('Sat scene', selected_scene);
    print('Snow cover:', snow_cover);
    print('SCoutline + stats', snow_cover_outline);
    // print('TSL_range_min', TSL_range_min);
    // print('TSL_range_max', TSL_range_max);
    print('DEM_TSLrange_area', DEM_TSLrange_area);
    print('Clouds', cloud_cover);
    print('cc_TSLrange_area', cc_TSLrange_area);
    print('cc_TSLrange_percent', cc_TSLrange_percent);
    */
    
    
    // Prepare visualizations
    var debrisViz = {min: 0, max: 1, palette: ['FF6600', '855208']};
    var iceViz = {min: 0, max: 1, palette: ['990000', '084f85']};
    var cloudViz = {min: 0, max: 2, palette: ['999999', '999999']};
    var snowViz = {min: 0, max: 1, palette: ['000000', 'FFFFFF']};
    // var topoViz = {min: 4000, max: 7000, palette: ['33AAAA', 'FF6666']};
    // var DEM_TSL_Viz = {min: 0, max: 9000, opacity: 0.7, palette: ['990000', '990000']}; 
    // var topoViz = {min: 4000, max: 7000, palette: ['376a4e', 'fae394', '8a5117', '7c7772', 'ffffff']};
    var topoViz = {min: 0, max: 9000, palette: ['2fc3c8', '2fc3c8']};
    
    
    // Add layers to map 
    Map.addLayer(selected_glacier, {color: '000000'}, 'RGI outline');
    Map.addLayer(selected_scene, {bands: ['SWIR1', 'NIR', 'R']}, 'LS ' + value);
    Map.addLayer(SC_masked, snowViz, 'Snow cover raster ' + value);
    Map.addLayer(IC_masked, iceViz, 'Ice cover raster ' + value);
    Map.addLayer(CC_masked, cloudViz, 'Cloud cover raster ' + value);
    Map.addLayer(DC_masked, debrisViz, 'Debris cover raster ' + value);
    // Map.addLayer(DEM_TSLrange_outline, {color: '990000'}, 'DEM_TSL_range ' + value);
    Map.addLayer(cc_TSL_range, {color: '333333'}, 'CC TSL Range ' + value);
    Map.addLayer(masked_DEM, topoViz, 'Masked DEM ' + value);
    
    // join different layers into one
    var class_collection = ee.ImageCollection([SC_masked.select('SWIR1').remap([0, 1],[0, 1]).rename('Classification').toInt(),
                                               IC_masked.select('SWIR1').remap([0, 1], [0, 2]).rename('Classification').toInt(),
                                               DC_masked.select('SWIR1').remap([0, 1], [0, 3]).rename('Classification').toInt(),
                                               CC_masked.select('covered').remap([0, 1], [0, 4]).rename('Classification').toInt()])
    var mosaic = class_collection.mosaic()
                                 .set('RGIId', rgi_id)
                                 .set('GLACIER_AREA', glacier_area);
    print(mosaic)
    Map.addLayer(mosaic, {min: 1, max: 4, palette: ['FFFFFF', '084f85', '855208', '999999']})
    
    var ls_version = scene_id.getInfo().split('_')[0]
    var projection = selected_scene.select('SWIR1').projection().getInfo();
    print(projection)
    
    var filename = ee.String('MANTRA_SurfClass_').cat(ls_version).cat('_').cat(rgi_id.replace('.', '-')).cat('_').cat(scene_date);
    print(filename)
    
    // Export the image, specifying the CRS, transform, and region.
    Export.image.toDrive({
      image: mosaic.clip(selected_glacier),
      description: filename.getInfo(),
      region: selected_glacier,
      crs: projection.crs,
      fileFormat: 'GeoTIFF',
      scale: 30 //this needs to be changed according to the Landsat resolution which should be 30 for most cases
    });
    
    /*  // Debugging: Uncomment to show all generated vector outlines
    Map.addLayer(total_class_cov_outline, {color: '123456'}, 'Coverage outline ' + value);
    Map.addLayer(snow_cover_outline, {color: 'AAAAAA'}, 'Snow outline ' + value);
    Map.addLayer(ice_cover_outline, {color: '0000FF'}, 'Ice outline ' + value);
    Map.addLayer(debris_cover_outline, {color: '988523'}, 'Debris outline ' + value);
    Map.addLayer(cloud_cover_outline, {color: '666666'}, 'Cloud outline ' + value);
    Map.addLayer(debris_ice_cover_outline, {color: '987234'}, 'Debris+Ice outline ' + value);
    */
    
    // (Re-)Populate the info-panels
    panel.remove(2);  panel.remove(3); panel.remove(4); panel.remove(5); panel.remove(6); panel.remove(7); panel.remove(8); panel.remove(9); panel.remove(10); panel.remove(11); panel.remove(12); panel.remove(13); panel.remove(14); panel.remove(15);
    
    var output_glacier_area = Math.round(glacier_area.getInfo() * 100) / 100;
    var output_mean = Math.round(snow_cover_outline.getInfo().properties.STATS_MEAN * 100) / 100;
    var output_stddev = Math.round(snow_cover_outline.getInfo().properties.STATS_STDEV * 100) / 100;
    var output_tsl_class_cov = Math.round(tcc_TSLrange_percent.getInfo() * 100) / 100;
    var output_cc_port = Math.round(CC_total_port.getInfo() * 100) / 100;
    var output_tsl_cc_port = Math.round(cc_TSLrange_percent.getInfo() * 100) / 100;
    
    var output_class_coverage = class_coverage.getInfo();
    if (output_class_coverage > 100) {
      output_class_coverage = "100.00";
    } else {
      output_class_coverage = Math.round(output_class_coverage * 100) / 100;
    }
    
    panel.widgets().set(1, ui.Label('RGI ID: ' + selected_glacier.get('RGIId').getInfo()));
    panel.widgets().set(2, ui.Label('LS ID: ' + scene_id.getInfo()));
    panel.widgets().set(3, ui.Label('Glacier area: ' + output_glacier_area + ' km²'));
    panel.widgets().set(4, ui.Label('Glacier min. elev.: ' + glacier_DEM_min.getInfo() + ' m a.s.l.'));
    panel.widgets().set(5, ui.Label('Glacier max. elev.: ' + glacier_DEM_max.getInfo() + ' m a.s.l.'));
    panel.widgets().set(6, ui.Label('TSLA median: ' + snow_cover_outline.getInfo().properties.STATS_MEDIAN  + ' m a.s.l.'));
    panel.widgets().set(7, ui.Label('TSLA mean: ' + output_mean  + ' m a.s.l.'));
    panel.widgets().set(8, ui.Label('TSLA std. dev.: ' + output_stddev  + ' m'));
    panel.widgets().set(9, ui.Label('TSLA min: ' + snow_cover_outline.getInfo().properties.STATS_MIN + ' m a.s.l.'));
    panel.widgets().set(10, ui.Label('TSLA max: ' + snow_cover_outline.getInfo().properties.STATS_MAX  + ' m a.s.l.')); 
    panel.widgets().set(11, ui.Label('Cloud-c. portion total area: ' + output_cc_port + '%'));   
    panel.widgets().set(12, ui.Label('Cloud-c. in TSL range: ' + output_tsl_cc_port+ '%'));
    panel.widgets().set(13, ui.Label('Class. coverage: ' + output_class_coverage+ '%'));
    panel.widgets().set(14, ui.Label('Class. c. in TSL range: ' + output_tsl_class_cov+ '%')); 
    panel.widgets().set(15, ui.Label('Ice/debris max. elev.: ' + DIC_max.getInfo() + ' m a.s.l.')); 
    // panel.widgets().set(15, ui.Label('Ice+Debris Max (masl): ' + debris_ice_cover_outline.getInfo()['properties']['STATS_MAX'])); 
    
    /* // Debugging: show all area measurements
    panel.remove(15); panel.remove(16); panel.remove(17); panel.remove(18);
    panel.widgets().set(15, ui.Label('Ice-covered (km²): ' + IC_area.getInfo()));
    panel.widgets().set(16, ui.Label('Snow-covered (km²): ' + SC_area.getInfo()));
    panel.widgets().set(17, ui.Label('Debris-covered (km²): ' + DC_area.getInfo()));
    panel.widgets().set(18, ui.Label('Cloud-covered (km²): ' + CC_area.getInfo()));
    */
    
    
  }
});


// Set a place holder.
imageSelect.setPlaceholder('Choose Landsat scene to display ...');

panel.add(imageSelect);
ui.root.add(panel); 

// set position of panel
var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px'
  }
});
 
// Create legend title
var legendTitle = ui.Label({
  value: 'Glacier surface classification',
  style: {
    fontWeight: 'bold',
    fontSize: '14px',
    margin: '0 0 4px 0',
    padding: '0'
    }
});
 
// Add the title to the panel
legend.add(legendTitle);
 
// Creates and styles 1 row of the legend.
var makeRow = function(fillColor, outlineColor, name) {
 
      // Create the label that is actually the colored box.
      var colorBox = ui.Label({
        style: {
          backgroundColor: '#' + fillColor,
          //color: '#' + outlineColor,
          border: '1px solid #' + outlineColor,
          // Use padding to give the box height and width.
          padding: '8px',
          margin: '0 0 4px 0'
        }
      });

 
      // Create the label filled with the description text.
      var description = ui.Label({
        value: name,
        style: {margin: '0 0 4px 6px'}
      });
 
      // return the panel
      return ui.Panel({
        widgets: [colorBox, description],
        layout: ui.Panel.Layout.Flow('horizontal')
      });
};
 
//  Palette with the colors
var fillPalette     = ['2fc3c8', 'FFFFFF', '084f85', '855208', '999999', '666666'];
var outlinePalette  = ['CCCCCC', 'CCCCCC', 'CCCCCC', 'CCCCCC', 'CCCCCC', '000000'];
 
// name of the legend
var names = ee.List(['TSL pixels', 'Snow','Ice','Debris, rock', 'Clouds', 'Clouds in TSL elev. range']);

// Add color and and names
for (var i = 0; i < names.length().getInfo(); i++) {
  legend.add(makeRow(fillPalette[i], outlinePalette[i], names.get(i).getInfo()));
}  
 
// add legend to map (alternatively you can also print the legend to the console)
Map.add(legend);




/***************************************
// Modify the following to enable selecting glaciers by clicking:
// -> After click action, filter the RGI dataset to the clicked coords


// Load and display an NDVI image.
var ndvi = ee.ImageCollection('LANDSAT/LC8_L1T_8DAY_NDVI')
    .filterDate('2014-01-01', '2015-01-01');
var vis = {min: 0, max: 1, palette: ['99c199', '006400']};
Map.addLayer(ndvi.median(), vis, 'NDVI');

// Configure the map.
Map.setCenter(-94.84497, 39.01918, 8);
Map.style().set('cursor', 'crosshair');

// Create a panel and add it to the map.
var inspector = ui.Panel([ui.Label('Click to get mean NDVI')]);
Map.add(inspector);

Map.onClick(function(coords) {
  // Show the loading label.
  inspector.widgets().set(0, ui.Label({
    value: 'Loading...',
    style: {color: 'gray'}
  }));

  // Determine the mean NDVI, a long-running server operation.
  var point = ee.Geometry.Point(coords.lon, coords.lat);
  var meanNdvi = ndvi.reduce('mean');
  var sample = meanNdvi.sample(point, 30);
  var computedValue = sample.first().get('NDVI_mean');

  // Request the value from the server.
  computedValue.evaluate(function(result) {
    // When the server returns the value, show it.
    inspector.widgets().set(0, ui.Label({
      value: 'Mean NDVI: ' + result.toFixed(2),
    }));
  });
});

*/




