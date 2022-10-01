/****************************************************

Tool to classify glacier surfaces and measure Transient Snowline Altitude (TSLA) calculation 
results from the "MountAiN glacier Transient snowline Retrieval Algorithm" (MANTRA).

Please note that individual results provided by this tool may be flawed, e.g. by cloud cover. 
It is advised to check the results in the evaluate tool.

FEATURES
- Classify surface materials based on Landsat band-ratios (snow, ice, debris/rock, clouds)
- Process data for relevant Landsat scenes from multiple missions (Landsat 4 to 8)
- Merge scenes taken at the same date that are covered by the glacier's outline
- Mask pixels with insufficient illumination
- Calculate metrics and statistics
- Write the results to a CSV in Google Drive
- Do all this in parallel for (very) many glaciers and Landsat scenes.


ACKNOWLEDGEMENTS:
Routine to merge multiple Landsat sensors into one collection intensively inspired by 
Matthias Baumann's (HU Berlin) scripts. Inge Grünberg (AWI Potsdam) developed the two 
ratios additional to the NDSI and helped finding adequate thresholds.


REQUIREMENTS AND DATA:
The script requires an RGI v6 shapefile of containing all glaciers that should be processed. 
The shapefile should be published as an "Asset" in GEE.

Input:  A list of RGI IDs for the glaciers to process.
Output: A CSV file that will be written to the current user's Google Drive.

Usage:  Set up the list of RGI IDs, start date and end according to your study aims.
        Click "Run". You should see the number of glaciers in the Console tab.
        Switch to the Tasks tab and launch the job.
        The processing may take up to several days.
        For comprehensive tasks, it is recommended to split the list of RGI IDs into multiple 
        packages and launch many jobs in parallel.
        

David Loibl, 2022

*****************************************************/


/* USER CONFIGURATION */

var package_number = "00001"                    // Will define the name of the output file.

var glacier_list = ee.List(['RGI60-13.53885']); // RGI IDs to process. 

var start    = ee.Date('2020-01-01');           // First date to include in the search for Landsat scenes
var finish   = ee.Date('2021-12-31');           // Last date to include in the search for Landsat scenes

// Path to RGI shapefile containing the glacier to investigate
var glaciers = ee.FeatureCollection("users/___your_GEE_account___/_YOUR-RGI-SHAPEFILE_");

// Path to MATNRA core library
var core = require('users/___your_GEE_account___/tsl-retrieval:core');

/* USER CONFIGURATION END */



// DO NOT MODIFY ANYTHING BELOW HERE UNLESS 
// YOU REALLY KNOW WHAT YOU ARE DOING



/* LOAD DATASETS */


var tsl_median_binsize  = 2;                    // Percentile of altitude range to obtain TSL from
var cf_threshold        = 20;                   // Minimum cloudfree portion [%] of glacier surface to conduct TSL analysis

var alos_DEM = ee.Image("JAXA/ALOS/AW3D30_V1_1").select('AVE');

var exclude_glaciers = ['RGI-15.03544'];
var tool_version = '0.8.1'

glaciers = glaciers.filter(ee.Filter.inList("RGIId", glacier_list));


/*
// For testing and debugging, select only one glacier.
// Make sure the following line is deactivated for production!
var glaciers = glaciers.filter(ee.Filter.eq('RGIId', 'RGI60-14.00005'));
*/

print(glaciers.size());


/* PREPARE LANDSAT DATA */

//Build a joined Landsat 4-8 collection ()
var bands         = ['B1', 'B2', 'B3', 'B4',  'B5',    'B7',    'B6',        'BQA'];
var band_names    = ['B',  'G',  'R',  'NIR', 'SWIR1', 'SWIR2', 'T',         'BQA'];
var l7bands       = ['B1', 'B2', 'B3', 'B4',  'B5',    'B7',    'B6_VCID_2', 'BQA'];
var l7band_names  = ['B',  'G',  'R',  'NIR', 'SWIR1', 'SWIR2', 'T',         'BQA'];
var l8bands       = ['B2', 'B3', 'B4', 'B5',  'B6',    'B7',    'B10',       'BQA'];
var l8band_names  = ['B',  'G',  'R',  'NIR', 'SWIR1', 'SWIR2', 'T',         'BQA'];

var ls4 = ee.ImageCollection('LANDSAT/LT04/C01/T1_TOA').map(function(image){return image.select(bands).rename(band_names)});
var ls5 = ee.ImageCollection('LANDSAT/LT05/C01/T1_TOA').map(function(image){return image.select(bands).rename(band_names)});
var ls7 = ee.ImageCollection('LANDSAT/LE07/C01/T1_TOA').map(function(image){return image.select(l7bands).rename(l7band_names)}).filterMetadata('GAIN_BAND_5', 'equals', 'L')  // Remove scenes where TOA was processed with high gain;
var ls8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA').map(function(image){return image.select(l8bands).rename(l8band_names)});

var landsatCollection = ee.ImageCollection(ls4.merge(ls5).merge(ls7).merge(ls8))
                  .filterDate(start, finish);

// Scene mosaic parameterization
var diff = finish.difference(start, 'day');
var temporalResolution = 1;  // days
var range = ee.List.sequence(0, diff.subtract(1), temporalResolution).map(function(day){return start.advance(day,'day')});


/* LOCAL FUNCTIONS */

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

// Prepare empty data container for non-results
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


/* PROCESSING ROUTINES */

var TSL_dataset = ee.FeatureCollection(glaciers.map(function(glacier) {
  
  var rgi_id = glacier.get('RGIId');
  var glacier_area = ee.Feature(glacier).area(0.1).divide(1000 * 1000);
  
  var LS_collection_glacier = landsatCollection.filterBounds(glacier.geometry());
  
  // Funtion for iteraton over the range of dates
  var day_list = function(date, newlist) {
    // Cast
    date = ee.Date(date);
    newlist = ee.List(newlist);
  
    // Filter collection between date and the next day
    var filtered = LS_collection_glacier.filterDate(date, date.advance(temporalResolution, 'day'));
    var LS_DATE = filtered.first().get('DATE_ACQUIRED'); //date.format('YYYY-MM-dd');
  
    return ee.List(ee.Algorithms.If(filtered.size(), newlist.add(LS_DATE), newlist));
  };
  
  // Check if multiple images exist for individual days. If so, create mosaics.
  var get_day_images = function(day) {
    var filtered = LS_collection_glacier.filterDate(ee.Date(day), ee.Date(day).advance(1, 'day'));
    var ls_id = filtered.first().get('LANDSAT_PRODUCT_ID');
    var sensor = filtered.first().get('SPACECRAFT_ID');
    var sun_azi = filtered.first().get('SUN_AZIMUTH');
    var sun_ele = filtered.first().get('SUN_ELEVATION');
    var image = ee.Algorithms.If(filtered.size().gt(1), 
                  ee.Image(filtered.mosaic().clip(glacier))
                          .set('LS_DATE', day)
                          .set('LS_ID', ls_id)
                          .set('SENSOR', sensor)
                          .set('RGIId', rgi_id)
                          .set('SUN_AZIMUTH', sun_azi)
                          .set('SUN_ELEVATION', sun_ele),
                  filtered.first().clip(glacier)
                          .set('LS_DATE', day)
                          .set('LS_ID', ls_id)
                          .set('SENSOR', sensor)
                          .set('RGIId', rgi_id)
                          .set('SUN_AZIMUTH', sun_azi)
                          .set('SUN_ELEVATION', sun_ele)
                  );
    return image;
  };
  
  // Create a new image collection of clipped and merged Landsat scenes
  var sceneList = ee.List(range.iterate(day_list, ee.List([])));
  LS_collection_glacier = ee.ImageCollection(sceneList.map(get_day_images));
  LS_collection_glacier = LS_collection_glacier.map(core.calculate_illumination);
  
  var glacier_results = LS_collection_glacier.map(function(ls_scene){
    var scene_date  = ee.Number(ls_scene.get('LS_DATE'));
    var scene_id    = ee.Number(ls_scene.get('LS_ID')); 
    var scene_sat   = ls_scene.get('SENSOR');
    
    // Apply band ratio to classify different surface cover types
    var snow_cover = ee.Image(core.snow_identification(ls_scene))
      .set('RGIId', rgi_id)
      .set('GLACIER_AREA', glacier_area);
    var SC_masked = snow_cover.updateMask(snow_cover.eq(1));    
    
    var ice_cover = ee.Image(core.ice_identification(ls_scene))
      .set('RGIId', rgi_id)
      .set('GLACIER_AREA', glacier_area);    
    var IC_masked = ice_cover.updateMask(ice_cover.eq(1));
    
    var debris_cover = ee.Image(core.debris_identification(ls_scene))
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
    //var total_area      = ee.Number(DIC_area.add(SC_area).add(CC_area));
    
    
    // For LS7 SLC off scenes: calculate coverage of total glacier (%)
    var class_coverage  = ee.Number(TCC_area.divide(glacier_area).multiply(100));
    //var class_coverage  = ee.Number(total_area.divide(glacier_area).multiply(100));
    
    // Calculate the portion of classfied surfaces areas that is cloud-covered (%)
    var CC_total_port   = ee.Number(CC_area.divide(TCC_area).multiply(100));
    
    
    // Obtain DEM stats for snow-covered and ice-plus-debris-covered regions
    snow_cover_outline  = ee.Feature(ee.Algorithms.If(SC_area.gt(0),
                            obtain_DEM_stats(snow_cover_outline),
                            create_nan_dict(snow_cover_outline)
                        ));
    
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
        
    var SC_stdev        = ee.Number(snow_cover_outline.get('STATS_STDEV'));
    var SC_mean         = ee.Number(snow_cover_outline.get('STATS_MEAN'));
    var SC_median       = ee.Number(snow_cover_outline.get('STATS_MEDIAN'));
    var SC_min          = ee.Number(snow_cover_outline.get('STATS_MIN'));
    var SC_max          = ee.Number(snow_cover_outline.get('STATS_MAX'));
    
    

    
    // Identify cloud-covered regions close to TSL (mean +/- 1sigma)
    var alos_glacier    = alos_DEM.clip(glacier);
    var alos_min_max    = core.get_DEM_min_max(alos_glacier);
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
    
    /*    
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
    */
    
    var result_dataset = ee.Feature(null, {
        'RGI_ID': rgi_id,
        'LS_ID: ': scene_id,
        'LS_DATE': scene_date,
        'LS_SAT': scene_sat,
        'glacier_area': glacier_area,
        'glacier_DEM_min': glacier_DEM_min,
        'glacier_DEM_max': glacier_DEM_max,
        'DIC_area': DIC_area,
        'DIC_max': DIC_max,
        'SC_stdev': SC_stdev,
        'SC_mean': SC_mean,
        'SC_median': SC_median,
        'SC_min': SC_min,
        'SC_max': SC_max,
        'SC_area': SC_area,
        'IC_area': IC_area,
        'DC_area': DC_area,
        'CC_area': CC_area,
        'CC_total_port': CC_total_port,
        'cc_TSLrange_percent': cc_TSLrange_percent,
        'class_coverage': class_coverage,
        'tool_version': tool_version
    });
    return result_dataset;
  });

  return glacier_results;

}));

// print('main.TSL_dataset', TSL_dataset.flatten())

// Flatten the FeatureCollection, so that all results are features (no sub-FeauterCollections)
var TSLs = TSL_dataset.flatten();

// Remove features that contain no TSL result
// var TSLs = TSLs.filter(ee.Filter.gt('TSL_MEDIAN', 0))

// Write the results to Google Drive
Export.table.toDrive({
    collection: TSLs,
    folder: "TSL-results", 
    //description: 'TSL-results-'+ lat_min * 10 +'-'+lat_max * 10+'-'+ lon_min * 10 +'-'+lon_max * 10,
    description: 'TSL-results-'+ package_number,
    fileFormat: 'CSV'
  });


