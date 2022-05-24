/****************************************************

Threshold definitions and functions used by the "MountAiN glacier Transient snowline Retrieval Algorithm" (MANTRA).

2022, David Loibl, Inge Grünberg

*****************************************************/

/* THRESHOLD DEFINITIONS  */

// Thresholds SNOW:                 NDSI min,  NDSI max,   Ratio1 min,   Ratio1 max,  Ratio2 min,   Ratio2 max
var th_snow_l1      = ee.List([     -99,       -0.6,       -0.9,         -0.15,       -0.15,         99]);  
var th_snow_l2      = ee.List([     -99,       -0.6,       -0.9,         -0.15,       -0.15,         99]);  
var th_snow_l3      = ee.List([     -99,       -0.6,       -0.9,         -0.15,       -0.15,         99]);  
var th_snow_l4      = ee.List([     -99,       -0.6,       -0.9,         -0.15,       -0.15,         99]);  
var th_snow_l5      = ee.List([     -99,       -0.6,       -0.9,         -0.15,       -0.15,         99]);
var th_snow_l7      = ee.List([     -99,       -0.6,       -1.5,          0,          -0.15,         99]);
var th_snow_l8      = ee.List([     -99,       -0.6,       -0.7,         -0.15,       -0.15,         99]);
var th_unknown      = ee.List([      99,        99,         99,           99,          99,           99]);

// Thresholds ICE:                  NDSI min,  NDSI max,   Ratio1 min,   Ratio1 max,  Ratio2 min,   Ratio2 max
var th_ice_l1       = ee.List([     -99,       -0.45,      -0.8,          0.3,        -0.6,         -0.15]); 
var th_ice_l2       = ee.List([     -99,       -0.45,      -0.8,          0.3,        -0.6,         -0.15]); 
var th_ice_l3       = ee.List([     -99,       -0.45,      -0.8,          0.3,        -0.6,         -0.15]); 
var th_ice_l4       = ee.List([     -99,       -0.45,      -0.8,          0.3,        -0.6,         -0.15]); 
var th_ice_l5       = ee.List([     -99,       -0.45,      -0.8,          0.3,        -0.6,         -0.15]);
var th_ice_l7       = ee.List([     -99,       -0.5,       -0.4,          0.55,       -0.7,         -0.15]);
var th_ice_l8       = ee.List([     -99,       -0.6,       -0.4,          0.55,       -0.7,         -0.15]);
var th_unknown      = ee.List([      99,        99,         99,           99,          99,           99]);  

// Thresholds DEBRIS:               NDSI min,  NDSI max,   Ratio1 min,   Ratio1 max,  Ratio2 min,   Ratio2 max
var th_debris_l1    = ee.List([     -0.2,       99,         0.1,          1,          -0.8,         -0.2]); 
var th_debris_l2    = ee.List([     -0.2,       99,         0.1,          1,          -0.8,         -0.2]); 
var th_debris_l3    = ee.List([     -0.2,       99,         0.1,          1,          -0.8,         -0.2]); 
var th_debris_l4    = ee.List([     -0.2,       99,         0.1,          1,          -0.8,         -0.2]);  
var th_debris_l5    = ee.List([     -0.2,       99,         0.1,          1,          -0.8,         -0.2]);
var th_debris_l7    = ee.List([     -0.1,       99,         0,            0.6,        -0.8,         -0.2]);
var th_debris_l8    = ee.List([     -0.2,       99,         0,            0.6,        -0.8,         -0.3]);
var th_unknown      = ee.List([      99,        99,         99,           99,          99,           99]);


/* FUNCTIONS */

// Preprocess LS bands and calculate band ratios
var get_band_ratios = function(image) {
  // Minimum temperature threshold [K]
  var T_threshold       = 263.15;

  var green           = image.select('G');
  var red             = image.select('R');
  var nir             = image.select('NIR');
  var swir            = image.select('SWIR1');
  var bt              = image.select('T');
  var T_mask          = bt.lt(T_threshold);
  var bt_no_freeze    = T_mask.multiply(T_threshold).add(bt.multiply(T_mask.not()));
  // var     = ee.Algorithms.If(bt.gt(273.15), bt, 273.15)
  var bt_corr         = bt_no_freeze.subtract(230).divide(70);
  
  var ndsi_term1      = swir.subtract(green);
  var ndsi_term2      = swir.add(green);
  var ndsi            = ndsi_term1.divide(ndsi_term2);
  
  var ratio_1_term1   = bt_corr.subtract(green.add(nir));
  var ratio_1_term2   = bt_corr.add(green.add(nir));
  var ratio_1         = ratio_1_term1.divide(ratio_1_term2);
  
  var ratio_2_term1   = green.subtract(bt_corr);
  var ratio_2_term2   = green.add(bt_corr);
  var ratio_2         = ratio_2_term1.divide(ratio_2_term2);
  
  return ee.List([ndsi, ratio_1, ratio_2]);
};

exports.get_band_ratios = get_band_ratios;


// Identify debris-coverd regions using band thresholds
exports.debris_identification = function(image) {
  var band_ratios   = get_band_ratios(image);
  var ndsi          = ee.Image(band_ratios.get(0));
  var ratio_1       = ee.Image(band_ratios.get(1));
  var ratio_2       = ee.Image(band_ratios.get(2));

  var thresholds      = ee.List(ee.Algorithms.If(ee.Algorithms.IsEqual(image.get('SENSOR'), 'LANDSAT_8'),
                            th_debris_l8,
                            ee.Algorithms.If(ee.Algorithms.IsEqual(image.get('SENSOR'), 'LANDSAT_7'),
                              th_debris_l7,
                              ee.Algorithms.If(ee.Algorithms.IsEqual(image.get('SENSOR'), 'LANDSAT_5'),
                                th_debris_l5,
                                ee.Algorithms.If(ee.Algorithms.IsEqual(image.get('SENSOR'), 'LANDSAT_4'),
                                  th_debris_l4,
                                  ee.Algorithms.If(ee.Algorithms.IsEqual(image.get('SENSOR'), 'LANDSAT_3'),
                                    th_debris_l3,
                                    ee.Algorithms.If(ee.Algorithms.IsEqual(image.get('SENSOR'), 'LANDSAT_2'),
                                      th_debris_l2,
                                      ee.Algorithms.If(ee.Algorithms.IsEqual(image.get('SENSOR'), 'LANDSAT_1'),
                                        th_debris_l1,
                                        th_unknown
                                      )
                                    )
                                  )
                                )
                              )
                            )
                        ));

  var debris          = ndsi.gt(ee.Number(thresholds.get(0)))
                            .and(ndsi.lt(ee.Number(thresholds.get(1))))
                            .and(ratio_1.gt(ee.Number(thresholds.get(2))))
                            .and(ratio_1.lt(ee.Number(thresholds.get(3))))
                            .and(ratio_2.gt(ee.Number(thresholds.get(4))))
                            .and(ratio_2.lt(ee.Number(thresholds.get(5))));
  
  var rgi_id          = image.get('RGIId');
  var scene_date      = image.get('LS_DATE');
  
  debris              = debris.set('RGIId', rgi_id)
                              .set('LS_DATE', scene_date);
  return debris;
};


// Identify clean ice regions using band thresholds
exports.ice_identification = function(image) {
  var band_ratios   = get_band_ratios(image);
  var ndsi          = ee.Image(band_ratios.get(0));
  var ratio_1       = ee.Image(band_ratios.get(1));
  var ratio_2       = ee.Image(band_ratios.get(2));
  
  var thresholds      = ee.List(ee.Algorithms.If(ee.Algorithms.IsEqual(image.get('SENSOR'), 'LANDSAT_8'),
                            th_ice_l8,
                            ee.Algorithms.If(ee.Algorithms.IsEqual(image.get('SENSOR'), 'LANDSAT_7'),
                              th_ice_l7,
                              ee.Algorithms.If(ee.Algorithms.IsEqual(image.get('SENSOR'), 'LANDSAT_5'),
                                th_ice_l5,
                                ee.Algorithms.If(ee.Algorithms.IsEqual(image.get('SENSOR'), 'LANDSAT_4'),
                                  th_ice_l4,
                                  ee.Algorithms.If(ee.Algorithms.IsEqual(image.get('SENSOR'), 'LANDSAT_3'),
                                    th_ice_l3,
                                    ee.Algorithms.If(ee.Algorithms.IsEqual(image.get('SENSOR'), 'LANDSAT_2'),
                                      th_ice_l2,
                                      ee.Algorithms.If(ee.Algorithms.IsEqual(image.get('SENSOR'), 'LANDSAT_1'),
                                        th_ice_l1,
                                        th_unknown
                                      )
                                    )
                                  )
                                )
                              )
                            )
                        ));
  
  var ice             = ndsi.gt(ee.Number(thresholds.get(0)))
                            .and(ndsi.lt(ee.Number(thresholds.get(1))))
                            .and(ratio_1.gt(ee.Number(thresholds.get(2))))
                            .and(ratio_1.lt(ee.Number(thresholds.get(3))))
                            .and(ratio_2.gt(ee.Number(thresholds.get(4))))
                            .and(ratio_2.lt(ee.Number(thresholds.get(5))));
                            
  var rgi_id          = image.get("RGIId");
  var scene_date      = image.get('LS_DATE');
  
  ice                 = ice.set("RGIId", rgi_id)
                            .set('LS_DATE', scene_date);
  return ice;
};

// Identify snow-covered regions using band thresholds
exports.snow_identification = function(image) {
  var band_ratios   = get_band_ratios(image);
  var ndsi          = ee.Image(band_ratios.get(0));
  var ratio_1       = ee.Image(band_ratios.get(1));
  var ratio_2       = ee.Image(band_ratios.get(2));
  
  var thresholds      = ee.List(ee.Algorithms.If(ee.Algorithms.IsEqual(image.get('SENSOR'), 'LANDSAT_8'),
                            th_snow_l8,
                              ee.Algorithms.If(ee.Algorithms.IsEqual(image.get('SENSOR'), 'LANDSAT_7'),
                              th_snow_l7,
                              ee.Algorithms.If(ee.Algorithms.IsEqual(image.get('SENSOR'), 'LANDSAT_5'),
                                th_snow_l5,
                                ee.Algorithms.If(ee.Algorithms.IsEqual(image.get('SENSOR'), 'LANDSAT_4'),
                                  th_snow_l4,
                                  ee.Algorithms.If(ee.Algorithms.IsEqual(image.get('SENSOR'), 'LANDSAT_3'),
                                    th_snow_l3,
                                    ee.Algorithms.If(ee.Algorithms.IsEqual(image.get('SENSOR'), 'LANDSAT_2'),
                                      th_snow_l2,
                                      ee.Algorithms.If(ee.Algorithms.IsEqual(image.get('SENSOR'), 'LANDSAT_1'),
                                        th_snow_l1,
                                        th_unknown
                                      )
                                    )
                                  )
                                )
                              )
                            )
                        ));
  
  var snow            = ndsi.gt(ee.Number(thresholds.get(0)))
                            .and(ndsi.lt(ee.Number(thresholds.get(1))))
                            .and(ratio_1.gt(ee.Number(thresholds.get(2))))
                            .and(ratio_1.lt(ee.Number(thresholds.get(3))))
                            .and(ratio_2.gt(ee.Number(thresholds.get(4))))
                            .and(ratio_2.lt(ee.Number(thresholds.get(5))));                            
  
  var rgi_id          = image.get("RGIId");
  var scene_date      = image.get('LS_DATE');
  
  snow                = snow.set("RGIId", rgi_id)
                            .set('LS_DATE', scene_date);
  return snow;
};



// Create a vector outline from boolean classified raster
exports.boolraster_to_outline = function(boolraster) {
  var area_true = boolraster.eq(1);
  var area_false = boolraster.eq(0);
  var band_names = boolraster.bandNames();
  var band_name1 = band_names.get(0);
  var band_name2 = band_names.get(-1);
  
  var boolraster_zones = ee.Image.cat([area_true, area_false]).select(
    [band_name1, band_name2], ['covered', 'free']);
  boolraster_zones = boolraster_zones.updateMask(boolraster_zones.neq(0));
  
  var vector_outline = boolraster_zones.reduceToVectors({
    reducer: ee.Reducer.max(), 
    geometryType: 'polygon',
    labelProperty: 'RGIId',
    scale: 30,
    bestEffort: true
  });

  var dict = {
    RGIId: boolraster.get('RGIId'), 
    LS_ID: boolraster.get('LS_ID'),
    LS_DATE: boolraster.get('LS_DATE'),
    GLACIER_AREA: boolraster.get('GLACIER_AREA')
    
  };
  
  vector_outline = ee.Feature(vector_outline.geometry(), dict);
  
  return vector_outline;
};


// Get minimum and maximum elevation from DEM data
exports.get_DEM_min_max = function(DEM) {
  var reducer_min_max = ee.Reducer.min();
  reducer_min_max = reducer_min_max.combine({reducer2: ee.Reducer.max(), sharedInputs: true});

  var DEM_min_max_stats = DEM.reduceRegion({reducer: reducer_min_max,
                                  geometry: DEM.geometry(),
                                  bestEffort: true});
  return DEM_min_max_stats;                                  
};

// Calculate illumination based on DEM-based topography and Sun position
exports.calculate_illumination = function(image){
  var extent = image.geometry(); //.bounds();
  extent = extent.buffer(1000);
  
  //Select DEM, calculate aspect and slope          
  var DEM = ee.Image("JAXA/ALOS/AW3D30_V1_1").select('AVE');
  var DEMs = DEM.clip(extent);
  var SLP_deg = ee.Terrain.slope(DEMs);
  var SLP = to_radians(SLP_deg);
  var ASP_deg = ee.Terrain.aspect(DEMs);
  var ASP = to_radians(ASP_deg);
  var cos_SLP = SLP.cos();
  
  var AZ;     
  var ZE;
  
  AZ = ee.Number(image.get('SUN_AZIMUTH'));
  ZE = ee.Number(ee.Number(90).subtract(image.get('SUN_ELEVATION')));

  var AZ_R = to_radians(ee.Image(ee.Number(AZ))).clip(extent);
  var ZE_R = to_radians(ee.Image(ee.Number(ZE))).clip(extent);
  
  var cos_ZE = ZE_R.cos();
  var cos_ZE_SLP = (ZE_R.cos()).multiply((SLP).cos());
  var cos_VA = ee.Image(0).clip(extent).cos();
  
  // Calculate illumination as cos(Z)⋅cos(s) + sin(Z)⋅sin(s)⋅cos(a - a')
  var IL1 = ((AZ_R.subtract(ASP)).cos()).multiply(SLP.sin()).multiply(ZE_R.sin())
           .add((ZE_R.cos()).multiply(SLP.cos()));
  
  var il_threshold = 0.35;
  var IL2 = IL1.where(IL1.lte(il_threshold), 0); // set negative illumination to 0
  var IL3 = IL2.where(IL2.gt(il_threshold), 1);
  var IL = IL1.mask(IL3).select([0], ['IL']); 
  var masked_image = image.mask(IL3); 

  masked_image = masked_image.updateMask(masked_image.select('T').gt(0)); 
  
  return masked_image;
};


// Convert image from degree to radians
var to_radians = function(image) {
  return image.toFloat().multiply(Math.PI).divide(180);
}
