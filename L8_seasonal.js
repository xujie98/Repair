//Method: Jean-François Pekel1, Nature Letter 2016 and NIR+Otsu
//Author: XU Jie, ITP CAS
// Load Landsat 8 top-of-atmosphere reflectance image.
// 3-5,6-8,9-11,12-2 was seperated to create seasonal permanent water surface

var Hankou  = ee.FeatureCollection("users/usrname/Hankou");
var L8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA')
        .filterBounds(Hankou)
        .filterDate('2013-04-01','2019-12-22')
        .select(["B1","B2","B3","B4","B5","B7","BQA"],['blue','green','red','nir','SWIR1','SWIR2','qa'])
        .map(function(image){
          var date = image.date().millis();
          return image.clip(Hankou).set('system:time_start', date);
        });

var qualityfunction = function(image){
  var nir = image.select('nir').gt(0)
  var blue = image.select('blue').gt(0)
  var red = image.select('red').gt(0)
  var green = image.select('green').gt(0)
  var quality = nir.multiply(blue).multiply(red).multiply(green)
  var qualitymin = quality.reduceRegion({
            reducer: ee.Reducer.min(),
            geometry: Hankou,
            scale: 30,
            });
  return image.set(qualitymin)
}
var L8 = L8.map(qualityfunction).filter(ee.Filter.greaterThan('nir',0));

var L8_Pekel = L8.filter(ee.Filter.calendarRange(5,9,'month'))

// parameter function:
var maskfunction_Pekel = function(image){
  //------------1. mask out water/cloud pixel -----------------------
  var NDVI = image.normalizedDifference(['red','nir']).rename('nd');
  var MNDWI = image.normalizedDifference(['green','nir']).rename('NDWI');
  var hsv = image.select(['SWIR2', 'nir', 'red']).rgbToHsv();
  image = image.addBands(NDVI).addBands(MNDWI).addBands(hsv)
  var mask_water_ori = image.expression(
    "(b('nd') > 0&&b('value') < 0.62&&(((b('hue')<((-9.867784585617413*b('nd'))+238.26034242940045))"+
      "&&(b('hue')>((-12960.000000000335*b('nd'))-12714.048607819708))||(b('hue')>((23.627546071775214*b('nd'))+255.53176874753507)))"+
      "||((b('hue')<((-54.685799109352004*b('nd'))+215.15052322834936))&&(b('hue')<((23.627546071775214*b('nd'))+255.53176874753507))"+
      "&&(b('hue')>((-7.321079389910027*b('nd'))+224.6166270396205)))||((b('hue')<((-172.0408163265306*b('nd'))+191.69646750224035))"+
      "&&(b('hue')<((-7.321079389910027*b('nd'))+224.6166270396205))&&(b('hue')>((-38.11764705882351*b('nd'))+193.8533786110101)))"+
      "||((b('hue')>((-52.06378986866776*b('nd'))+179.92232432949075))&&(b('hue')<((-879.6226415094455*b('nd'))+180.3004476242325))"+
      "&&(b('hue')<((-38.11764705882351*b('nd'))+193.8533786110101)))))? 1" + 
      //"(b('hue')>((23.627546071775214*b('nd'))+255.53176874753507))?1"+ this condition makes the image full balck.
      // So change && to || in front of it.

        ":(b('nd') < 0&&b('value')<0.62&&(((b('hue')<((-119.15098406819945*b('nd'))+180.0533162435398))"+
        "&&(b('hue')>((-994.2857142867327*b('nd'))+180.04805813312743))&&(b('hue')>((-116.5000234173271*b('nd'))+179.9633248496054)))"+
        "||((b('hue')<((-2368.4258422651174*b('nd'))+256.40879883589054))&&(b('hue')<((-116.5000234173271*b('nd'))+179.9633248496054))"+
        "&&(b('hue')>((-267.6720052547653*b('nd'))+179.97791758964533)))"+
        "||((b('hue')<((-108.07947019867622*b('nd'))+179.67747476669464))&&(b('hue')>((-2368.4258422651174*b('nd'))+256.40879883589054))"+
        "&&(b('hue')>((58.99660016815455*b('nd'))+168.09286521078695)))||((b('hue')<((-104.45621862799788*b('nd'))+179.4262481567021))"+
        "&&(b('hue')<((58.99660016815455*b('nd'))+168.09286521078695))&&(b('hue')>((-52.1565190088889*b('nd'))+172.13690440390852)))"+
        "||((b('hue')<((-52.1565190088889*b('nd'))+172.13690440390852))&&(b('hue')>((-204.2258047185466*b('nd'))+177.66958001421082))"+
        "&&(b('hue')>((37.74894387447151*b('nd'))+159.60620482085795)))))? 1" + //(b('nd') < 0&&
          ":(b('SWIR1')<0||b('SWIR2')<0)? 1"+ 
          ": 0")   
  var mask_water = ee.Image(0).clip(Hankou).where(mask_water_ori.eq(1),1).rename('mask_water') 
   
  //--------------------mask out cloud ---------------------------------
  //get the Landsat 8 Pixel Quality Assessment(pixel_qa) Bit Index
  var image_qa = image.select('qa');
 
  // Create a mask for the dual QA bit "Cloud Confidence".
  // Bits 5-6: Cloud Confidence  
  // 0: Not Determined / Condition does not exist. 
  // 1: Low, (0-33 percent confidence)
  // 2: Medium, (34-66 percent confidence)
  // 3: High, (67-100 percent confidence)
  var RADIX = 2;  // Radix for binary (base 2) data.
  
  var extractQABits = function (qaBand, bitStart, bitEnd) {
    var numBits = bitEnd - bitStart + 1;
    var qaBits = qaBand.rightShift(bitStart).mod(Math.pow(RADIX, numBits));
    return qaBits;
  };
  
  var bitStartCloudConfidence = 5;
  var bitEndCloudConfidence = 6;
  var qaBitsCloudConfidence = extractQABits(image_qa, bitStartCloudConfidence, bitEndCloudConfidence);
  // Test for clouds, based on the Cloud Confidence value.
  var mask_cloud_ori = qaBitsCloudConfidence.gte(2)     
  var mask_cloud     = ee.Image(0).select('constant').clip(Hankou).where(mask_cloud_ori.eq(1),1).rename('mask_cloud')
  //--------------------------------------------------------------------

//mask those pixels from the image
image = image.addBands(mask_water).addBands(mask_cloud);

  return image;
};

var L8_mask_Pekel = L8_Pekel.map(maskfunction_Pekel);
print('L8_mask_Pekel',L8_mask_Pekel);

//=========================================================================================
var L8_NIR = L8.filter(ee.Filter.calendarRange(10,12,'month')).merge(L8.filter(ee.Filter.calendarRange(1,4,'month')))
                .sort('system:time_start')  // Sort chronologically in ascending order.

var histofunction = function(image){
    var histogram = image.select('nir').reduceRegion({
      reducer: ee.Reducer.histogram(255, 2)
        .combine('mean', null, true)
        .combine('variance', null, true), 
      geometry: Hankou, 
      scale: 30,
      bestEffort: true
    });

  image = image.set(histogram)
  return image;
}
var L8_histo = L8_NIR.map(histofunction);
print('L8_histo',L8_histo)

var maskfunction_NIR = function(image){
  // Return the DN that maximizes interclass variance in B5 (in the region).
  var otsu = function(histogram) {
    var counts = ee.Array(ee.Dictionary(histogram).get('histogram'));
    var means = ee.Array(ee.Dictionary(histogram).get('bucketMeans'));
    var size = means.length().get([0]);
    var total = counts.reduce(ee.Reducer.sum(), [0]).get([0]);
    var sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]);
    var mean = sum.divide(total);
  
    var indices = ee.List.sequence(1, size);
  
    // Compute between sum of squares, where each mean partitions the data.
    var bss = indices.map(function(i) {
      var aCounts = counts.slice(0, 0, i);
      var aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0]);
      var aMeans = means.slice(0, 0, i);
      var aMean = aMeans.multiply(aCounts)
        .reduce(ee.Reducer.sum(), [0]).get([0])
        .divide(aCount);
      var bCount = total.subtract(aCount);
      var bMean = sum.subtract(aCount.multiply(aMean)).divide(bCount);
      return aCount.multiply(aMean.subtract(mean).pow(2)).add(
           bCount.multiply(bMean.subtract(mean).pow(2)));
    });
  
  // Return the mean value corresponding to the maximum BSS.
  return means.sort(bss).get([-1]);
};

var threshold = otsu(image.get('nir_histogram'));

var mask_water = image.select('nir').lt(threshold).rename('mask_water')

  //--------------------mask out cloud ---------------------------------
  //get the Landsat 8 Pixel Quality Assessment(pixel_qa) Bit Index
  var image_qa = image.select('qa');
 
  // Create a mask for the dual QA bit "Cloud Confidence".
  // Bits 5-6: Cloud Confidence  
  // 0: Not Determined / Condition does not exist. 
  // 1: Low, (0-33 percent confidence)
  // 2: Medium, (34-66 percent confidence)
  // 3: High, (67-100 percent confidence)
  var RADIX = 2;  // Radix for binary (base 2) data.
  
  var extractQABits = function (qaBand, bitStart, bitEnd) {
    var numBits = bitEnd - bitStart + 1;
    var qaBits = qaBand.rightShift(bitStart).mod(Math.pow(RADIX, numBits));
    //Map.addLayer(qaBits, {min:0, max:(Math.pow(RADIX, numBits)-1)}, 'qaBits');
    return qaBits;
  };
  
  var bitStartCloudConfidence = 5;
  var bitEndCloudConfidence = 6;
  var qaBitsCloudConfidence = extractQABits(image_qa, bitStartCloudConfidence, bitEndCloudConfidence);
  // Test for clouds, based on the Cloud Confidence value.
  var mask_cloud = qaBitsCloudConfidence.gte(2).rename('mask_cloud')
  //--------------------------------------------------------------------

  // //mask those pixels from the image
  image = image.addBands(mask_water).addBands(mask_cloud);

  return image;
};

var L8_mask_NIR = L8_histo.filter(ee.Filter.gt('nir_mean', 0)).map(maskfunction_NIR);

// merge collections
var L8_mask = L8_mask_Pekel.merge(L8_mask_NIR).sort("system:time_start")
print('L8_mask',L8_mask)

// // Export the water mask image
// // Get the number of images.
// var L8_mask_id = L8_mask.size().getInfo()-1;
// var list = L8_mask.toList(L8_mask_id);
//   for (var i = 0;i<L8_mask_id;i++){
//     var image = ee.Image(list.get(i));
//     var date = image.date().format('yyyy-MM-dd').getInfo()
//     var label= '9-Hknir_'+i.toString()+'_L8_'+date 
//       // Export NDVI        
//       Export.image.toDrive({ 
//             image: image.select('nir','red','green'),
//             description: label,
//             fileNamePrefix: label,
//             scale: 30,
//             region:Hankou,
//             //crs : 'EPSG:32649'  //  WGS 84 / UTM zone 49N
//             });
// }

// // Export the water mask image
// // Get the number of images.
// var L8_mask_id = L8_mask.size().getInfo()-1;
// var list = L8_mask.toList(L8_mask_id);
//   for (var i = 0;i<L8_mask_id;i++){
//     var image = ee.Image(list.get(i));
//     var date = image.date().format('yyyy-MM-dd').getInfo()
//     var label= '9-Hkcolor_'+i.toString()+'_L8_'+date 
//       // Export NDVI        
//       Export.image.toDrive({ 
//             image: image.select('red','green','blue'),
//             description: label,
//             fileNamePrefix: label,
//             scale: 30,
//             region:Hankou,
//             //crs : 'EPSG:32649'  //  WGS 84 / UTM zone 49N
//             });
// }
  
// // Export the water mask image
// // Get the number of images.
// var L8_mask_id = L8_mask.size().getInfo()-1;
// var list = L8_mask.toList(L8_mask_id);
//   for (var i = 0;i<L8_mask_id;i++){
//     var image = ee.Image(list.get(i));
//     var date = image.date().format('yyyy-MM-dd').getInfo()
//     var label= '9-Hkori_'+i.toString()+'_L8_'+date 
//       // Export NDVI        
//       Export.image.toDrive({ 
//             image: image.select('mask_water'),
//             description: label,
//             fileNamePrefix: label,
//             scale: 30,
//             region:Hankou,
//             crs : 'EPSG:32649'  //  WGS 84 / UTM zone 49N
//             });
// }

// // Export the water mask image
// // Get the number of images.
// var L8_mask_id = L8_mask.size().getInfo()-1;
// var list = L8_mask.toList(L8_mask_id);
//   for (var i = 0;i<L8_mask_id;i++){
//     var image = ee.Image(list.get(i));
//     var date = image.date().format('yyyy-MM-dd').getInfo()
//     var label= '9-Hkcloud_'+i.toString()+'_L8_'+date 
//       // Export NDVI        
//       Export.image.toDrive({ 
//             image: image.select('mask_cloud'),
//             description: label,
//             fileNamePrefix: label,
//             scale: 30,
//             region:Hankou,
//             crs : 'EPSG:32649'  //  WGS 84 / UTM zone 49N
//             });


//------------  2. reduce the image collection to one image by summing all the rasters  ----
//--------------------  and select out the minum count in water area ------------------
var Hankou_water  = ee.FeatureCollection("users/christinaluintel/Hankou100751water");
var spring_pixelcounts = ee.Image(L8_mask.filter(ee.Filter.calendarRange(3,5,'month')).select('mask_water').sum().rename('counts')).toInt()
var spring_min =  spring_pixelcounts                       
                        .reduceRegion({
                        reducer: ee.Reducer.min(),
                        geometry: Hankou_water,
                        scale: 30,
                        });
print(spring_min)                       
var image =  spring_pixelcounts
var label= 'Hk_counts_L8_spring'
Export.image.toDrive({ 
  image: image.toFloat(),
  description: label,
  fileNamePrefix: label,
  scale: 30,
  crs : 'EPSG:32649',  //  WGS 84 / UTM zone 49N // original projection
  region: Hankou,
})

var count = ee.Number(spring_min.get('counts'))
var image_permanent = spring_pixelcounts.gt(count)
var label= 'Hk_L8_'+count.getInfo()+'fill_spring'
Export.image.toDrive({ 
  image: image_permanent,
  description: label,
  fileNamePrefix: label,
  scale: 30,
  region:Hankou,
  crs : 'EPSG:32649'  //  WGS 84 / UTM zone 49N // original projection
}); 

Map.centerObject(Hankou,15);
Map.addLayer(image,{"min":0,"max":count.getInfo(),"palette":["white","blue"]},'spring_count');
Map.addLayer(image_permanent,{"min":0,"max":1,"palette":["white","grey"]},'spring_permanent')

var summer_pixelcounts = ee.Image(L8_mask.filter(ee.Filter.calendarRange(6,8,'month')).select('mask_water').sum().rename('counts')).toInt()
var summer_min =  summer_pixelcounts                        
                        .reduceRegion({
                            reducer: ee.Reducer.min(),
                            geometry: Hankou_water,
                            scale: 30,
                        }); 

var image =  summer_pixelcounts
var label= 'Hk_counts_L8_summer'
Export.image.toDrive({ 
  image: image.toFloat(),
  description: label,
  fileNamePrefix: label,
  scale: 30,
  crs : 'EPSG:32649',  //  WGS 84 / UTM zone 49N // original projection
  region: Hankou,
})

var count = ee.Number(summer_min.get('counts'))
var image_permanent = summer_pixelcounts.gt(count)
var label= 'Hk_L8_'+count.getInfo()+'fill_summer'
Export.image.toDrive({ 
  image: image_permanent,
  description: label,
  fileNamePrefix: label,
  scale: 30,
  region:Hankou,
  crs : 'EPSG:32649'  //  WGS 84 / UTM zone 49N // original projection
});

var vector = image_permanent.addBands(summer_pixelcounts).reduceToVectors({
      geometry: Hankou,
      crs : 'EPSG:32649',
      //geometryInNativeProjection: ture,
      scale: 30,
      geometryType: 'polygon',
      eightConnected: false,
      reducer: ee.Reducer.sum()
      })
var vectorlabel= '9-Hkpm_vector_summer_'+count.getInfo()+'_L8'
Export.table.toDrive({ 
        collection: vector,
        description: vectorlabel,
        fileNamePrefix: vectorlabel,
        fileFormat:"SHP",
        selectors: "label",
        //region:Hankou,
});
vector = null

var autumn_pixelcounts = ee.Image(L8_mask.filter(ee.Filter.calendarRange(9,11,'month')).select('mask_water').sum().rename('counts')).toInt()
var autumn_min =  autumn_pixelcounts                        
                        .reduceRegion({
                            reducer: ee.Reducer.min(),
                            geometry: Hankou_water,
                            scale: 30,
                        });
                        
var image =  autumn_pixelcounts
var label= 'Hk_counts_L8_autumn'
Export.image.toDrive({ 
  image: image.toFloat(),
  description: label,
  fileNamePrefix: label,
  scale: 30,
  crs : 'EPSG:32649',  //  WGS 84 / UTM zone 49N // original projection
  region: Hankou,
})

var count = ee.Number(autumn_min.get('counts'))
var image_permanent = autumn_pixelcounts.gt(count)
var label= 'Hk_L8_'+count.getInfo()+'fill_autumn'
Export.image.toDrive({ 
  image: image_permanent,
  description: label,
  fileNamePrefix: label,
  scale: 30,
  region:Hankou,
  crs : 'EPSG:32649'  //  WGS 84 / UTM zone 49N // original projection
}); 
                    
var winter_pixelcounts = ee.Image(L8_mask.filter(ee.Filter.calendarRange(1,2,'month')).merge(L8_mask.filter(ee.Filter.calendarRange(12,12,'month'))).select('mask_water').sum().rename('counts')).toInt()
var winter_min =  winter_pixelcounts                        
                        .reduceRegion({
                            reducer: ee.Reducer.min(),
                            geometry: Hankou_water,
                            scale: 30,
                        });
                        
var image =  winter_pixelcounts
var label= 'Hk_counts_L8_winter'
Export.image.toDrive({ 
  image: image.toFloat(),
  description: label,
  fileNamePrefix: label,
  scale: 30,
  crs : 'EPSG:32649',  //  WGS 84 / UTM zone 49N // original projection
  region: Hankou,
})

var count = ee.Number(winter_min.get('counts'))
var image_permanent = winter_pixelcounts.gt(count)
var label= 'Hk_L8_'+count.getInfo()+'fill_winter'
Export.image.toDrive({ 
  image: image_permanent,
  description: label,
  fileNamePrefix: label,
  scale: 30,
  region:Hankou,
  crs : 'EPSG:32649'  //  WGS 84 / UTM zone 49N // original projection
});                        

var L8_mask_spring = L8_mask.filter(ee.Filter.calendarRange(3,5,'month')).map(function(image){
  var spring_threshold = ee.Image(ee.Number(spring_min.get('counts'))).select('constant').rename('threshold').clip(Hankou);
  return image.addBands(spring_pixelcounts).addBands(spring_threshold)
})

var L8_mask_summer = L8_mask.filter(ee.Filter.calendarRange(6,8,'month')).map(function(image){
  var summer_threshold = ee.Image(ee.Number(summer_min.get('counts'))).select('constant').rename('threshold').clip(Hankou);
  return image.addBands(summer_pixelcounts).addBands(summer_threshold)
})

var L8_mask_autumn = L8_mask.filter(ee.Filter.calendarRange(9,11,'month')).map(function(image){
  var autumn_threshold = ee.Image(ee.Number(autumn_min.get('counts'))).select('constant').rename('threshold').clip(Hankou);
  return image.addBands(autumn_pixelcounts).addBands(autumn_threshold)
})

var L8_mask_winter = L8_mask.filter(ee.Filter.calendarRange(1,2,'month')).merge(L8_mask.filter(ee.Filter.calendarRange(12,12,'month'))).map(function(image){
  var winter_threshold = ee.Image(ee.Number(winter_min.get('counts'))).select('constant').rename('threshold').clip(Hankou);
  return image.addBands(winter_pixelcounts).addBands(winter_threshold)
})

var L8_mask = L8_mask_spring.merge(L8_mask_summer).merge(L8_mask_autumn).merge(L8_mask_winter).sort("system:time_start")
print('L8_mask',L8_mask)

//--------------- 3. fill the water pixels under cloud -----------------------------------
var fillfunction = function(image){
  var mask_water_whole = image.expression("(b('mask_water')<1&&b('mask_cloud')>0&&b('counts')>b('threshold'))?1"+
                                        ":b('mask_water')").rename('mask_water_whole')
  // residual non-water pixels 
  var mask = mask_water_whole.eq(0)
  var masked = mask.updateMask(mask).rename('masked')
  var patchsize = masked.connectedPixelCount({
    maxSize: 25, eightConnected: false
  }).rename('patchsize');  
  mask_water_whole = mask_water_whole.addBands(patchsize)
  var fill_mask = mask_water_whole.expression("((b('patchsize')<25&&b('mask_water_whole')<1)?1"+ 
                                            ": b('mask_water_whole'))").rename('fill_mask')  //block parts gap filling
                                            
  //-------------------- get the 900m2 area band---------------------------------
  var area = ee.Image.pixelArea().reproject({crs:'EPSG:2163',scale:30}); 
  var water01 = fill_mask.eq(1);
  var fillArea = water01.multiply(area).rename('fillArea').reproject({crs:'EPSG:32649',scale:30});  //block parts gap filling
  var stats = fillArea.reduceRegion({
    reducer: ee.Reducer.sum(), 
    geometry: Hankou, 
    scale: 30,
  });
  return image.addBands(fill_mask).addBands(fillArea).set(stats); //
}
var L8_fill = L8_mask.map(fillfunction).filter(ee.Filter.greaterThan('fillArea',900))
print('L8_fill',L8_fill)

  Map.centerObject(Hankou,15);
  Map.addLayer(L8_fill.first().select('fill_mask'),{opacity:0.5},'fill_mask'); //image_water
  
// // Export the water mask image
// // Get the number of images.
// var L8_mask_id = L8_mask.size().getInfo()-1;
// var list = L8_mask.toList(L8_mask_id);
//   for (var i = 0;i<L8_mask_id;i++){
//     var image = ee.Image(list.get(i));
//     var date = image.date().format('yyyy-MM-dd').getInfo()
//     var label= '9_Hkfill_nir-otsu_'+i.toString()+'_ls8_'+date 
//     // Export NDVI        
//     Export.image.toDrive({ 
//         image: image.select('fill_mask'),
//         description: label,
//         fileNamePrefix: label,
//         scale: 30,
//         region:Hankou,
//         crs : 'EPSG:32649'  //  WGS 84 / UTM zone 49N // original projection
//         });
//   }

//-------------3.2 add the maxium area to the image collection------------
//Convert the zones of water to vectors.
var maxareafunction = function(image){
  var vectors = image.reduceToVectors({
      geometry: Hankou,
      scale: 30,
      geometryType: 'polygon',
      eightConnected: false,
      labelProperty: 'zone',
      reducer: ee.Reducer.sum()
  });
  var maxwater = vectors.aggregate_array('sum')
  var lista = ee.List(maxwater)
  var list  = lista.sort()
  var counts = vectors.size().subtract(1)
  var highest =  list.get(counts);
  return image.set('fillmaxArea', highest);
}
    
var L8_maxarea = L8_fill.select('fill_mask','fillArea').map(maxareafunction);
print('L8_maxarea',L8_maxarea)  
// // Export the water vector
// var L8_maxarea_id = L8_maxarea.size().getInfo()-1;
// var list = L8_maxarea.toList(L8_maxarea_id);
//   for (var i = 0;i<L8_maxarea_id;i++){
//     var image = ee.Image(list.get(i));
//     var vector = image.reduceToVectors({
//         geometry: Hankou,
//         crs : 'EPSG:32649',
//         //geometryInNativeProjection: ture,
//         scale: 30,
//         geometryType: 'polygon',
//         eightConnected: false,
//         reducer: ee.Reducer.sum()
//     })

//     var date = image.date().format('yyyy-MM-dd').getInfo()
//     var label= '9_Hkfill_'+i.toString()+'_ls7_'+date 
//     // Export NDVI        
//       Export.image.toDrive({ 
//           image: image.select('fill_mask'),
//           description: label,
//           fileNamePrefix: label,
//           scale: 30,
//           region:Hankou,
//           crs : 'EPSG:32649'
//           });
      
//     var vectorlabel= '9-vector_'+i.toString()+'_L8_'+date 
//       Export.table.toDrive({ 
//           collection: vector,
//           description: vectorlabel,
//           fileNamePrefix: vectorlabel,
//           fileFormat:"SHP",
//           selectors: "label",
//           //region:Hankou,
//           });
//   }

Map.centerObject(Hankou,15);
Map.addLayer(L8_fill.first().select('fillArea'),{opacity:0.5},'fillArea'); //image_water
Map.addLayer(L8_fill.first().select('mask_cloud'), {opacity:0.5},'mask_cloud'); //image_water
Map.addLayer(L8_fill.first().select('mask_water'), {opacity:0.5},'mask_water'); //image_water

  
var vector = L8_maxarea.first().reduceToVectors({
      geometry: Hankou,
      //crs: 'EPSG:32649',//crs: image.select("fill_mask").projection(),
      scale: 30,  
      geometryType: 'polygon',
      eightConnected: false,
      bestEffort: true,
      //labelProperty: 'zone',
      reducer: ee.Reducer.sum()
    });
var display = ee.Image(0).updateMask(0).paint(vector, '000000', 3);
Map.addLayer(display, {palette: '000000'}, 'vectors');