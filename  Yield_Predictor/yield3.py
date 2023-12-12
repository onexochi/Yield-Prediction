import geemap, ee
ee.Authenticate()
ee.Initialize()
aoi = ee.FeatureCollection("users/onexpeters/Kenya_County_Boundaries").\
filter(ee.Filter.eq('COUNTY', 'Kericho')).geometry()
def getEVI(image):
    # Compute the EVI using an expression.
    EVI = image.expression(
        '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
            'NIR': image.select('B8').divide(10000),
            'RED': image.select('B4').divide(10000),
            'BLUE': image.select('B2').divide(10000)
        }).rename("EVI")

    image = image.addBands(EVI)

    return(image)
def maskQuality(image):
    # Select the QA band.
    QA = image.select('pixel_qa')
    # Get the internal_cloud_algorithm_flag bit.
    sombra = getQABits(QA,2,2,'cloud_shadow')
    nubes = getQABits(QA,4,4,'cloud')
    #  var cloud_confidence = getQABits(QA,6,7,  'cloud_confidence')
    cirrus_detected = getQABits(QA,8,8,'cirrus_detected')
    #var cirrus_detected2 = getQABits(QA,8,8,  'cirrus_detected2')
    #Return an image masking out cloudy areas.
    return image.updateMask(sombra.eq(0)).updateMask(nubes.eq(0).updateMask(cirrus_detected.eq(0)))
def addDate(image):
    img_date = ee.Date(image.date())
    img_date = ee.Number.parse(img_date.format('YYYYMMdd'))
    return image.addBands(ee.Image(img_date).rename('date').toInt())
def getLAI(image):
    LAI = image.expression(
        '(3.618*EVI - 0.118)', {
            'EVI': image.select('EVI')
        }).rename("LAI")
    image = image.addBands(LAI)

    return(image)
def getQABits(image, start, end, mascara):
    # Compute the bits we need to extract.
    pattern = 0
    for i in range(start,end+1):
        pattern += 2**i
    # Return a single band image of the extracted QA bits, giving the     band a new name.
    return image.select([0], [mascara]).bitwiseAnd(pattern).rightShift(start)
#A function to mask out cloudy pixels.
def main():
    Sentinel_data = ee.ImageCollection('COPERNICUS/S2').filterDate("2022-03-01","2022-03-31").filterBounds(aoi).filterMetadata('CLOUDY_PIXEL_PERCENTAGE', 'less_than',10).map(getEVI)
    Lai_image = Sentinel_data.map(getLAI).map(addDate).median().clip(aoi)
    palett = ['FFFFFF', 'CE7E45', 'DF923D', 'F1B555', 'FCD163', '99B718','74A901', '66A000', '529400', '3E8601', '207401', '056201','004C00', '023B01', '012E01', '011D01', '011301']
    pall = {"min":0.5, "max":3.5, 'palette':palett}
    map1 = geemap.Map()
    map1.centerObject(aoi, 8)
    map1.addLayer(Lai_image.select('LAI'), pall, "LAI")
    map1.addLayerControl()
    print(map1)
   
if __name__ == "__main__":
	main()