from PyQt5.QtWidgets import QApplication, QFileDialog
from PyQt5 import uic
from osgeo import gdal
import numpy as np
import cv2
import os, sys, re


def getRasterData(inras):
    gdalData = gdal.Open(inras, gdal.GA_ReadOnly)
    wkt = gdalData.GetProjection()
    wktList = wkt.split('EPSG')
    lastEPSG = wktList[len(wktList) - 1]
    epsg = "EPSG:" + re.findall(r'\d+', lastEPSG)[0]
    rasinfo = gdalData.GetGeoTransform()
    print ("reading chanel 1 ...")
    #gdalBand = gdalData.GetRasterBand(1)
    xsize = gdalData.RasterXSize
    ysize = gdalData.RasterYSize
    print ("width: " + str(xsize) + ", height: " + str(ysize))
    xul = rasinfo[0]
    yul = rasinfo[3]
    px_x = rasinfo[1]
    px_y = rasinfo[5]
    print ("x: " + str(px_x) + ", y: " + str(px_y))
    xdr = xul + xsize*px_x
    ydr = yul + ysize*px_y
    ohvat = str(xul) + ', ' + str(ydr) + ', ' + str(xdr) + ', ' + str(yul)
    print ("Ohvat: " + ohvat)
    return(gdalData.ReadAsArray(), epsg, wkt, xsize, ysize, rasinfo, ohvat)


def outputRes(Arr, file_name):
    dst_ds = driver.Create(file_name, col, row, nBands, gdal.GDT_Byte)
    dst_ds.SetGeoTransform(rasinfo)
    dst_ds.SetProjection(wkt)
    for i in range(nBands):
        out = dst_ds.GetRasterBand(i+1)
        out.WriteArray(Arr[i])
    out.FlushCache()
    dst_ds = None

    outDir = os.path.split(file_name)[0] + "/"
    basename = os.path.basename(file_name)
    band_names = []
    for i in range(nBands):
        band_names.append(outDir + basename[:basename.rfind('.')] + "_b" + str(i+1) + ".tif")

    gdalData = gdal.Open(file_name, gdal.GA_ReadOnly)
    for i in range(nBands):
        gdalBand = gdalData.GetRasterBand(i+1)
        band = gdalBand.ReadAsArray()
        dst_ds = driver.Create(band_names[i], col, row, 1, gdal.GDT_Byte)
        dst_ds.SetGeoTransform(rasinfo)
        dst_ds.SetProjection(wkt)
        out = dst_ds.GetRasterBand(1)
        out.WriteArray(band)
        out.FlushCache()
        dst_ds = None

    return


def isDigit(val):
    try:
        x = float(val)
        return(True)
    except ValueError:
        return(False)


def onBoard_algo():
    global col, row, rasinfo, wkt

    homeDir = QFileDialog.getExistingDirectory (None, "Выберите каталог с зашумлёнными изображениями:", Init_path, options = QFileDialog.ShowDirsOnly)
    dirUpLev = homeDir[0:homeDir.rfind("/")] + "/"
    homeDir = homeDir + "/"
    print("Home dir: " + homeDir)
    if homeDir == "/":return

    outputDir = QFileDialog.getExistingDirectory (None, "Куда сохраняем результат:", dirUpLev, options = QFileDialog.ShowDirsOnly)
    outputDir = outputDir + "/"
    print ("Output dir: " + outputDir)
    if outputDir == "/":return

    noiseDirs = os.listdir(homeDir)
    imageFiles = dict()

    for noiseDir in noiseDirs:
        imfiles = os.listdir(homeDir + noiseDir)
        if len(imfiles) == 0: noiseDirs.pop(noiseDirs.index(noiseDir))

    noIm = len(noiseDirs)        # noIm - kol-wo izobr., 1 noise = 1 izobr.

    for noiseDir in noiseDirs:
        imfiles = os.listdir(homeDir + noiseDir)
        #if len(imfiles) == 0: continue
        imageFiles[noiseDir] = imfiles

    imfile1 = homeDir + noiseDir + "/" + imageFiles[noiseDir][0]
    _, epsg, wkt, col, row, rasinfo, ohvat = getRasterData(imfile1)

    Ia_fname = outputDir + "I_average.tif"
    Iet_fname = outputDir + "I_etalon-average.tif"

    Ia = np.array([],dtype=np.uint8)
    for i in range(nBands):
        arr1Bands = np.array([])
        for noiseDir in noiseDirs:
            img_fname = homeDir + noiseDir + "/" + imageFiles[noiseDir][i]
            #print("img_fname: " + img_fname)
            img = cv2.imread(img_fname)
            img = img[:,:,0]
            arr1Bands = np.append(arr1Bands, img)
        arr1Bands = arr1Bands.reshape(noIm, row, col)
        srar = np.mean(arr1Bands, axis=0).astype(np.uint8)
        Ia = np.append(Ia, srar)
    Ia = Ia.reshape(nBands, row, col)
    outputRes(Ia, Ia_fname)
    Ea = np.sum(Ia**2)
    #Ia_tr = np.transpose(Ia, (1,2,0))

    Iglob = list()          # array of images
    for key in imageFiles:
        #print("image: " + key)
        I = np.array([],dtype=np.uint8)
        for elem in imageFiles[key]:
            if 'RGB' in elem: continue
            img_fname = homeDir + key + "/" + elem
            #print("img_fname: " + img_fname)
            img = cv2.imread(img_fname)
            img = img[:,:,0]
            I = np.append(I, img)
        I = I.reshape(nBands, row, col)
        Iglob.append(I)

    E = np.array([])
    for i in range(noIm):
        Im = Iglob[i]
        Energy = np.sum(Im**2)
        E = np.append(E, Energy)

    r = np.array([])
    for i in range(noIm):
        #I_tr = np.transpose(Iglob[i], (1,2,0))
        #kkor = (1/np.sqrt(Ea*E[i])) * np.sum(I_tr * Ia_tr)
        kkor = (1/np.sqrt(Ea*E[i])) * np.sum(Iglob[i] * Ia)
        r = np.append(r, kkor)
    print("kkor: " + str(r))

    porogIsDigit = False
    while not porogIsDigit:
        KKporog_form.exec_()
        r_porog_txt = KKporog_form.lineEdit.text()
        if isDigit(r_porog_txt):
            #print("digit")
            porogIsDigit = True
        else:
            #print("no digit")
            porogIsDigit = False
            KKporog_form.lineEdit.setStyleSheet('color: #820404')

    r_porog = float(r_porog_txt)
    print("kkor porog = " + str(r_porog))

    ixs = np.where(r >= r_porog)
    ixs = np.array(ixs[0])
    M = ixs.shape[0]                    # num of images with gt porog
    print("отобрано " + str(M) + " изображений")

    if M > 0:
        Isum = Iglob[ixs[0]].astype(np.uint16)
        if M > 1:
            for i in range(1, M):
                Isum += Iglob[ixs[i]]
        Iet = Isum / M
        Iet = Iet.astype(np.uint8)
        outputRes(Iet, Iet_fname)

    return


app = QApplication(sys.argv)

#scriptPath = "D:/work/FCO/pycode/"
scriptPath = os.path.dirname(os.path.abspath(__file__))
scriptPath = scriptPath.replace('\\', "/") + "/"
print("script path: " + scriptPath)
form_path = scriptPath + "Kkor_porog_form.ui"
KKporog_form = uic.loadUi(form_path)

driver = gdal.GetDriverByName("GTiff")
Init_path = "D:/"
nBands = 4
#r_porog = 1.05

onBoard_algo()
