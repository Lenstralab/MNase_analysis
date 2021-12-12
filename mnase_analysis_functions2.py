import pysam
import numpy as np
import os
from tqdm.auto import tqdm,trange
import json
import subprocess
import pandas as pd
import copy
import scipy
from scipy import ndimage, signal
import matplotlib.pyplot as plt
import seaborn as sb
import math
##########
def calcNDRpropsPerGene(geneCovPerBamWT,geneInfo,outPath):
    NDRpropsPerGene = []

    for gene in tqdm(geneInfo):
        systname = gene['systname']
        name = gene['name']

        #get sum of coverage - not normalized for reads
        sumCovNotNorm = np.zeros(4000)
        for bam in geneCovPerBamWT:
            geneInfoDictBam = bam['geneInfoDictBam']

            nbrGenes = len(geneInfoDictBam)
            for i in range(0,nbrGenes):
                if geneInfoDictBam[i]['systname'] != systname:
                    continue
                try: 
                    covNotNorm = geneInfoDictBam[i]['covNotNorm']
                    sumCovNotNorm = sumCovNotNorm + covNotNorm
                except KeyError: 
                    continue

        try: plusOnePeak, minusOnePeak, smoothCov, peaks = findPlusMinusOne(sumCovNotNorm,40)
        except TooFewPeaks: 
            continue

        geneNDRprops = {}
        geneNDRprops['name'] = name
        geneNDRprops['systname'] = systname
        geneNDRprops['plusOneNucPos'] = plusOnePeak-2000
        geneNDRprops['widthNDR'] = plusOnePeak - minusOnePeak

        NDRpropsPerGene.append(geneNDRprops)

        text_file = open(outPath, 'w')
        text_file.write(json.dumps(NDRpropsPerGene))
        text_file.close()
##########
def filterLength(inBam,outBam,minLength,maxLength):
    for read in inBam:
        length = abs(read.template_length)
        if length < minLength or length > maxLength: continue
        outBam.write(read)
##########
def filterLengthDir(bamDir,outDir,minLength,maxLength):

    fileNames = os.listdir(bamDir)
    bamFiles = []
    for bamFile in fileNames:
        if bamFile.split('.')[-1]=='bam' and bamFile[0]!='.': bamFiles.append(bamFile)
            
    for bamFile in tqdm(bamFiles):
        inBam = pysam.AlignmentFile(bamDir+bamFile)
        outBamFile = bamFile.split('.bam')[0]+'_filtered_95-225bp.bam'
        outBam = pysam.AlignmentFile(outDir+outBamFile,'wb',template = inBam)
        filterLength(inBam,outBam,minLength,maxLength)
        inBam.close()
        outBam.close()
##########
def indexBamDir(bamDir):
    fileNames = os.listdir(bamDir)
    bamFiles = []
    for fileName in fileNames:
        if fileName.split('.')[-1] == 'bam': bamFiles.append(fileName)

    for fileName in tqdm(bamFiles):
        inpath = bamDir+fileName
        baiFile = inpath.split('.bam')[0]+'.bai'
        subprocess.call(['samtools','index',inpath,baiFile])
##########
def getBamList(bamDir):
    filelist = os.listdir(bamDir)
    bamFileList = []
    for fileName in filelist:
        if fileName.split('.')[-1] == 'bam': bamFileList.append(fileName)

    return bamFileList
##########
def calcCoverageNormCovChr(inPath,outPath):
    idxfilePath = inPath.split('.bam')[0]+'_idxstats.txt'

    #extract chromosome names in bam files
    text_file = open(idxfilePath)
    lines = text_file.read().split('\n')
    chromcontigs = [];
    end = len(lines[:-2])
    for i in range(0,end):
        chrom = lines[i].split('\t')[0]
        chromcontigs.append(chrom)
    nbrChroms = len(chromcontigs)

    ChromLengths = []
    bamfile = pysam.AlignmentFile(inPath)

    for chromNr in range(0,len(chromcontigs)):
        ChromLengths.append(bamfile.get_reference_length(chromcontigs[chromNr]))

    chromCovDicts = []
    for chromNr in range(0,17):
        chromCovDict = {}
        
        length = ChromLengths[chromNr]
        chromCov = np.zeros(length)
        fetchedChrom = bamfile.fetch(chromcontigs[chromNr])
        nbrReads = 0

        for read in fetchedChrom:
            if read.is_read2: continue
            nbrReads = nbrReads + 1
            readSide1 = read.reference_start
            readSide2 = read.reference_end
            readLength = read.template_length
                
            if readLength > 0:
                readStart = readSide1
                readEnd = readSide1 + readLength
            else:
                readStart = readSide2 + readLength
                readEnd = readSide2
                
            if readStart < 0: readStart = 0
            if readEnd > length: readEnd = length
            for base in range(readStart,readEnd):
                chromCov[base] = chromCov[base]+1

        normVal = np.sum(chromCov)
        covNorm = chromCov/normVal

        chromCovDict['filename'] = inPath
        chromCovDict['chromNr'] = chromNr
        chromCovDict['chromcontig'] = chromcontigs[chromNr]
        chromCovDict['nbrReads'] = nbrReads
        chromCovDict['covNotNorm'] = list(chromCov)
        chromCovDict['coverage'] = list(covNorm)

        chromCovDicts.append(chromCovDict)
        
    bamfile.close()
    text_file = open(outPath, 'w')
    text_file.write(json.dumps(chromCovDicts))
    text_file.close()
##########
def extractCoverageListNormCovChr(bamList,bamDir,outDir):
    for filename in tqdm(bamList):
        inPath = bamDir+filename
        outPath = outDir+(inPath.split('/')[-1]).split('.bam')[0]+'_chromosome_coverages.txt'
        calcCoverageNormCovChr(inPath,outPath)
##########
def idxstats(inPath,outPath):
    outfile = open(outPath,'w')
    outfile.writelines(subprocess.check_output(['samtools','idxstats',inPath]))
    outfile.close()
##########
def idxstatsList(bamList,bamDir):
    for filename in tqdm(bamList):
        outname = filename.split('.bam')[0]+'_idxstats.txt'
        idxstats(bamDir+filename,bamDir+outname)
##########
def buildGeneInfo(geneFilePath,chromcontigs,ChromLengths):
    geneinfo = pd.read_excel(geneFilePath,header=None)
    genecount = geneinfo.shape[0]

    #building list based on gene annotations
    geneInfoDicts = []
    for i in range(0,genecount):
        if geneinfo.iloc[i,2] != 'Verified': continue
        chromname = geneinfo.iloc[i,6]
        try: chromnr = int(chromname.split(' ')[1])-1
        except: continue
        contig = chromcontigs[chromnr]
        chromlen = ChromLengths[chromnr]
        TSS = geneinfo.iloc[i,9]
        TSE = geneinfo.iloc[i,10]
        name = geneinfo.iloc[i,4]
        strand = geneinfo.iloc[i,11]
        if strand == 'W': direction = 1
        elif strand == 'C': direction = -1
        systname = geneinfo.iloc[i,3]
        geneDict = {"chromosome": chromname,"chromnr":chromnr,"contig":contig,"TSS":TSS,"TSE":TSE,"name":name,
                    "chromlen":chromlen,"strand":strand,'systname':systname}
        geneInfoDicts.append(geneDict)
        
    return geneInfoDicts
##########
def addTATAInfo(TATAFilePath,geneInfo):
    TATAinfo = pd.read_excel(TATAFilePath,header=0)

    #Make dictionary of TATA and TATA-mismatch locations per gene
    TATADicts = []
    TATAmisDicts = []
    lenTATAinfo = len(TATAinfo)
    for i in range(0,lenTATAinfo):
        if TATAinfo.iloc[i,17] == 'TATA-containing':
            TATADict = {'systname':TATAinfo.iloc[i,0],'TATAstart':TATAinfo.iloc[i,6]}
            TATADicts.append(TATADict)
        elif TATAinfo.iloc[i,17] == 'TATA-less':
            TATAmisDict = {'systname':TATAinfo.iloc[i,0],'TATAmisStart':TATAinfo.iloc[i,6]}
            TATAmisDicts.append(TATAmisDict)


    #Add tatabox and tatamismatch locations to geneInfo
    for element in geneInfo:
        systname = element['systname']
        for TATADict in TATADicts:
            if TATADict['systname'] == systname:
                element['TATAstart'] = TATADict['TATAstart']
        for TATAmisDict in TATAmisDicts:
            if TATAmisDict['systname'] == systname:
                element['TATAmisStart'] = TATAmisDict['TATAmisStart']

    return geneInfo
# ##########
def getChromcontigs(inPath):
    idxfilePath = inPath.split('.bam')[0]+'_idxstats.txt'

    #extract chromosome names in bam files
    text_file = open(idxfilePath)
    lines = text_file.read().split('\n')
    chromcontigs = [];
    end = len(lines[:-2])
    for i in range(0,end):
        chrom = lines[i].split('\t')[0]
        chromcontigs.append(chrom)

    return chromcontigs
##########
def getChromLengths(inPath,chromcontigs):

    ChromLengths = []
    bamfile = pysam.AlignmentFile(inPath)

    for chromNr in range(0,len(chromcontigs)):
        ChromLengths.append(bamfile.get_reference_length(chromcontigs[chromNr]))
        
    return ChromLengths
##########
class TooFewPeaks(Exception):
    pass
##########
def findPlusMinusOne(genecov,smoothWindow):
    smoothCov =  scipy.ndimage.gaussian_filter(genecov,smoothWindow)
    peaks = scipy.signal.find_peaks(smoothCov)
    peakBeforeTSS = -1
    peakAfterTSS = -1
    nbrPeaks = len(peaks[0])
    if nbrPeaks < 2: raise TooFewPeaks
    if peaks[0][-1]-peaks[0][0] < 1000: raise TooFewPeaks
    minLoc = np.argmin(smoothCov[1000:2500])+1000
    if minLoc > peaks[0][-1]: minLoc = np.argmin(smoothCov[1000:peaks[0][-1]])+1000
    for i in range(1,len(peaks[0])):
        if peaks[0][i] > minLoc:
            plusOnePeak = peaks[0][i]
            minusOnePeak = peaks[0][i-1]
            break

    return plusOnePeak, minusOnePeak, smoothCov, peaks
##########
def calcCovBam(covPath, geneInfo, chromLengths):
    geneInfoBam = copy.deepcopy(geneInfo)
    
    #loading coverage file
    text_file = open(covPath, 'r')
    text = text_file.read()
    chromCovDicts = json.loads(text)
    text_file.close()

    #defining over which chroms to calculate
    chromnrs = []
    for geneInfoDict in geneInfo:
        chromnrs.append(geneInfoDict["chromnr"])
    chromnrsUnique = list(set(chromnrs))

    #loop over all chromosomes
    for chrom in chromnrsUnique:
        coverage = chromCovDicts[chrom]['coverage']
        coverageNotNorm = chromCovDicts[chrom]['covNotNorm']
        bases = chromLengths[chrom]
            
        for geneInfoDict in geneInfoBam:
            if geneInfoDict['chromnr'] != chrom: continue
            TSS = int(geneInfoDict["TSS"])
            if TSS-2000<0 or TSS+2000>bases: continue

            #determine coverage in window around TSS
            genecov = coverage[TSS-2000:TSS+2000]
            genecovNotNorm = coverageNotNorm[TSS-2000:TSS+2000]
            if geneInfoDict["strand"] == 'C': 
                genecov = np.flip(genecov)
                genecovNotNorm = np.flip(genecovNotNorm)
            geneInfoDict['coverage'] = genecov
            geneInfoDict['covNotNorm'] = genecovNotNorm

            #total coverage at TSS
            covTSS = np.sum(coverage[TSS-10:TSS+10])
            geneInfoDict['covTSS'] = covTSS

            #total coverage at TSS
            covTSS = np.sum(coverage[TSS-10:TSS+10])
            geneInfoDict['covTSS'] = covTSS

            #total coverage at TATA box
            try:
                TATAstart = int(geneInfoDict['TATAstart'])
                if TATAstart < TSS: covTATA = np.sum(coverage[TATAstart:TATAstart+8])
                else: covTATA = np.sum(coverage[TATAstart-8:TATAstart])
                geneInfoDict['covTATA'] = covTATA
            except KeyError:
                pass
            
            #total coverage at TATAmismatch
            try:
                TATAmisStart = int(geneInfoDict['TATAmisStart'])
                if TATAmisStart < TSS: covTATAmis = np.sum(coverage[TATAmisStart:TATAmisStart+8])
                else: covTATAmis = np.sum(coverage[TATAmisStart-8:TATAmisStart])
                geneInfoDict['covTATAmis'] = covTATAmis
            except KeyError:
                pass

            #find +/-1 nucleosomes and NDR width
            try: plusOnePeak, minusOnePeak, smoothCov, peaks = findPlusMinusOne(genecov,40)
            except TooFewPeaks: continue

    return geneInfoBam
##########
def importAlignedCovs(inDir):
    filelist = os.listdir(inDir)
    alignedCovs = []
    for fileName in filelist:
        text_file = open(inDir+fileName, 'r')
        text = text_file.read()
        alignedCov = json.loads(text)
        text_file.close()
        alignedCovs.append(alignedCov)
    return alignedCovs
##########
def importNDRprops(inPath):
    text_file = open(inPath, 'r')
    text = text_file.read()
    output = json.loads(text)
    text_file.close()

    return output
##########
def calcCovPlusOneControlAligned(geneCov,plusOneList,outPath):
    output = {}
    geneInfo = geneCov['geneInfoDictBam']
    sumcov = np.zeros(2000)
    nbrGenes = 0
    ###
    for element in geneInfo:
        plusOnePosWT = 999999999
        systName = element['systname']
        for gene in plusOneList:
            if gene['systname'] != systName: continue
            plusOnePosWT = gene['plusOneNucPos']
        if abs(plusOnePosWT) > 1000: continue
        
        startCov = int(1000+plusOnePosWT)
        endCov = int(3000+plusOnePosWT)
        coverage = element['coverage']
        
        try:
            sumcov = sumcov+element['coverage'][startCov:endCov]
            nbrGenes = nbrGenes + 1
        except KeyError: pass
    sumcov = sumcov/nbrGenes

    output['metadata'] = copy.deepcopy(geneCov['metadata'])
    output['metadata']['alignment'] = 'plusOneControl'
    output['covAligned'] = list(sumcov)

    text_file = open(outPath, 'w')
    text_file.write(json.dumps(output))
    text_file.close()
##########
##Per Bam, make a dict with all TSS Coverages and TATA coverages
def extractCovFeatures(geneCovPerBam):

    AllCovsPerBam = []
    for geneCov in tqdm(geneCovPerBam):
        geneInfo = geneCov['geneInfoDictBam']
        TSSCovs = []
        TATACovs = []
        TATAmisCovs = []
        # NDRwidths = []
        # NDRminCovs = []
        # minOneHeigths = []
        # plusOneHeigths = []
        for element in geneInfo:
            try: TSSCovs.append(element['covTSS'])
            except KeyError: pass

            try: TATACovs.append(element['covTATA'])
            except KeyError: pass

            try: TATAmisCovs.append(element['covTATAmis'])
            except KeyError: pass

            # try: NDRwidths.append(element['widthNDR'])
            # except KeyError: pass

            # try: NDRminCovs.append(element['minCoverageNDR'])
            # except KeyError: pass

            # try: minOneHeigths.append(element['minusOneNucHeigth'])
            # except KeyError: pass

            # try: plusOneHeigths.append(element['plusOneNucHeigth'])
            # except KeyError: pass

        output = {}
        output['metadata'] = copy.deepcopy(geneCov['metadata'])
        output['TSSCovs'] = TSSCovs
        output['TATACovs'] = TATACovs
        output['TATAmisCovs'] = TATAmisCovs
        # output['NDRwidths'] = NDRwidths
        # output['NDRminCovs'] = NDRminCovs
        # output['minOneHeigths'] = minOneHeigths
        # output['plusOneHeigths'] = plusOneHeigths
        AllCovsPerBam.append(output)
    return AllCovsPerBam
# ##########
def calcHistNeg(binWidth, data):
    minData = min(data)
    maxData = max(data)
    edge = 0
    binList = [0]
    while edge > minData:
        edge = edge - binWidth
        binList.insert(0,edge)

    edge = 0
    while edge < maxData:
        edge = edge + binWidth
        binList.append(edge)  

    histData = np.histogram(data,binList)
    histSum = np.sum(histData[0])
    histVals = [float(i)/histSum for i in histData[0]]
    binCenters = [i+0.5*binWidth for i in binList[:-1]]
    return [histVals,binCenters]
##########
def importCovFeatures(inDir):
    filelist = os.listdir(inDir)
    covFeatures = []
    for fileName in filelist:
        if len(fileName.split('.')[0]) == 0: continue
        text_file = open(inDir+fileName, 'r')
        text = text_file.read()
        element = json.loads(text)
        text_file.close()
        covFeatures.append(element)
    return covFeatures
##########
def writeHist(covFeatures,covGenRegion,binWidth,outDir):
    for element in covFeatures:
        
        metadata = copy.deepcopy(element['metadata'])
        
        histList = element[covGenRegion]
        histData = calcHist(binWidth,histList)
        
        outFileName = metadata['filename'].split('.bam')[0]+'_'+covGenRegion+'.json'
        outPath = outDir+outFileName
        output = {}
        output['metadata'] = metadata
        output['histData'] = histData
        output['metadata']['covFeature'] = covGenRegion
        
        text_file = open(outPath, 'w')
        text_file.write(json.dumps(output))
        text_file.close()
##########
def sortByNDR(NDRprops):
    nbrGenes = len(NDRprops)
    geneNames = []
    avNDRWidthList = []
    for geneNr in range(0,nbrGenes):
        avNDRWidthList.append(NDRprops[geneNr]['widthNDR']) 
        geneNames.append(NDRprops[geneNr]['systname'])

    sortList = np.flip(np.argsort(avNDRWidthList))
    sortedGeneNames = []
    for element in sortList:
        sortedGeneNames.append(geneNames[element])
    return sortedGeneNames
##########
def calcHeatmapMatrix(geneCov, alignmentMethod, plusOneList, sortedGeneNames,outPath):
    #get coverage in appropriate range and normalize
    geneInfo = geneCov['geneInfoDictBam']
    metadata = geneCov['metadata']
    covMatrix = []
    geneList = []
    nbrGenes = len(sortedGeneNames)
    
    for geneNbr in trange(0,nbrGenes):
        systName = sortedGeneNames[geneNbr]
        for element in geneInfo:
            if element['systname'] != systName:
                continue
                
            if alignmentMethod == 'TSS':
                labelPartOne = 'TSS'
                delta = 0
            # elif alignmentMethod == 'PlusOneLoc':
            #     labelPartOne = '+1 nucleosome'
            #     try:
            #         delta = element['plusOneNucPos']
                # except KeyError: continue
            elif alignmentMethod == 'TATA':
                labelPartOne = 'TATA-box'
                try:
                    TATApos = int(element['TATAstart'])
                    TSS = int(element['TSS'])
                    delta = (-1)*abs(TSS-TATApos)
                except KeyError: continue
            elif alignmentMethod == 'plusOneWT':
                labelPartOne = 'wt +1 nucleosome'
                plusOnePosWT = 999999999
                for gene in plusOneList:
                    if gene['systname'] != systName: continue
                    plusOnePosWT = gene['plusOneNucPos']
                if plusOnePosWT == 999999999: continue
                delta = int(plusOnePosWT)

            else: raise Error

            try:
                if abs(delta) > 1000: continue
                covRange = list(element['coverage'][1000+delta:3000+delta])
                covMatrix.append(covRange)
                geneList.append(systName)
            except KeyError: continue
    output = list([covMatrix,geneList,metadata])
    np.save(outPath,output)
##########
def plotHeatmap(pathMatrix, vmaxVal, xLabelTxt, titleTxt, showFig, saveFig, outPath):

    calcHeatMap = np.load(pathMatrix,allow_pickle=True)
    heatMap = sb.heatmap(calcHeatMap[0], yticklabels = False, cmap = 'Greys',vmax = vmaxVal)
    plt.xlabel(xLabelTxt)
    plt.ylabel('Genes')
    plt.xticks([0,500,1000,1500,2000],[-1000,-500,0,500,1000])
    plt.title(titleTxt)
    if showFig: plt.show()
    if saveFig:
        fig = heatMap.get_figure()
        fig.savefig(outPath, dpi=400, bbox_inches="tight")
    plt.close()
##########
def combineReps(heatMapArray, sortedGeneList, outPath):
    if len(heatMapArray) == 1:
        covMatrix = heatMapArray[0][0]
        updatedGeneList = sortedGeneList
    
    else:
        covMatrix = []
        updatedGeneList = []
        for gene in sortedGeneList:
            covs = []
            for heatMap in heatMapArray:
                try:
                    geneNr = heatMap[1].index(gene)
                    covs.append(heatMap[0][geneNr])
                except ValueError: pass
            
            avCov = []
            try:
                for j in range(0,len(covs[0])):
                    sumValue = 0
                    sumCounts = 0
                    for i in range(0,len(covs)):
                        sumValue = sumValue + covs[i][j]
                        sumCounts = sumCounts + 1
                    avCov.append(sumValue/sumCounts)
                covMatrix.append(avCov)
                updatedGeneList.append(gene)
            except IndexError:
                pass
    
    metadata = {}
    metadata['sugar'] = heatMapArray[0][2]['sugar']
    metadata['mnaseCond'] = heatMapArray[0][2]['mnaseCond']
    try:
        metadata['factorDepl'] = heatMapArray[0][2]['factorDepl']
    except KeyError: pass

    output = list([covMatrix,updatedGeneList,metadata])
    np.save(outPath,output)
##########
def resizeHeatmaps(heatMap1,heatMap2):
    geneList = []
    covMatrix = []
    metadata = heatMap1[2]
    
    if len(heatMap1[0]) > len(heatMap2[0]):
        geneListAdj = []
        covMatrixAdj = []
        for i in range(0,len(heatMap1[1])):
            if heatMap1[1][i] not in heatMap2[1]: continue
            covMatrixAdj.append(heatMap1[0][i])
            geneListAdj.append(heatMap1[1][i])
        heatMap1 = [covMatrixAdj,geneListAdj,metadata]

    if len(heatMap1[0]) < len(heatMap2[0]):
        geneListAdj = []
        covMatrixAdj = []
        for i in range(0,len(heatMap2[1])):
            if heatMap2[1][i] not in heatMap1[1]: continue
            covMatrixAdj.append(heatMap2[0][i])
            geneListAdj.append(heatMap2[1][i])
        heatMap2 = [covMatrixAdj,geneListAdj,metadata]
    
    return heatMap1,heatMap2
##########
def calcHeatMapLog2FoldChange(heatMap1,heatMap2,outPath):
    shapeToCalc = np.shape(heatMap1[0])
    changeHeatmap = np.zeros(shapeToCalc)
    table1 = heatMap1[0]
    table2 = heatMap2[0]

    for row in range(0,shapeToCalc[0]):
        for column in range(0,shapeToCalc[1]):
            if table1[row][column] == 0 or table2[row][column] == 0: changeHeatmap[row][column] = 0
            else: changeHeatmap[row][column] = math.log(table1[row][column]/table2[row][column],2)
    
    output = [list(changeHeatmap),heatMap1[1],heatMap1[2]]
    np.save(outPath,output)
##########
def plotHeatmapChange(pathMatrix, vmaxVal, vminVal, xLabelTxt, titleTxt, showFig, saveFig, outPath):

    calcHeatMap = np.load(pathMatrix,allow_pickle = True)
    
    heatMap = sb.heatmap(calcHeatMap[0], yticklabels = False, cmap = 'bwr', center = 0,
                         vmax = vmaxVal, vmin = vminVal)
    plt.xlabel(xLabelTxt)
    plt.ylabel('Genes')
    plt.xticks([0,500,1000,1500,2000],[-1000,-500,0,500,1000])
    plt.title(titleTxt)
    if showFig: plt.show()
    if saveFig:
        fig = heatMap.get_figure()
        fig.savefig(outPath, dpi=400, bbox_inches="tight")
    plt.close()
##########
def selectGenes(geneInfo,geneList):
    geneInfoSelected = []
    for element in geneInfo:
        if element['name'] not in geneList: continue
        geneInfoSelected.append(element)
    return geneInfoSelected
##########
def addGal4Info(gal4Path,geneInfoGAL):
    gal4bsinfo = pd.read_excel(gal4Path,header=1)

    for element in geneInfoGAL:
        name = element['name']
        for geneNr in range(0,len(gal4bsinfo)):
            gal4UASs = []
            if gal4bsinfo.iloc[geneNr,0] != name: continue
            for i in range(1,5):
                if gal4bsinfo.iloc[geneNr,i] != 0:
                    gal4UASs.append(gal4bsinfo.iloc[geneNr,i])
            element['gal4UASs'] = gal4UASs

    return geneInfoGAL
##########
def addTATAInfoGal(TATAGalLocs,geneInfoGAL):
    for element in geneInfoGAL:
        name = element['name']
        element['TATAGalStart'] = TATAGalLocs[name]
    return geneInfoGAL
##########
def calcCovGene(geneCov,geneName,outPath):
    output = {}
    geneInfo = geneCov['geneInfoDictBam']
    metadata = copy.deepcopy(geneCov['metadata'])
    for gene in geneInfo:
        if gene['name'] != geneName: continue
        coverage = gene['coverage'][1000:3000]
        TATAGalStart = gene['TATAGalStart']
        gal4UASs = gene['gal4UASs']
    
    output['metadata'] = metadata
    output['metadata']['geneName'] = geneName
    output['coverage'] = list(coverage)
    output['TATAGalStart'] = TATAGalStart
    output['gal4UASs'] = gal4UASs

    text_file = open(outPath, 'w')
    text_file.write(json.dumps(output))
    text_file.close()
##########
def importGeneCovs(inDir):
    filelist = os.listdir(inDir)
    geneCovs = []
    for filename in filelist:
        if filename.split('.')[-1]!='json': continue
        text_file = open(inDir+filename,'r')
        text = text_file.read()
        geneCov = json.loads(text)
        text_file.close()
        geneCovs.append(geneCov)
    return geneCovs
##########
def calcChangePlusOne(geneCov,plusOneList,outPath):
    geneInfo = geneCov['geneInfoDictBam']
    metadata = geneCov['metadata']
    outputListAll = []
    outputListTATA = []
    outputListTATAmis = []
    for gene in geneInfo:
        systname = gene['systname']
        for gene2 in plusOneList:
            if gene2['systname'] != systname: continue
            plusOnePos = gene2['plusOneNucPos']

        try: covNorm = gene['coverage']
        except KeyError: pass
        
        try: plusOnePeak, minusOnePeak, smoothCov, peaks = findPlusMinusOne(covNorm,40)
        except TooFewPeaks: 
            continue

        changePlusOne = (plusOnePeak-2000) - plusOnePos
        outputListAll.append(changePlusOne)
        try:
            gene['TATAstart']
            outputListTATA.append(changePlusOne)
        except KeyError: pass
        try:
            gene['TATAmisStart']
            outputListTATAmis.append(changePlusOne)
        except  KeyError: pass

    output = [{'outputAll': outputListAll,'outputTATA': outputListTATA, 'outputTATAmis': outputListTATAmis, 'metadata':metadata}]
    np.save(outPath,output)
##########
def cumToPlot(input):
    sortedvals = np.sort(input)
    xToPlot = []
    yToPlot = []
    nbr = len(input)

    for i in range(0,nbr):
        xToPlot.append(sortedvals[i])
        yToPlot.append((i+1)/(nbr*1.))
    
    return [xToPlot,yToPlot]
##########
def writeCumDistribution(covFeatures,covGenRegion,outDir):
    for element in covFeatures:
        metadata = element['metadata']
        outFileName = metadata['filename'].split('.bam')[0]+'_'+covGenRegion+'_cumulative.json'
        outPath = outDir + outFileName
        cumDistr = cumToPlot(element[covGenRegion])
    
        output = {}
        output['metadata'] = metadata
        output['cumDistr'] = cumDistr
        output['metadata']['covFeature'] = covGenRegion
        text_file = open(outPath, 'w')
        text_file.write(json.dumps(output))
        text_file.close()
##########

