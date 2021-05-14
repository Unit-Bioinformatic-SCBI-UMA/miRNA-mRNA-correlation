#!/usr/bin/env python
import requests
import json
import re 
import sys
import os
import csv
from xml.etree import ElementTree
from collections import defaultdict
import pandas as pd
import mygene
import itertools

def read_genes(file):
    genes=''
    lista = []
    genName=''
    dicc={}
    if os.path.exists(file):
        if os.stat(file).st_size != 0:
            with open(file) as f:
            	for line in f:
                    print(line)
                    gen,name=line.split()
                    if not name=="gene":
                        if genName!='':
                            genName=genName + ',' + name
                        else:
                            genName=genName + name
            dicc={'genes':genName}
            return(dicc)
        else:
            raise Exception("El fichero se encuentra vacio")
    else:
        raise Exception("El fichero no existe")	
def iteration_genes_drugs(parameters,file):
	output = open(file, 'w')
	url='http://62.217.127.8/DianaTools/tarbaseApi?'
	r = requests.get(url, params=parameters)
	root = ElementTree.fromstring(r.content)
	for child in root.iter('*'):
		items = child.attrib
		if items.get('geneName') != None and items.get('mirnaName')!=None:
			output.write(str(items.get('geneName'))+ "\t"+str(items.get('mirnaName'))+ "\n")
def interaction_DIANA_microT_CDS(parameters, file):
    output = open(file, 'w')
    url='http://62.217.127.8/DianaTools/microT_CDSApi?'
    for genesID in parameters:
        r = requests.get(url, params='genes='+str(genesID))
        if r.status_code == 200:
            root = ElementTree.fromstring(r.content)
        for child in root.iter('*'):
            items = child.attrib
            if items.get('geneName') != None and items.get('mirnaName')!=None:
                genName=items.get('geneName').split(' ')[1].replace("("," ").replace(")"," ")
                output.write(str(genName)+ "\t"+str(items.get('mirnaName'))+ "\n")
def orderFile(file,fileOutput):
	lines = open(file, 'r').readlines()
	output = open(fileOutput, 'w')
	for line in sorted(lines, key=lambda line: line.split()[0]):
		output.write(str(line))
def buildDiccionary(file):
    l= defaultdict(list)
    if os.path.exists(file):
        if os.stat(file).st_size != 0:
            with open(file) as f:
                for line in f:
                    gen,name=line.split()
                    if not name=="gene":
                        l[gen].append(name)
            return(l)
        else:
            raise Exception("El fichero se encuentra vacio")
    else:
        raise Exception("El fichero no existe")
def filterDiccionary(dicc1,dicc2):
    l= defaultdict(list)
    listafinal=[]
    for key, value in dicc1.items():
        if key in dicc2.keys():            
            for i in value:
                if i in dicc2.get(key):
                    l[key].append(i)
    return(l)
def convertToDataframe(dicc, filename):
    file = open(filename, 'w')
    for key, values in dicc.items():
        for value in values:
            file.write(str(key)+"\t"+str(value)+"\n")
def getEmsemblIDs(genes, file):
    idsEmsembl=open(file, 'w')
    #idsEmsembl={}
    mg = mygene.MyGeneInfo()
    for gene in genes:
        result = mg.query(gene, scopes="symbol", fields=["ensembl"], species="human", verbose=False)
        hgnc_name = gene
        for hit in result["hits"]:
            if "ensembl" in hit and "gene" in hit["ensembl"]:
                #idsEmsembl[hgnc_name]=hit["ensembl"]["gene"]
                idsEmsembl.write("%s\t%s\n" % (hgnc_name, hit["ensembl"]["gene"]))
def buildInteraction(file, output):
    lista=[]
    IDs=buildDiccionary(file)
    for x in IDs.values():
        lista=lista+x
    interaction_DIANA_microT_CDS(lista, output)
def writeCommonDataBases(list1,list2, filename):
    if type(list1) is {}.values().__class__:
        results = list(itertools.chain(*(list(list1)+list(list2))))
    else:
        results=list(list1)+list(list2)

    with open(filename, mode='wt', encoding='utf-8') as myfile:
        print(set(results))
        myfile.write('\n'.join(set(results)))
        myfile.write('\n')
def writeCommonsInteraction(dic, filename):
    f = open(filename, "w")
    for key, value in dic.items():
        for micro in value:
            f.write(key+'\t'+micro+'\n')

if __name__ == '__main__':
    ######unfiltered interactions
    parameters=read_genes(sys.argv[1])
    #prediccion
    geneList=parameters['genes'].split(',')
    getEmsemblIDs(geneList, './idsEmsembl.dat')

    buildInteraction('./idsEmsembl.dat', './predictionInteraction.dat')

    iteration_genes_drugs(parameters,'./experimentDatabaseInteraction.dat')
    orderFile('./experimentDatabaseInteraction.dat', './experimentDatabaseInteractionOrder.dat')

    ############################
    ##########commons
    interaccionesExperimento=buildDiccionary(sys.argv[2])
    interactionPredict= buildDiccionary('./predictionInteraction.dat')
    interactionExperimentDataBase=buildDiccionary('./experimentDatabaseInteractionOrder.dat')
    commonPredict=filterDiccionary(interaccionesExperimento, interactionPredict)
    commonExperiment=filterDiccionary(interaccionesExperimento, interactionExperimentDataBase)

    writeCommonDataBases(commonPredict.keys(),commonExperiment.keys(), 'commonGenesfilter.dat')
    writeCommonDataBases(commonPredict.values(), commonExperiment.values(),'commonMicrosfilter.dat')

    writeCommonsInteraction(commonPredict, "commonPredict.txt")
    writeCommonsInteraction(commonExperiment, "commonExperiment.txt")








