#!/usr/bin/env python

import argparse
import re
import sys
import logging

############################ Arguments #########################################

parser = argparse.ArgumentParser(description='converts genome index file' + 
                                ' to ron file.')
parser.add_argument('-i', required=True, help='input')
parser.add_argument('-o', required=True, help='output')
parser.add_argument('-l', required=False, help='log file')

args = parser.parse_args()

####################### Logging Configuration ##################################

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(name)s:' + 
                            '%(funcName)s:%(message)s')

if args.l:
    file_handler = logging.FileHandler(args.l, 'w')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
else:
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

############################ Functions #########################################

def checkIndexFormat(myLines):

    '''
    Validates the fai index format.
    '''

    for line in myLines:
        ll = line.strip().split("\t")
        try:
            int(ll[1])
        except:
            return False

    return True

def index2ronList(filename=args.i):

    '''
    Reads fai index file and converts it to a list of info in ron format.
    '''

    try:
        with open(filename, "r") as f:
            lines = f.readlines()

        if checkIndexFormat(lines):
            logger.info(' Fai index format is valid. ' + 
                        'Converting fai index to ron...')

            chrList = []
            for line in lines:
                ll = line.strip().split("\t")
                newline = '(name: "' + ll[0] + '", start: 0, end: ' + \
                    ll[1] + ')'
                chrList.append(newline)

        else:
            raise ValueError

    except FileNotFoundError:
        logger.error(' File not found:', filename)

    except ValueError:
        logger.error(' Values of 2. column must be numeric.')

    else:
        logger.info(' Done!')
        return chrList

def list2ronFile(ronList, outputFile=args.o):

    '''
    Write list info into output file.
    '''
    logger.info(' Writing to file...')
    with open(outputFile, 'w') as fout:
        fout.write("[\n")
        for item in ronList:
            fout.write("%s,\n" % item)
        fout.write("]")
    logger.info(' Done!')

############################### Main ###########################################

myList = index2ronList()
list2ronFile(myList)
