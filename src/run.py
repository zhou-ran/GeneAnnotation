#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/21 2:11 PM
__author__ = 'Zhou Ran'

from multiprocessing import Pool
import click
from mRNA import mRNA

import logging

from logging import handlers

logger = logging.getLogger()
handler = logging.StreamHandler()
fh = handlers.RotatingFileHandler(
    'plot.log',
    mode='a+',
    encoding="utf-8",
    maxBytes=5 * 1024 * 1024,
    backupCount=2,
    delay=0
)
fh.setLevel(logging.DEBUG)
formatter = logging.Formatter(
    '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.addHandler(fh)
logger.setLevel(logging.DEBUG)


def process(arg):
    geneid, file, painfo = arg
    chr, pa, strand = painfo.split(':')
    try:
        minfo = mRNA.gene(
            geneid,
            file,
            10
        )

        interinfo = minfo.siteintersect(minfo.txlst, int(pa), strand)
        for isoform, info in interinfo.items():
            for type, interval in info.items():
                tmpinterval = []
                for subinter in interval:
                    subinter = list(subinter)
                    subinter[0] = list(map(lambda x: ','.join(map(str, x)), zip(*subinter[0])))
                    tmpinterval.append(tuple(subinter))
                info[type] = tmpinterval
        interinfo['geneid'] = geneid
        interinfo['pa'] = painfo
    except:
        return ""

    return dict(interinfo)


def wrapperprocess(arg):
    return process(arg)


@click.command()
@click.option('--gtf',
              help='gtf file')
# @click.option('--geneid',
#               help='the ensembl geneid')
@click.option('--cores',
              type=int,
              default=1,
              help='multiple process the request')
@click.option('--pafile',
              help='The pasite filename.')
@click.option('--output',
              help='The output filename.')
def main(gtf, cores, pafile, output):
    """

    :return:
    """
    pool = Pool(processes=cores)
    res = []
    failed = open('failed.txt', 'w')
    with open(pafile) as painfo:
        for line in painfo.readlines():
            pa, geneid = line.strip().split('\t')
            arg = [geneid, gtf, pa]
            try:
                res.append(pool.apply_async(wrapperprocess, (arg,)))
            except:
                failed.write('\t'.join(line) + '\n')

    pool.close()
    pool.join()

    failed.close()
    writeorder = ['include',
                  'impair',
                  'exlude']

    outfile = open(output, 'w')
    outfile.write('\t'.join(['pA',
                             'geneid',
                             'transid',
                             'include_s',
                             'include_e',
                             'include_feature',
                             'impair_s',
                             'impair_e',
                             'impair_feature',
                             'exclude_s',
                             'exclude_e',
                             'exclude_feature']))

    for pa in res:
        pa = pa.get()
        if not pa:
            continue
        for id, info in pa.items():
            if id == 'pa': continue
            if id == 'geneid': continue
            tmp = []
            for type in writeorder:
                try:
                    tmpinfo = info[type]
                    if tmpinfo:
                        nameinfo = '|'.join(list(map(lambda x: x[1], tmpinfo)))
                        coordinfo = list(map(lambda x: '|'.join(x), zip(*map(lambda x: x[0], tmpinfo))))
                        tmp += coordinfo + [nameinfo]
                    else:
                        tmp += ['NA', 'NA']
                except KeyError:
                    tmp += ['NA', 'NA']
            outfile.write('\t'.join([pa['pa'], pa['geneid'], id, '\t'.join(tmp)]) + '\n')
    outfile.close()


if __name__ == '__main__':
    main()
