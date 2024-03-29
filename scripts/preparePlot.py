#!/usr/bin/env python

import sys
import os
import argparse
import logging

def getWhInfo(wh_log_file):
    info = {'score' : '-', 'downs' : '-', 'time' : '-'}
    with open(wh_log_file, 'r') as whlog:
        for line in whlog.readlines():
            if line.startswith("MEC cost:"):
                info["score"] = int(line.strip().split(':')[1])
            if line.startswith("Time spent selecting reads:"):
                info["downs"] = float(line.strip().replace(" ", "").split(':')[1][:-1])
            if line.startswith("Time spent phasing:"):
                info["time"] = float(line.strip().replace(" ", "").split(':')[1][:-1])
        return info

def getScore(log_file):
    with open(log_file, 'r') as log:
        for line in log.readlines():
            line = line.strip()
            if line.startswith("* ") and " | " in line:
                cleanline = line.split(" | ")[1]
                if cleanline.startswith("OPTIMUM:"):
                    return int(cleanline.split(":")[1])

def getMEC(mec_file):
    with open(mec_file, 'r') as mf:
        line = mf.readline().rstrip()
        return line.split(":")[1].replace(" ", "") if line else '-'

def getQUAL(diff_file):
    qual = {'ser' : '-', 'ham' : '-'}
    with open(diff_file, 'r') as d:
        nl = 0
        for l in d.readlines():
            nl += 1
            if nl == 23:
                qual["ser"]=(l.rstrip().replace(" ", "").split(":")[1]).replace("%", "")
            if nl == 27:
                qual["ham"]=(l.rstrip().replace(" ", "").split(":")[1]).replace("%", "")
    return qual

def getTimeMem(log_file):
    info = {'time' : '-', 'mem' : '-'}
    with open(log_file, 'r') as lf:
        for line in lf.readlines():
            line = line.strip()
            if line.startswith("Elapsed (wall clock) time"):
                wallv = line.split(": ")[1].split(":")
                hours = 0 if len(wallv) == 2 else int(wallv[-3])
                mins = int(wallv[-2])
                secs = float(wallv[-1])
                info["time"] = ((hours*60) + mins)*60+secs
            if line.startswith("Maximum resident set size"):
                info["mem"] = int(line.split(":")[1])/(1024)#*1024)
        return info

def getFeasibility(logfile) :
    with open(logfile,'r') as lines :
        for line in lines :
            if line.startswith(' --- no solution at alpha = 1e-2') :
                return 'no'
    return 'yes'

def main():
    coverages = ['cov{}'.format(c) for c in range(15, 65, 5)]
    mergings = ['merged_e15_m25_t{}_n3'.format(t) for t in '6 17'.split()]

    parser = argparse.ArgumentParser(prog = "preparePlot",
                                     description = "Prepare Haplotyping data.",
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--inck', help = "Increase k HapChat dir.",
                        required = True, dest = 'inck_dir')
    parser.add_argument('-p', '--pre-input', help = "Input WIF file dir.",
                        required = True, dest = 'wif_dir')
    parser.add_argument('-w', '--whatshap', help = "WhatsHap dir.",
                        required = False, dest = 'wh_dir')
    parser.add_argument('-c', '--hapcol', help = 'HapCol dir.',
                        required = False, dest = 'hc_dir')
    parser.add_argument('-u', '--hapcut2', help = 'HapCUT2 dir.',
                        required = False, dest = 'hapcut2_dir')
    parser.add_argument('-o', '--out-file', help = "CSV output file.",
                        required = True, dest = 'out_file')
    parser.add_argument('-v', '--verbose',
                        help='increase output verbosity',
                        action='count', default=0)
    args = parser.parse_args()

    if args.verbose == 0:
        log_level = logging.INFO
    elif args.verbose == 1:
        log_level = logging.DEBUG
    else:
        log_level = logging.DEBUG
    
    logging.basicConfig(level=log_level,
                        format='%(levelname)-8s [%(asctime)s]  %(message)s',
                        datefmt="%y%m%d %H%M%S")

    logging.info("Program Started")
    with open(args.out_file, 'w') as out:
        out.write("Tool,Data,Technology,Individual,Chr,MeanCov,RawRealigned,WhDowns,WhDownsTimeSec")
        out.write(",Merge,MergeE,MergeM,MergeT,MergeN,MergeTimeSec,MergeMaxMemMB")
        out.write(",RndDowns,RndDownsSeed,RndDownsMaxCov,RndDownsTimeSec,RndDownsMaxMemMB")
        out.write(",FurtherMerging,Epsilon,Alpha,BalThr,BalRatio,IndelMode")
        out.write(",SwErrRatePerc,HamDistPerc,MecScore,PhasTimeSec,PhasMaxMemMB")
        out.write(",CleanFinish,FeasibleSoln")
        out.write("\n")
        # Hapchat
        num_inck = 0
        logging.info("Parsing Inck files")
        for df in os.listdir(args.inck_dir):
            if df.endswith(".bN_0.sum"):
                cleanfinish = 'yes' if os.stat(args.inck_dir+df).st_size else 'no'
                ds = df.rstrip().split(".")[:-1]
                dataset = ".".join(ds)
                # filter
                maxcovs = [15, 20, 25, 30]
                hs = ['h{}'.format(m) for m in maxcovs]
                downs = ['downs_s1_m{}'.format(m) for m in maxcovs]
                if ds[4] not in coverages :
                    continue
                if ds[6] not in hs + ['hN'] :
                    continue
                if ds[7] not in mergings + ['no_merging'] :
                    continue
                if ds[8] not in downs + ["no_downs"] :
                    continue

                whdataset = ".".join(ds[:7])+".wif"
                whdownstime = "NA"
                if ds[6] != "hN" :
                    whlfile = args.wif_dir + whdataset + ".log"
                    if not os.path.isfile(whlfile) :
                        logging.error("File not found: %s", whlfile)
                        exit()
                    whdownstime = str(getWhInfo(whlfile)["downs"])

                out.write("HapChat,")                
                num_inck += 1
                qual = getQUAL(args.inck_dir + dataset + '.diff')
                tfile = args.inck_dir + dataset + ".time"
                if not os.path.isfile(tfile):
                    logging.error("File not found: %s", tfile)
                    exit()
                info = getTimeMem(tfile)
                lfile = args.inck_dir + dataset + ".mec"
                if not os.path.isfile(lfile):
                    logging.error("File not found: %s", lfile)
                    exit()
                score = getMEC(lfile)
                # print(dataset + " " + str(score) +
                #       qual["ser"] + " " +
                #       str(info["time"]) + " " +
                #       str(info["mem"]))
                merge = "no,NA,NA,NA,NA,NA,NA"
                if(ds[7] != "no_merging"):
                    mfile = args.wif_dir + ".".join(ds[0:8]) + ".wif.time"
                    if not os.path.isfile(mfile):
                        logging.error("File not found: %s", mfile)
                        exit()
                    merge_info = getTimeMem(mfile)
                    mfields = ds[7].split('_')
                    merge = "yes" + "," + \
                            "0." + mfields[1][1:] + "," + \
                            "0." + mfields[2][1:] + "," + \
                            "1e+" + mfields[3][1:] + "," + \
                            "1e+" + mfields[4][1:] + "," + \
                            str(merge_info["time"]) + "," + \
                            str(merge_info["mem"])
                downs = "no,NA,NA,NA,NA"
                if(ds[8] != "no_downs"):
                    dfile = args.wif_dir + ".".join(ds[0:8])
                    # if(ds[7] != "no_merging"):
                    #     dfile += "." + ds[7]
                    dfile += ".wif.sample_" + "_".join(ds[8].split('_')[1:]) + ".time"
                    if not os.path.isfile(dfile):
                        logging.error("File not found: %s", dfile)
                        exit()
                    downs_info = getTimeMem(dfile)
                    dfields = ds[8].split('_')
                    downs = "yes" + "," + \
                            dfields[1][1:] + "," + \
                            dfields[2][1:] + "," + \
                            str(downs_info["time"]) + "," + \
                            str(downs_info["mem"])
                if(ds[9] != "no_merging"):
                    logging.error("Unrecognized further merging step!")
                    exit()
                out.write(",".join(ds[0:4]) + "," +
                          ds[4].replace("cov", "") + "," +
                          ds[5] + "," +
                          ds[6].replace("h", "") + "," +
                          whdownstime + "," +
                          merge + "," +
                          downs + "," +
                          "no" + "," +
                          "0." + ds[10].split("_")[0] + "," +
                          "0." + ds[10].split("_")[1] + "," +
                          ds[11].split("_")[0][1:] + "," +
                          ds[11].split("_")[1] + "," +
                          'NA,' +
                          qual["ser"] + "," +
                          qual["ham"] + "," +
                          str(score) + "," +
                          str(info["time"]) + "," +
                          str(info["mem"]) + ',' +
                          cleanfinish + ',yes')
                out.write("\n")
        logging.info("Parsed %d Inck files.", num_inck)
        # WhatsHap
        num_wh = 0
        logging.info("Parsing WhatsHap files")
        for df in os.listdir(args.wh_dir):
            if df.endswith(".sum"):
                cleanfinish = 'yes' if os.stat(args.wh_dir+df).st_size else 'no'
                ds = df.rstrip().split(".")[:-1]
                dataset = ".".join(ds)
                # filter
                maxcovs = [15, 20]
                hs = ['h{}'.format(m) for m in maxcovs]
                if ds[4] not in coverages :
                    continue
                if ds[6] not in hs :
                    continue

                out.write("WhatsHap,")
                num_wh += 1
                qual = getQUAL(args.wh_dir + dataset + '.diff')
                tfile = args.wh_dir + dataset + ".time"
                if not os.path.isfile(tfile):
                    logging.error("File not found: %s", tfile)
                    exit()
                info = getTimeMem(tfile)
                whlfile = args.wh_dir + dataset + ".log"
                if not os.path.isfile(whlfile):
                    logging.error("File not found: %s", whlfile)
                    exit()
                wh_info = getWhInfo(whlfile)
                msfile = args.wh_dir + dataset + ".mec"
                if not os.path.isfile(msfile):
                    logging.error("File not found: %s", msfile)
                    exit()
                score = getMEC(msfile)
                out.write(",".join(ds[0:4]) + "," +
                          ds[4].replace("cov", "") + "," +
                          ds[5] + "," +
                          ds[6].replace("h", "") + "," +
                          str(wh_info["downs"]) + "," +
                          "no,NA,NA,NA,NA,NA,NA" + "," +
                          "no,NA,NA,NA,NA" + "," +
                          "no,NA,NA" + "," +
                          "NA,NA,NA" + "," +
                          qual["ser"] + "," +
                          qual["ham"] + "," +
                          score + "," +
                          str(wh_info["time"]) + "," +
                          str(info["mem"]) + ',' +
                          cleanfinish + ',yes')
                out.write("\n")
                
        logging.info("Parsed %d WhatsHap files.", num_wh)

        # HapCol
        num_hc = 0
        logging.info('Parsing HapCol files')
        for df in os.listdir(args.hc_dir) :
            if df.endswith('.sum') :
                cleanfinish = 'yes' if os.stat(args.hc_dir+df).st_size else 'no'
                ds = df.rstrip().split('.')[:-1]
                dataset = '.'.join(ds)

                # filter
                maxcovs = [15, 20, 25, 30]
                hs = ['h{}'.format(m) for m in maxcovs]
                if ds[4] not in coverages :
                    continue
                if ds[6] not in hs :
                    continue

                whdataset = '.'.join(ds[:7]) + '.wif'
                whlfile = args.wif_dir + whdataset + '.log'
                if not os.path.isfile(whlfile) :
                    logging.error('File not found: %s', whlfile)
                    exit()
                whdownstime = str(getWhInfo(whlfile)['downs'])

                out.write('HapCol,')
                num_hc += 1
                qual = getQUAL(args.hc_dir + dataset + '.diff')
                tfile = args.hc_dir + dataset + '.time'
                if not os.path.isfile(tfile) :
                    logging.error('File not found: %s', tfile)
                    exit()
                info = getTimeMem(tfile)
                lfile = args.hc_dir + dataset + '.mec'
                if not os.path.isfile(lfile) :
                    logging.error('File not found: %s', lfile)
                    exit()
                score = getMEC(lfile)
                logfile = args.hc_dir + dataset + '.log'
                if not os.path.isfile(logfile) :
                    logging.error('File not found: %s', logfile)
                    exit()
                feasible = getFeasibility(logfile)
                out.write(','.join(ds[0:4]) + ',' +
                          ds[4].replace('cov', '') + ',' +
                          ds[5] + ',' +
                          ds[6].replace('h', '') + ',' +
                          whdownstime + ',' +
                          'no,NA,NA,NA,NA,NA,NA,' +
                          'no,NA,NA,NA,NA,' +
                          'no,NA,NA,' +
                          'NA,NA,NA,' +
                          qual['ser'] + ',' +
                          qual['ham'] + ',' +
                          score + ',' +
                          str(info['time']) + ',' +
                          str(info['mem']) + ',' +
                          cleanfinish + ',' + feasible)
                out.write('\n')

        logging.info('Parsed %d HapCol files.', num_hc)

        # HapCUT2
        num_hapcut2 = 0
        logging.info('Parsing HapCUT2 files')
        for df in os.listdir(args.hapcut2_dir) :
            if df.endswith('.sum') :
                cleanfinish = 'yes' if os.stat(args.hapcut2_dir+df).st_size else 'no'
                ds = df.rstrip().split('.')[:-1]
                dataset = '.'.join(ds)

                # filter
                if ds[4] not in coverages :
                    continue

                out.write('HapCUT2,')
                num_hapcut2 += 1

                qual = getQUAL(args.hapcut2_dir + dataset + '.diff')
                tfile = args.hapcut2_dir + dataset + '.time'
                if not os.path.isfile(tfile) :
                    logging.error('file not found: %s', tfile)
                    exit()
                info = getTimeMem(tfile)
                msfile = args.hapcut2_dir + dataset + '.mec'
                if not os.path.isfile(msfile) :
                    logging.error('file not found: %s', msfile)
                    exit()
                score = getMEC(msfile)
                indelmode = 'yes' if ds[10] == 'indels' else 'no'
                out.write(','.join(ds[0:4]) + ',' +
                          ds[4].replace('cov','') + ',' +
                          ds[5] + ',' +
                          ds[6].replace('h','') + ',' +
                          'NA,' +
                          'no,NA,NA,NA,NA,NA,NA,' +
                          'no,NA,NA,NA,NA,' +
                          'no,NA,NA,' +
                          'NA,NA,' +
                          indelmode + ',' +
                          qual['ser'] + ',' + qual['ham'] + ',' +
                          str(score) + ',' +
                          str(info['time']) + ',' + str(info['mem']) + ',' +
                          cleanfinish + ',yes')
                out.write('\n')
        logging.info('Parsed %d HapCUT2 files.', num_hapcut2)

    logging.info("Program Finshed")

if __name__ == "__main__":
    main()
