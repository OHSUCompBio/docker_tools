#!/usr/bin/env python

# DESCRIPTION: Given a panel of normals VCF, as created by
# panel_of_normals.py, annotate or remove entries in a sample VCF.
# USAGE: use_pon.py <pon> <infile> <outfile>
# CODED BY: John Letaw

from __future__ import print_function
import argparse

VERSION = '0.4.1'

def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('pon', help='Panel of normals VCF.')
    parser.add_argument('infile', help='Input VCF.')
    parser.add_argument('hotspots', help='Hotspots VCF.')
    parser.add_argument('outfile', help='Output VCF.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def parse_pon(args):
    """
    Run through the PON and grab what we need.
    :return:
    """
    drop_list = []
    annot_list = []
    all_pon = {}
    with open(args.pon, 'rU') as pon:
        for line in pon:
            if not line.startswith('#'):
                line = line.rstrip('\n').split('\t')
                seen = line[7].split(';')[5].split('=')[1]
                cosmic = line[7].split(';')[4].split('=')[1]
                avg_af = float(line[7].split(';')[0].split('=')[1])
                stdev_af = float(line[7].split(';')[1].split('=')[1])
                uniq_key = (line[0], line[1], line[3], line[4])
                if uniq_key not in all_pon:
                    all_pon[uniq_key] = (avg_af, stdev_af)
#                if avg_af < 0.15 and stdev_af < 0.06:
                    
                # if int(seen) >= 3 and cosmic == 'F':
                #     drop_list.append((line[0], line[1], line[3], line[4]))
                # nixed
                # elif int(seen) >= 5 and cosmic == 'T':
                #     drop_list.append((line[0], line[1], line[3], line[4]))
                # else:
                #     annot_list.append((line[0], line[1], line[3], line[4]))

    return all_pon


def replace_filter(filt, anno):
    """
    Remove the dot if it's there, otherwise just append something with a
    semicolon.
    :return:
    """
    if filt == '.' or filt == 'PASS':
        return anno
    else:
        return ';'.join([filt, anno])

def write_new_vcf(args, all_pon, hotspots):
    """
    Write a new VCF with the panel_of_normals tag in FILTER section.
    :return:
    """
    outfile = open(args.outfile, 'w')
    check = True
    with open(args.infile, 'rU') as infile:
        for line in infile:
            if not line.startswith('#'):
                sline = line.rstrip('\n').split('\t')
                labels = sline[8].split(':')
                info = sline[9].split(':')
                vaf = float(info[labels.index('AF')])
                uniq_key = (sline[0], sline[1], sline[3], sline[4])
                write_me = False
                if uniq_key in all_pon:
                    avg_af = float(all_pon[uniq_key][0])
                    stdev_af = float(all_pon[uniq_key][1])

                    if stdev_af == 0.0:
                        stdev_af = 0.001
                    if avg_af < 0.2 and stdev_af < 0.06:
                        stdev_af = max(stdev_af, .01)
                        if vaf >= (3*stdev_af + avg_af):
                            sline[6] = replace_filter(sline[6], 'PON_OV')
                            write_me = True
                    else:
                        sline[6] = replace_filter(sline[6], 'PON')
                        write_me = True

                else:
                    write_me = True

                if write_me:
                    outfile.write('\t'.join(sline))
                    outfile.write('\n')

                    # if uniq_key in hotspots:
                    #     sline[6] = replace_filter(sline[6], 'PON_hotspot')
                    #     outfile.write('\t'.join(sline))
                    #     outfile.write('\n')

            else:
                outfile.write(line)
                if line.startswith("##FILTER") and check:
                    outfile.write("##FILTER=<ID=PON,Description=\"Variant found in panel of normals in 3 or more samples.\">\n")
                    outfile.write("##FILTER=<ID=PON_OV,Description=\"Variant found in panel of normals at low frequency, but call is at significantly higher frequency.\">\n")
                    # outfile.write("##FILTER=<ID=PON_hotspot,Description=\"Variant would have been filtered but appears in the hotspot list.\">\n")
                    check = False

    outfile.close()


def hotspots_grab(filename):
    """

    :param filename:
    :return:
    """
    hotspots = []
    with open(filename, 'rU') as hots:
        for line in hots:
            if not line.startswith('#'):
                line = line.rstrip('\n').split('\t')
                chrom = line[0]
                pos = line[1]
                ref = line[3]
                alt = line[4]
                uniq_key = (chrom, pos, ref, alt)
                if uniq_key not in hotspots:
                    hotspots.append(uniq_key)
    return hotspots


def main():

    args = supply_args()
    hotspots = hotspots_grab(args.hotspots)
    all_pon = parse_pon(args)
    write_new_vcf(args, all_pon, hotspots)


if __name__ == "__main__":
    main()
