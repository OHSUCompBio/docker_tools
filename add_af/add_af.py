#!/usr/bin/env python3

# Add an allele frequency and total depth entry, useful for Pindel VCF's.
# USAGE: python add_af.py <input VCF> <output VCF>
# TODO: Split filtering functionality out of this script in to with another tool
# or a separate script.  This script was originally put together for one specific
# use case (Pindel) but now may be more widely applicable.

import argparse
VERSION = '0.4.6'


def create_header_entry(attrib):
    """
    This goes in the header of the VCF, and looks like this:
    ##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Allele depth, how many reads support this allele">
    Modify to allow for different Number, Type, and Description field to be passed.
    """

    header = None
    if attrib == "AF":
        header = "##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Variant allele frequency\">\n"
    elif attrib == "DP":
        header = "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth at this locus.\">\n"
    return header


def find_index(splitting, delim, to_find):
    """
    Find the index of the field in the FORMAT column you care about.  Usually will be AD.
    """

    for entry in splitting.split(delim):
        if entry == to_find:
            curr_index = splitting.split(delim).index(to_find)
            return curr_index

    return None


def calc_af(ref, alt):
    """
    Calculate the allele frequency and format the output.
    """
    
    if ref + alt != 0:
        af = alt / (ref + alt + 0.0)
        if af != 0.0:
            af = "{:.7f}".format(alt / (ref + alt + 0.0)).rstrip('0')
    else:
        af = 0.0

    return af


def calc_dp(ref_count, alt_count):
    """
    Pindel VCF's do not contain a DP value in the SAMPLE field.  Calculate it here.
    """

    return str(ref_count + alt_count)


def check_genos(geno1, geno2):
    """
    If the genotypes are the same, do not write to output.
    """
    if geno1 == geno2:
        return True

    return False


def compare_af(af1, af2, margin=0.2, thresh=0.1):
    """
    Make a comparison between the allele frequencies, and remove variants that don't
    differ by at least margin and aren't greater than thresh.
    The idea is to keep those with af close to 0, so we don't lose low af somtic variants.
    """

    if af1 >= af2:
        diff = af1-af2
    else:
        diff = af2-af1    

    if diff <= margin:
        if af1 <= thresh or af2 <= thresh:
            return False
    elif diff > margin:
        return False

    return True


def check_format(frmt, field):
    """
    Check to see if the field you want to add is already there.
    """
    return field in frmt


def gather_header_terms(header_line):
    """
    Gather everything from the header, and determine which fields are already there.
    Can look like:
    ##INFO=<ID=NTLEN,Number=.,Type=Integer,Description="Number of bases inserted in place of deleted code">
    ##FORMAT=<ID=PL,Number=3,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined
    in the VCF specification">
    """

    return header_line.split('=')[2].split(',')[0]


def supply_args():
    """
    Populate args.                                         
    https://docs.python.org/2.7/library/argparse.html
    """

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('infile', help='Input VCF')
    parser.add_argument('outfile', help='Output VCF')
    parser.add_argument('format_label', type=str, help='Label to search for in the format field, such as AD or DPR.')
    parser.add_argument('--margin', type=float, help='Allele frequencies with spread lower than this value will be '
                                                     'considered the same if they are also not below the threshold '
                                                     'defined below.')
    parser.add_argument('--thresh', type=float, help='If allele frequencies are below this value, this record will '
                                                     'be included regardless.')
    parser.add_argument('--filtering', action="store_true", help='Should additional filtering be performed based on '
                                                                 '--margin and --thresh')
    parser.add_argument('--vaf', type=float, help='Provide a VAF cutoff for '
                                                  'final output.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()

    return args


def main():
    args = supply_args()
    handle_vcf = open(args.infile, 'r')
    handle_out = open(args.outfile, 'w')
    margin = args.margin
    thresh = args.thresh
    # We will usually target ':' delimeters, as that is what the VCF uses.
    delim = ':'
    # Generally should be looking for AD's (allele depth) here.
    to_find = args.format_label
    header_terms = []

    with handle_vcf as myvcf:
        for variant in myvcf:
            if variant[0] != "#":

                curr_index = None
                split_variant = variant.rstrip('\n').split('\t')
                frmt = split_variant[8]

                if not curr_index:
                    curr_index = find_index(frmt, delim, to_find)
                if not curr_index:
                    continue

                if len(split_variant) == 11 and split_variant[10] == '':
                    split_variant = split_variant[:10]

                if len(split_variant) == 11:
                    normal = split_variant[9]
                    tumor = split_variant[10]
                    normal_ref = int(normal.split(':')[curr_index].split(',')[0])
                    normal_alt = int(normal.split(':')[curr_index].split(',')[1])
                    tumor_ref = int(tumor.split(':')[curr_index].split(',')[0])
                    tumor_alt = int(tumor.split(':')[curr_index].split(',')[1])

                    normal_af = calc_af(normal_ref, normal_alt)
                    tumor_af = calc_af(tumor_ref, tumor_alt)

                    if not check_format(frmt, "AF"):
                        # Create FORMAT and SAMPLE strings with AF included.
                        frmt = frmt + delim + "AF"
                        normal = normal + delim + str(normal_af)
                        tumor = tumor + delim + str(tumor_af)

                    if not check_format(frmt, "DP"):
                        # Create FORMAT and SAMPLE strings with DP included.
                        frmt = frmt + delim + "DP"
                        normal = normal + delim + calc_dp(normal_ref, normal_alt)
                        tumor = tumor + delim + calc_dp(tumor_ref, tumor_alt)

                    normal_geno = normal.split(':')[0]
                    tumor_geno = tumor.split(':')[0]

                    if args.filtering:
                        if not check_genos(normal_geno, tumor_geno) and \
                                not compare_af(float(normal_af), float(tumor_af), margin, thresh) \
                                and int(calc_dp(tumor_ref, tumor_alt)) > 2 and normal_af < .05:
                            handle_out.write('\t'.join(['\t'.join(split_variant[:8]), frmt, normal, tumor, '\n']))
                    else:
                        handle_out.write('\t'.join(['\t'.join(split_variant[:8]), frmt, normal, tumor, '\n']))

                # For situation where there is not a matched normal.
                elif len(split_variant) == 10:
                    tumor = split_variant[9]
                    if to_find == "DPR":
                        curr_index = find_index(frmt, delim, to_find)
                        dpr0 = int(tumor.split(':')[curr_index].split(',')[0])
                        try:
                            dpr1 = int(tumor.split(':')[curr_index].split(',')[1])
                        except IndexError:
                            dpr1 = 0
                        tumor_ref = dpr0 - dpr1
                        tumor_alt = dpr1
                        tumor_af = calc_af(tumor_ref, tumor_alt)
                    else:
                        try:
                            tumor_ref = int(tumor.split(':')[curr_index].split(',')[0])
                            tumor_alt = int(tumor.split(':')[curr_index].split(',')[1])
                            tumor_af = calc_af(tumor_ref, tumor_alt)
                        except IndexError:
                            tumor_af = 0

                    if not check_format(frmt, "AF"):
                        # Create FORMAT and SAMPLE strings with AF included.
                        frmt = frmt + delim + "AF"
                        tumor = tumor + delim + str(tumor_af)

                    if not check_format(frmt, "DP"):
                        # Create FORMAT and SAMPLE strings with DP included.
                        frmt = frmt + delim + "DP"
                        tumor = tumor + delim + calc_dp(tumor_ref, tumor_alt)

                    if args.filtering:
                        raise Exception("This option is not valid for tumor-only VCF's.")
                    elif args.vaf:
                        if float(tumor_af) >= args.vaf:
                            handle_out.write('\t'.join(['\t'.join(
                                split_variant[:8]), frmt, tumor, '\n']))
                    elif not args.vaf:
                        handle_out.write('\t'.join(['\t'.join(
                            split_variant[:8]), frmt, tumor, '\n']))

            else:
                if "CHROM" not in variant:
                    handle_out.write(variant)
                    if "FORMAT" in variant:
                        header_terms.append(gather_header_terms(variant))
                else:
                    if "AF" not in header_terms:
                        handle_out.write(create_header_entry("AF"))
                    if "DP" not in header_terms:
                        handle_out.write(create_header_entry("DP"))
                    handle_out.write(variant)


if __name__ == "__main__":
    main()
