# -*- coding: utf-8 *-*

import sys
import os
import errno
import time
import pprint

# global configuration
import vprimer.glv as glv
import vprimer.utils as utl

from vprimer.logging_config import LogConf
log = LogConf.open_log(__name__)

import vcfpy

from vprimer.allele_select import AlleleSelect

class Variant(object):

    def __init__(self):
        pass


    def pick_variant(self):
        """
        """

        proc_name = "variant"
        log.info("-------------------------------")
        log.info("Start processing {}\n".format(proc_name))

        # stop, action, gothrough
        ret_status = utl.decide_action_stop(proc_name)

        if ret_status == "stop":
            msg = "STOP. "
            msg += "Current process \'{}\' ".format(proc_name)
            msg += "has exceeded the User-specified stop point "
            msg += "\'{}', ".format(glv.conf.stop)
            msg += "so stop program. exit."
            log.info(msg)
            sys.exit(1)


        elif ret_status == "gothrough":
            msg = "SKIP \'{}\' proc, ".format(proc_name)
            msg += "glv.conf.progress = {}, ".format(glv.conf.progress)
            msg += "glv.conf.stop = {}, ".format(glv.conf.stop)
            msg += "so skip program."
            log.info(msg)
            return


        # open vcf through vcfpy
        reader = vcfpy.Reader.from_path(glv.conf.vcf_file_path)

        # for each distinguish_groups
        for proc_cnt, distin_dict in enumerate(glv.outlist.distin_files, 1):

            reg = distin_dict['region']
            # vcf_ittr for each distinguish groups
            vcf_ittr = reader.fetch(glv.conf.regions_dict[reg]['reg'])
            self._iterate_vcf(vcf_ittr, distin_dict, proc_cnt)


    def _iterate_vcf(self, vcf_ittr, distin_dict, proc_cnt):
        """
        """

        # basic informations
        gr_list = [distin_dict[0], distin_dict[1]]

        reg = distin_dict['region']
        reg_dict = glv.conf.regions_dict[reg]
        pick_mode = distin_dict['pick_mode']
        indel_size = distin_dict['indel_size']
        min_indel_len, max_indel_len = \
            [int(i) for i in indel_size.split('-')]

        # At first, we check difference of genotype between two sample
        # that described at the beginning of each group
        top_smpl_list = [
            glv.conf.group_members_dict[gr_list[0]][0],
            glv.conf.group_members_dict[gr_list[1]][0]]

        # logging current target
        utl.print_distin_info("variant", distin_dict, proc_cnt)

        start = time.time()

        # File name to export variant
        out_txt_file = distin_dict['variant']['out_path']

        utl.save_to_tmpfile(out_txt_file)

        #------------------------------------------------------
        # To add an allele_int column for all sample
        # Members of the specified group come first
        # gr0:s1 g0:s2 g0:s3 g1:s4 g1:s5 g1:s6 s7 s8 s9 s10

        sample_nickname_ordered_list, \
        sample_fullname_ordered_list = \
            utl.get_ordered_sample_list(gr_list)

        sample_added_header = "{}\t{}".format(
            distin_dict['variant']['hdr_text'],
            "\t".join(sample_nickname_ordered_list))

        # Can I parallelize here?
        with open(out_txt_file, mode='a') as f:

            # write sample added header
            f.write("{}\n".format(sample_added_header))

            # access to vcf using iterater
            for record in vcf_ittr:

                # 1. Skip same GT between top two sample
                if self._skip_same_GT_between_top2sample(
                    record, top_smpl_list) > 0:
                    continue

                # 2. Check GT in your own group
                if self._skip_different_GT_in_own_group(
                    record, top_smpl_list, gr_list) > 0:
                    continue

                # 3. Select different allele combination among 2x2 allele
                asel = AlleleSelect(min_indel_len, max_indel_len)
                asel.select_diff_allele(record, top_smpl_list, gr_list)

                # from record, construct allele_int of the member
                # who is paying attention
                allele_int_line = ""

                # 4. Save variant information as text file
                for var_type, line in zip(asel.var_types, asel.lines):
                    if utl.is_my_pick_mode(
                        var_type, distin_dict['pick_mode']) == True:

                        # make allele_int line
                        if allele_int_line == "":
                                #self._get_ai_line(
                            allele_int_line = \
                                self._get_allele_line(
                                    record, sample_fullname_ordered_list)

                        # add allele line
                        f.write("{}\t{}\n".format(line, allele_int_line))

        log.info("variant {} > {}.txt\n".format(
            utl.elapsed_time(time.time(), start),
            distin_dict['variant']['base_nam']))


    def _get_ai_line(self, record, sample_fullname_list):
        '''
        '''
        #line = [record.CHROM, record.POS, record.REF]
        #alt_list = [alt.value for alt in record.ALT]
        #line += [",".join(alt_list)]
        line = list()
        line += [AlleleSelect.allele_int("{}/{}".format(
            record.call_for_sample[fn].gt_alleles[0],
            record.call_for_sample[fn].gt_alleles[1]
            ), "int") for fn in sample_fullname_list]

        line_str = '\t'.join(map(str, line))
        return line_str


    def _get_allele_line(self, record, sample_fullname_list):
        '''
        '''
        #line = [record.CHROM, record.POS, record.REF]
        #alt_list = [alt.value for alt in record.ALT]
        #line += [",".join(alt_list)]
        line = list()
        line += [AlleleSelect.allele_int("{}/{}".format(
            record.call_for_sample[fn].gt_alleles[0],
            record.call_for_sample[fn].gt_alleles[1]
            ), "allele") for fn in sample_fullname_list]

        line_str = '\t'.join(map(str, line))
        return line_str

    def _skip_different_GT_in_own_group(self, record, tsl, gr_list):

        skip = glv.SKIP_DONT_SKIP

        # check twice, group0, and group1
        for gr_no in range(2):
            # pick sample name belong to a group
            for (sample_no, sample_name) in enumerate(
                glv.conf.group_members_dict[gr_list[gr_no]]):

                if sample_no == 0:
                    continue    # self

                sample0 = tsl[gr_no]
                sample1 = sample_name

                # もし、サンプル間でvariantが見つかった場合は、
                s0_0, s0_1, s1_0, s1_1 = \
                    AlleleSelect.record_call_for_sample(
                        record, sample0, sample1)

                # compare alleles with first sample
                if s0_0 == s1_0 and s0_1 == s1_1:

                    #log.debug("SKIP_SAME_HOMO {},({}){} {}{}/{}{}".format(
                    #    gr_list[gr_no],
                    #    sample_no, sample_name,
                    #    record.call_for_sample[tsl[gr_no]].gt_alleles[0],
                    #    record.call_for_sample[tsl[gr_no]].gt_alleles[1],
                    #    record.call_for_sample[sample_name].gt_alleles[0],
                    #    record.call_for_sample[sample_name].gt_alleles[1]))
                    pass

                else:
                    skip = glv.SKIP_DIFF_INGROUP
                    #log.debug("SKIP_SAME_HOMO {},({}){} {}{}/{}{}".format(
                    #    gr_list[gr_no],
                    #    sample_no, sample_name,
                    #    record.call_for_sample[tsl[gr_no]].gt_alleles[0],
                    #    record.call_for_sample[tsl[gr_no]].gt_alleles[1],
                    #    record.call_for_sample[sample_name].gt_alleles[0],
                    #    record.call_for_sample[sample_name].gt_alleles[1]))
                    return skip

        return skip


    def _skip_same_GT_between_top2sample(self, record, tsl):

        # for REF 20200708
        sample0 = tsl[0]
        sample1 = tsl[1]

        s0_0, s0_1, s1_0, s1_1 = \
            AlleleSelect.record_call_for_sample(record, sample0, sample1)

        skip = glv.SKIP_DONT_SKIP

        # ./. only 0
        if utl.is_None(s0_0, s0_1, s1_0, s1_1):
            skip = glv.SKIP_None
            #log.debug("SKIP_None {}{}/{}{}".format(s0_0,s0_1,s1_0,s1_1))
            return skip

        # same homo: AA,AA
        if utl.is_same_homo(s0_0, s0_1, s1_0, s1_1):
            skip = glv.SKIP_SAME_HOMO
            #log.debug("SKIP_SAME_HOMO {}{}/{}{}".format(s0_0,s0_1,s1_0,s1_1))
            return skip

        # same hetero: AB,AB
        if utl.is_same_hetero(s0_0, s0_1, s1_0, s1_1):
            skip = glv.SKIP_SAME_HETERO
            #log.debug("SKIP_SAME_HETERO {}{}/{}{}".format(
            #    s0_0,s0_1,s1_0,s1_1))
            return skip

        return skip


    def print_allele(self):

        ''' When show_genotype is specified, the genotype of the specified
        regions and members are output to a file.
        '''

        proc_name = "genotype"
        log.info("-------------------------------")
        log.info("Start processing {}\n".format(proc_name))

        # header
        header = list()
        header += ["CHROM", "POS", "REF", "ALT", "Rlen", "Alen", "diff"]
        header += glv.conf.group_members_dict['all']

        # reader
        reader = vcfpy.Reader.from_path(glv.conf.vcf_file_path)

        total_cnt = len(glv.conf.region_name_list)

        # Save to file for each region
        for proc_cnt, region_name in enumerate(glv.conf.region_name_list, 1):

            region = glv.conf.regions_dict[region_name]['reg']

            # Create a list of fullname for the specified members
            sample_fullname_list = list()
            for nickname in glv.conf.group_members_dict['all']:
                sample_fullname_list.append(utl.get_fullname(nickname))

            # if group priority
            #sample_fullname_list = \
            #    utl.get_sample_list_from_groupname(
            #        group_list, "fullname")

            # out file name
            outf_pref = "005_genotype"
            basename = "{}~{}~{}".format(
                outf_pref, region_name, glv.conf.show_genotype)
            out_file_path = "{}/{}.txt".format(
                glv.conf.out_dir_path, basename)

            # backup
            utl.save_to_tmpfile(out_file_path)

            log.info("")
            log.info("{} / {}, {}({}) > {}".format(
                proc_cnt, total_cnt, region_name, region, out_file_path ))

            start = time.time()
            with open(out_file_path, mode='w') as f:

                f.write("{}\n".format('\t'.join(map(str, header))))

                vcf_ittr = reader.fetch(region)
                for record in vcf_ittr:

                    # Main informations
                    line = [record.CHROM, record.POS, record.REF]
                    alt_list = [alt.value for alt in record.ALT]
                    line += [",".join(alt_list)]

                    # variant length and diff
                    len_ref = len(record.REF)
                    lens_alt_list = list()
                    for alt in alt_list:
                        lens_alt_list.append(len(alt))

                    diff_len = abs(len_ref - lens_alt_list[0])
                    lens_alt = ",".join(map(str, lens_alt_list))

                    line += [len_ref]
                    line += [lens_alt]
                    line += [diff_len]

                    line += [AlleleSelect.allele_int(
                        "{}/{}".format(
                            record.call_for_sample[fn].gt_alleles[0],
                            record.call_for_sample[fn].gt_alleles[1]
                        ), glv.conf.show_genotype
                        ) for fn in sample_fullname_list]

                    f.write("{}\n".format('\t'.join(map(str, line))))

            log.info("genotype {} > {}.txt\n".format(
                utl.elapsed_time(time.time(), start),
                out_file_path))


    def print_all_allele_int(self):
        '''
        '''

        header = list()
        header += ["CHROM", "POS", "REF", "ALT"]
        header += glv.conf.vcf_sample_nickname_list
        #print("#{}".format("\t".join(header)))

        reader = vcfpy.Reader.from_path(glv.conf.vcf_file_path)

        # all chromosome region
        for region in glv.conf.ref_fasta_chrom_region_list:

            # for members full name
            glv.conf.vcf_sample_fullname_list

            vcf_ittr = reader.fetch(region)
            for record in vcf_ittr:

                line = [record.CHROM, record.POS, record.REF]
                alt_list = [alt.value for alt in record.ALT]
                line += [",".join(alt_list)]
                line += [AlleleSelect.allele_int("{}/{}".format(
                    record.call_for_sample[fn].gt_alleles[0],
                    record.call_for_sample[fn].gt_alleles[1]
                    )) for fn in glv.conf.vcf_sample_fullname_list]

                #print('\t'.join(map(str, line)))


    def print_allele_int(self):
        '''
        '''

        header = list()
        header += ["CHROM", "POS", "REF", "ALT"]
        header += glv.conf.vcf_sample_nickname_list
        #print("#{}".format("\t".join(header)))

        reader = vcfpy.Reader.from_path(glv.conf.vcf_file_path)
        for distin_dict in glv.outlist.distin_files:

            # for region
            region_name = distin_dict['region']
            region = glv.conf.regions_dict[region_name]['reg']

            # for members full name
            group_list = [distin_dict[0], distin_dict[1]]
            sample_fullname_list = \
                utl.get_sample_list_from_groupname(
                    group_list, "fullname")

            vcf_ittr = reader.fetch(region)
            for record in vcf_ittr:

                line = [record.CHROM, record.POS, record.REF]
                alt_list = [alt.value for alt in record.ALT]
                line += [",".join(alt_list)]
                line += [AlleleSelect.allele_int("{}/{}".format(
                    record.call_for_sample[fn].gt_alleles[0],
                    record.call_for_sample[fn].gt_alleles[1]
                    )) for fn in sample_fullname_list]

                #print('\t'.join(map(str, line)))


