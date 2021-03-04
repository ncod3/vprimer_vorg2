# -*- coding: utf-8 -*-

import sys
import os
import errno
import time
import re

# global configuration
import vprimer.glv as glv
import vprimer.utils as utl

from vprimer.logging_config import LogConf
log = LogConf.open_log(__name__)

import pandas as pd
import vcfpy
import subprocess as sbp

from joblib import Parallel, delayed
#import dill as pickle

from vprimer.product import Product
from vprimer.allele_select import AlleleSelect
from vprimer.eval_variant import EvalVariant
from vprimer.variant import Variant
from vprimer.inout_primer3 import InoutPrimer3
from vprimer.blast import Blast

class Primer(object):

    def __init__(self):

        pass

    def construct_primer(self):

        proc_name = "primer"
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


        # for each distinguish_groups
        for proc_cnt, distin_dict in enumerate(glv.outlist.distin_files, 1):

            # logging current target
            utl.print_distin_info("primer", distin_dict, proc_cnt)

            marker_file = distin_dict['marker']['out_path']
            df_distin = pd.read_csv(
                marker_file, sep='\t', header=0, index_col=None)

            out_txt_file = distin_dict['primer']['out_path']
            utl.save_to_tmpfile(out_txt_file)

            with open(out_txt_file, mode='a') as f:
                # write header
                #f.write("{}\n".format(distin_dict['primer']['hdr_text']))

                start = time.time()

                if glv.conf.parallel == True:
                    log.info(
                        "do Parallel cpu {}, parallel {} blast {}".format(
                            glv.conf.thread,
                            glv.conf.parallel_blast_cnt,
                            glv.conf.blast_num_threads))

                    Parallel(
                        n_jobs=glv.conf.parallel_blast_cnt,
                        backend="threading")(
                        [
                            delayed(self._loop_primer3_check_blast) \
                                (distin_dict, marker_df_row, f) \
                                for marker_df_row in df_distin.itertuples()
                        ]
                    )

                else:
                    log.info("do Serial cpu {} / serial {} blast {}".format(
                        glv.conf.thread,
                        1,
                        glv.conf.blast_num_threads))

                    for marker_df_row in df_distin.itertuples():

                        self._loop_primer3_check_blast(
                            distin_dict, marker_df_row, f)

            utl.sort_file(
                'primer', distin_dict, out_txt_file,
                'chrom', 'pos', 'try_cnt', 'number')

            log.info("primer {} > {}.txt\n".format(
                utl.elapsed_time(time.time(), start),
                distin_dict['primer']['base_nam']))


    def _loop_primer3_check_blast(self, distin_dict, marker_df_row, f):

        prinfo = PrimerInfo()
        prinfo.prepare_from_marker_file(distin_dict, marker_df_row)

        prinfo.get_excluded_region()
        #print("stop")
        #sys.exit(1)


        prinfo.make_p3_info()
        blast_check_result_list = list()

        line = ''
        try_cnt = 0
        blast_check = ''

        while True:

            #log.debug("<try_cnt {}> SEQUENCE_EXCLUDED_REGION={}".format(
            #    try_cnt,
            #    iopr3.get_sequence_excluded_region()))
            try_cnt += 1

            # primer3
            self._do_primer3_pipe(prinfo.iopr3)

            #log.debug("")
            #log.debug("{}".format(prinfo.iopr3.get_primer3_out()))

            #log.debug("\nP3_COMMENT={}".format(iopr3.get_p3_comment()))

            # break if cann't generate primer
            if prinfo.iopr3.PRIMER_ERROR != '':
                #log.info("try_cnt {} {} PRIMER_ERROR={}".format(
                #    try_cnt,
                #    prinfo.iopr3.get_sequence_id(),
                #    prinfo.iopr3.PRIMER_ERROR))
                break

            elif prinfo.iopr3.PRIMER_PAIR_NUM_RETURNED == 0:
                #log.info("try_cnt {} {} PRIMER_PAIR_NUM_RETURNED={}".format(
                #    try_cnt,
                #    prinfo.iopr3.get_sequence_id(),
                #    prinfo.iopr3.PRIMER_PAIR_NUM_RETURNED))
                break

            # calc abs pos
            abs_left_stt, abs_left_end, \
            abs_right_stt, abs_right_end = \
                self._get_primer_pos_info(prinfo)

            # make fasta id
            left_fasta_id = self._make_primer_name(
                prinfo.chrom, abs_left_stt, abs_left_end, "plus")
            right_fasta_id = self._make_primer_name(
                prinfo.chrom, abs_right_stt, abs_right_end, "minus")

            prinfo.iopr3.set_primer_name(left_fasta_id, right_fasta_id)

            # search my blastn-short
            blast_check_result_list = \
                Blast.primer_blast_check(
                    left_fasta_id, right_fasta_id,
                    prinfo.iopr3.get_primer_left_seq(),
                    prinfo.iopr3.get_primer_right_seq())

            # break if complete
            if len(blast_check_result_list) == 0:

                #log.info("complete try_cnt {} {}".format(
                #    try_cnt, prinfo.iopr3.get_sequence_id()))

                #------------------------------
                complete = 1
                line, blast_check = self._primer_complete_to_line(
                    complete, blast_check_result_list, prinfo, try_cnt)
                f.write('{}\n'.format(line))
                break

            else:
                # go to next chance to add the primer pos to ex
                prinfo.iopr3.add_ex_region(prinfo.iopr3.get_primer_region())

                #log.info("next try_cnt {} {}".format(
                #    try_cnt,
                #    prinfo.iopr3.get_sequence_excluded_region()))

                #------------------------------
                complete = 0
                line, blast_check = self._primer_complete_to_line(
                    complete, blast_check_result_list, prinfo, try_cnt)
                f.write('{}\n'.format(line))


        if len(blast_check_result_list) != 0:
            log.info("skipped try_cnt {} {} {}".format(
                try_cnt,
                prinfo.iopr3.get_sequence_id(),
                blast_check))
           # log.info(iopr3.p3_out)

#        if line != '':
#            f.write('{}\n'.format(line))
#            # if you need
#            f.flush()


    def _primer_complete_to_line(
        self, complete, blast_check_result_list, prinfo, try_cnt):

        # blast_check
        blast_check = "-"
        add_cnt = len(blast_check_result_list) - 1
        if add_cnt == -1:
            pass
        elif add_cnt == 0:
            blast_check = "{}".format(blast_check_result_list[0])
        else:
            blast_check = "{}(+{})".format(
                blast_check_result_list[0], add_cnt)


        l_list = list()

        # to primer out file
        l_list += [prinfo.marker_id]
        # --------------
        l_list += [prinfo.chrom]
        l_list += [prinfo.pos]
        l_list += [prinfo.targ_grp]
        l_list += [prinfo.targ_ano]
        l_list += [prinfo.vseq_gno_str]
        l_list += [prinfo.gts_segr_lens]
        l_list += [prinfo.var_type]
        l_list += [prinfo.set_enz_cnt]
        # --------------
        l_list += [prinfo.marker_info]
        l_list += [prinfo.vseq_lens_ano_str]
        # --------------
        l_list += [prinfo.g0_seq_target_len]
        l_list += [prinfo.g0_seq_target]
        l_list += [prinfo.g1_seq_target_len]
        l_list += [prinfo.g1_seq_target]

        l_list += [prinfo.seq_template_ref_len]
        l_list += [prinfo.seq_template_ref_abs_pos]
        l_list += [prinfo.seq_template_ref_rel_pos]
        # --------------
        l_list += [try_cnt]
        l_list += [complete]
        l_list += [blast_check]
        # --------------

        l_list += [prinfo.iopr3.get_primer_product_size()]

        l_list += [prinfo.iopr3.get_primer_left()]
        l_list += [prinfo.iopr3.get_primer_left_id()]
        l_list += [prinfo.iopr3.get_primer_left_seq()]

        l_list += [prinfo.iopr3.get_primer_right()]
        l_list += [prinfo.iopr3.get_primer_right_id()]
        l_list += [prinfo.iopr3.get_primer_right_seq()]

        l_list += [prinfo.SEQUENCE_TARGET]
        # 
        l_list += [prinfo.iopr3.get_sequence_excluded_region()]
        l_list += [prinfo.seq_template_ref]

        return '\t'.join(map(str, l_list)), blast_check


    def _make_primer_name(
            self, chrom, abs_primer_stt_pos, abs_primer_end_pos, strand):

        # {NC_028450.1}44676.44700.plus
        #primer_name = "{{{}}}{}.{}.{}".format(
        #    chrom,
        #    abs_primer_stt_pos,
        #    abs_primer_end_pos,
        #    strand)

        # NC_028450.1:44676-44700:plus
        primer_name = "{}:{}-{}:{}".format(
            chrom,
            abs_primer_stt_pos,
            abs_primer_end_pos,
            strand)

        return primer_name


    def _do_primer3_pipe(self, iopr3):

        # exec primer3 through pipe
        primer3_in = iopr3.get_p3_input()

        primer3_out_p = sbp.Popen(
            ['primer3_core'],
            stdin=sbp.PIPE,
            stdout=sbp.PIPE)

        primer3_out = primer3_out_p.communicate(
            primer3_in.encode())[0].decode()

        iopr3.set_primer3_out(primer3_out)


    def _get_primer_pos_info(self, prinfo):
        '''
        '''

        # This value brings PRIMER_LEFT_0 value directly.
        # 353,25
        PRIMER_LEFT_0_stt, PRIMER_LEFT_0_len = \
            prinfo.iopr3.get_primer_left_info()

        # For PRIMER_LEFT_0, the starting point is reported
        # in the 5'-3'direction as usual.
        # Therefore, the absolute position is calculated
        # based on the length as usual.
        abs_PRIMER_LEFT_0_stt = self._get_abspos(
            PRIMER_LEFT_0_stt, prinfo.abs_frag_pad_pre_stt)
        abs_PRIMER_LEFT_0_end = abs_PRIMER_LEFT_0_stt + PRIMER_LEFT_0_len - 1

        # This value brings PRIMER_RIGHT_0 value directly.
        # 567,25
        # On the other hand, for PRIMER_RIGHT_0,
        # the position on the 3'side is reported.
        PRIMER_RIGHT_0_end, PRIMER_RIGHT_0_len = \
            prinfo.iopr3.get_primer_right_info()

        # Therefore, the starting point is obtained using the length.
        PRIMER_RIGHT_0_stt = PRIMER_RIGHT_0_end - PRIMER_RIGHT_0_len + 1
        # After that, just convert to absolute position.
        abs_PRIMER_RIGHT_0_stt = self._get_abspos(
            PRIMER_RIGHT_0_stt, prinfo.abs_frag_pad_pre_stt)
        abs_PRIMER_RIGHT_0_end = \
            abs_PRIMER_RIGHT_0_stt + PRIMER_RIGHT_0_len - 1

        return \
            abs_PRIMER_LEFT_0_stt, abs_PRIMER_LEFT_0_end, \
            abs_PRIMER_RIGHT_0_stt, abs_PRIMER_RIGHT_0_end


    def _get_abspos(self, ref_pos, abs_template_stt):

        # 11
        #  1234567890
        #        7 11+7=18 -1

        return abs_template_stt + ref_pos - 1


class PrimerInfo(object):

    def __init__(self):

        self.marker_id = ''

        self.chrom = ''
        self.pos = 0

        self.targ_grp = ''
        self.g0_name = ''
        self.g1_name = ''

        self.gts_segr_lens = ''

        self.targ_ano = ''
        self.g0_ano = -1
        self.g1_ano = -1


        self.vseq_gno_str = ""

        self.set_enz_cnt = ''
        self.var_type = ''
        self.marker_info = ''

        self.vseq_lens_ano_str = ''

        self.target_gno = -1
        self.target_len = 0
        self.enzyme_name = ''
        self.digest_pattern = ''

        self.g0_seq_target_len = 0
        self.g0_seq_target = ''
        self.g1_seq_target_len = 0
        self.g1_seq_target = ''
        self.seq_template_ref_len = 0
        self.seq_template_ref_abs_pos = ''
        self.seq_template_ref_rel_pos = ''
        self.SEQUENCE_TARGET = ''
        self.seq_template_ref = ''

        # abs template info
        self.abs_frag_pad_pre_stt = 0
        self.abs_frag_pad_pre_end = 0
        self.abs_around_seq_pre_stt = 0
        self.abs_around_seq_pre_end = 0
        self.abs_pos = 0
        self.abs_around_seq_aft_stt = 0
        self.abs_around_seq_aft_end = 0
        self.abs_frag_pad_aft_stt = 0
        self.abs_frag_pad_aft_end = 0

        # abs template info
        self.rel_frag_pad_pre_stt = 0
        self.rel_frag_pad_pre_end = 0
        self.rel_around_seq_pre_stt = 0
        self.rel_around_seq_pre_end = 0
        self.rel_pos = 0
        self.rel_around_seq_aft_stt = 0
        self.rel_around_seq_aft_end = 0
        self.rel_frag_pad_aft_stt = 0
        self.rel_frag_pad_aft_end = 0

        self.SEQUENCE_EXCLUDED_REGION = ''
   
        self.p3_comment = ''


    def prepare_from_marker_file(self, distin_dict, marker_df_row):

        hdr_dict = distin_dict['marker']['hdr_dict']

        # basic
        self.marker_id, \
        self.chrom, \
        self.pos, \
        self.targ_grp, \
        self.g0_name, \
        self.g1_name, \
        self.targ_ano, \
        self.g0_ano, \
        self.g1_ano, \
        self.vseq_gno_str, \
        self.gts_segr_lens, \
        self.var_type, \
        self.set_enz_cnt, \
        self.marker_info, \
        self.vseq_lens_ano_str, \
        self.enzyme_name, \
        self.digest_pattern, \
        self.target_gno, \
        self.target_len = \
            utl.get_basic_primer_info(marker_df_row, hdr_dict)

        self.g0_seq_target_len = \
            int(marker_df_row[hdr_dict['g0_seq_target_len']])
        self.g0_seq_target = \
            str(marker_df_row[hdr_dict['g0_seq_target']])
        self.g1_seq_target_len = \
            int(marker_df_row[hdr_dict['g1_seq_target_len']])
        self.g1_seq_target = \
            str(marker_df_row[hdr_dict['g1_seq_target']])

        self.seq_template_ref_len = \
            int(marker_df_row[hdr_dict['seq_template_ref_len']])
        self.seq_template_ref_abs_pos = \
            str(marker_df_row[hdr_dict['seq_template_ref_abs_pos']])

        #log.debug("{}".format(self.seq_template_ref_abs_pos))

        self.seq_template_ref_rel_pos = \
            str(marker_df_row[hdr_dict['seq_template_ref_rel_pos']])

        #log.debug("{}".format(self.seq_template_ref_rel_pos))

        self.SEQUENCE_TARGET = \
            str(marker_df_row[hdr_dict['SEQUENCE_TARGET']])
        self.seq_template_ref = \
            str(marker_df_row[hdr_dict['seq_template_ref']])

        # abs template info
        self.abs_frag_pad_pre_stt, \
        self.abs_frag_pad_pre_end, \
        self.abs_around_seq_pre_stt, \
        self.abs_around_seq_pre_end, \
        self.abs_pos, \
        self.abs_around_seq_aft_stt, \
        self.abs_around_seq_aft_end, \
        self.abs_frag_pad_aft_stt, \
        self.abs_frag_pad_aft_end = \
            Product.separate_seq_template_pos(
                self.seq_template_ref_abs_pos)

        # rel template info
        self.rel_frag_pad_pre_stt, \
        self.rel_frag_pad_pre_end, \
        self.rel_around_seq_pre_stt, \
        self.rel_around_seq_pre_end, \
        self.rel_pos, \
        self.rel_around_seq_aft_stt, \
        self.rel_around_seq_aft_end, \
        self.rel_frag_pad_aft_stt, \
        self.rel_frag_pad_aft_end = \
            Product.separate_seq_template_pos(
                self.seq_template_ref_rel_pos)


    def _get_relpos(self, abs_pos):

        #   | self.abs_frag_pad_pre_stt
        #                 self.abs_frag_pad_aft_end|
        #   <-------><========>P<=========><------->
        #            |self.abs_around_seq_pre_stt
        #      self.abs_around_seq_aft_end|
        #   |101 109|
        #   123456789
        #       |105
        #   105-101+1 = 5

        rel_pos = abs_pos - self.abs_frag_pad_pre_stt + 1

        return rel_pos


    def get_excluded_region(self):

        SEQUENCE_EXCLUDED_REGION = list()

        #logf_l = ["{} self.pos={} rel_pos={} rel_end_pos={} "]
        #logf_l += ["region_len={} template_len={}"]
        #logf = "".join(logf_l)

        region = "{}:{}-{}".format(
            self.chrom,
            self.abs_frag_pad_pre_stt,
            self.abs_frag_pad_aft_end)

        reader = vcfpy.Reader.from_path(glv.conf.vcf_file_path)
        vcf_ittr = reader.fetch(region)

        # access to vcf using iterater
        for record in vcf_ittr:
            sample0 = glv.conf.group_members_dict[self.g0_name][0]
            sample1 = glv.conf.group_members_dict[self.g1_name][0]

            # もし、サンプル間でvariantが見つかった場合は、
            s0_0, s0_1, s1_0, s1_1 = \
                AlleleSelect.record_call_for_sample(record, sample0, sample1)

            if utl.is_same_gt(s0_0, s0_1, s1_0, s1_1) == False:

                # 20200713 here
                if self.pos != record.POS:
                    rel_pos = self._get_relpos(record.POS) #
                    # そのポジションのrefのvseq分を登録する
                    # 長さがfragment長を超える場合は、
                    # 調整する。
                    # 見つかったのはPOS
                    # REFのlength

                    # self.pos|
                    #     1036|
                    # ATGCATGCA ref_len=1
                    #         T
                    #         C
                    #         1036 + 1 - 1
                    region_len = len(record.REF)
                    rel_end_pos = rel_pos + region_len - 1

                    #  pos   len     end
                    #  1036 (10)     1045
                    #            1041 temp_len

#                    log.debug(logf.format(
#                        1, self.pos, rel_pos, rel_end_pos, region_len,
#                        self.seq_template_ref_len))


                    #log.debug("{}, {}".format(
                    #    rel_end_pos, self.seq_template_ref_len))
                    #log.debug("{}, {}".format(
                    #    type(rel_end_pos), type(self.seq_template_ref_len)))


                    if rel_end_pos > self.seq_template_ref_len:
                        diff_len = rel_end_pos - self.seq_template_ref_len
                        region_len = region_len - diff_len

#                        log.debug(logf.format(
#                            2, self.pos, rel_pos, rel_end_pos, region_len,
#                            self.seq_template_ref_len))

                    SEQUENCE_EXCLUDED_REGION += [
                        "{},{}".format(rel_pos, region_len)]

        self.SEQUENCE_EXCLUDED_REGION = " ".join(SEQUENCE_EXCLUDED_REGION)
        #log.debug("SEQUENCE_EXCLUDED_REGION={}".format(
        #    self.SEQUENCE_EXCLUDED_REGION))


    def make_p3_info(self):

        # save information
        self.p3_comment = self._make_p3_comment()
        #log.debug("{}".format(self.p3_comment))

        self.iopr3 = InoutPrimer3()
        self.iopr3.set_p3_comment(self.p3_comment)
        self.iopr3.set_sequence_target(self.SEQUENCE_TARGET)
        self.iopr3.add_ex_region(self.SEQUENCE_EXCLUDED_REGION)
        self.iopr3.set_sequence_id(self.marker_id)
        self.iopr3.set_sequence_template(self.seq_template_ref)

        #log.debug("")
        #log.debug("{}".format(self.iopr3.get_p3_input()))


    def _make_p3_comment(self):

        p3_comment = list()
        p3_comment += ["{}:{}".format('marker_id', self.marker_id)]
        p3_comment += ["{}:{}".format('var_type', self.var_type)]
        p3_comment += ["{}:{}".format('g0_name', self.g0_name)]
        p3_comment += ["{}:{}".format('g1_name', self.g1_name)]
        p3_comment += ["{}:{}".format('marker_info', self.marker_info)]


        p3_comment += ["{}:{}".format(
            'g0_seq_target_len', self.g0_seq_target_len)]
        p3_comment += ["{}:{}".format(
            'g0_seq_target', self.g0_seq_target)]
        p3_comment += ["{}:{}".format(
            'g1_seq_target_len', self.g1_seq_target_len)]
        p3_comment += ["{}:{}".format(
            'g1_seq_target', self.g1_seq_target)]

        p3_comment += ["{}:{}".format(
            'seq_template_ref_len', self.seq_template_ref_len)]
        p3_comment += ["{}:{}".format(
            'seq_template_ref_abs_pos', self.seq_template_ref_abs_pos)]
        p3_comment += ["{}:{}".format(
            'seq_template_ref_rel_pos', self.seq_template_ref_rel_pos)]

        return ';'.join(map(str, p3_comment))




