# -*- coding: utf-8 -*-

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

import pandas as pd

from vprimer.eval_variant import EvalVariant
from vprimer.product import Product
from vprimer.primer import PrimerInfo

class FormTxt(object):

    def __init__(self):

        self.line = ''

        self.marker_id = ''

        self.chrom = ''
        self.pos = 0

        self.vseq_gno_str = ''
        self.g0_vseq = ''
        self.g1_vseq = ''

        self.g0_gt = ''
        self.g1_gt = ''

        self.targ_grp = ''
        self.g0_name = ''
        self.g1_name = ''

        self.g0_vseq = ''
        self.g1_vseq = ''

        self.gts_segr_lens = ''
        self.g0_genotype = ''
        self.g1_genotype = ''

        self.targ_ano = ''
        self.g0_ano = -1
        self.g1_ano = -1

        self.set_enz_cnt = ''
        self.var_type = ''
        self.marker_info = ''

        self.vseq_lens_ano_str = ''

        self.target_gno = -1
        self.target_len = 0
        self.enzyme_name = ''
        self.digest_pattern = ''

        self.g0_seq_target_len = ''
        self.g0_seq_target = ''
        self.g1_seq_target_len = ''
        self.g1_seq_target = ''
        self.seq_template_ref_len = ''
        self.seq_template_ref_abs_pos = ''
        self.seq_template_ref_rel_pos = ''
        self.PRIMER_PAIR_0_PRODUCT_SIZE = ''
        self.PRIMER_LEFT_0 = ''
        self.left_primer_id = ''
        self.PRIMER_LEFT_0_SEQUENCE = ''
        self.PRIMER_RIGHT_0 = ''
        self.right_primer_id = ''
        self.PRIMER_RIGHT_0_SEQUENCE = ''
        self.SEQUENCE_TEMPLATE = ''


    def format_text(self):
        '''
        '''

        # for each distinguish_groups
        for proc_cnt, distin_dict in enumerate(glv.outlist.distin_files, 1):

            # 
            # read primer file
            primer_file = distin_dict['primer']['out_path']

            # read variant file and set allele int informations
            # to a dictionary.
            variant_file = distin_dict['variant']['out_path']

            df_variant = pd.read_csv(
                variant_file, sep='\t', header=0, index_col=None)
            header_list = distin_dict['variant']['hdr_text'].split("\t")
            existing_column_cnt = len(header_list)
            # not including REF,ALT
            #alint_start = existing_column_cnt + 1 - 1
            alint_start = existing_column_cnt + 1

            variant_alint_dict = dict()
            alint_list = list()

            for variant_df_row in df_variant.itertuples():
                chrom_name = variant_df_row[1]
                pos = variant_df_row[2]
                alint_list = [variant_df_row[6]]
                alint_list += list(variant_df_row[alint_start:])

                if chrom_name not in variant_alint_dict.keys():
                    variant_alint_dict[chrom_name] = dict()

                variant_alint_dict[chrom_name][pos] = alint_list

            #--------------------------------------------------------

            df_distin = pd.read_csv(
                primer_file, sep='\t', header=0, index_col=None)

            # complete == 1 or == 0
            fail = 0
            safe = 1
            for complete, proc in zip([fail, safe], ['formfail', 'formsafe']):

                # stop, action, gothrough
                proc_name = proc
                ret_status = utl.decide_action_stop(proc_name)

                if ret_status == "stop":
                    msg = "STOP. "
                    msg += "Current process \'{}\' ".format(proc_name)
                    msg += "has exceeded the User-specified stop point "
                    msg += "\'{}', ".format(glv.conf.stop)
                    msg += "so stop program. exit."
                    log.info(msg)
                    #sys.exit(1)
                    continue


                elif ret_status == "gothrough":
                    msg = "SKIP \'{}\' proc, ".format(proc_name)
                    msg += "glv.conf.progress = {}, ".format(
                        glv.conf.progress)
                    msg += "glv.conf.stop = {}, ".format(glv.conf.stop)
                    msg += "so skip program."
                    log.info(msg)
                    continue


                log.info("-------------------------------")
                log.info("Start processing {} complete={}\n".format(
                    proc_name, complete))

                # logging current target
                sub_proc = "{}_{}".format(proc, complete)
                utl.print_distin_info(sub_proc, distin_dict, proc_cnt, True)

                df_distin_complete = \
                    df_distin[df_distin['complete'] == complete]

                #------------------------
                # check chrom-pos duplicate marker 
                df_chrom_pos = df_distin_complete.loc[:, ['chrom', 'pos']]
                df_chrom_pos_duplicated = \
                    df_chrom_pos[df_chrom_pos.duplicated()]

                duplicate_pos_dict = dict()
                for c_p_row in df_chrom_pos_duplicated.itertuples():

                    chrom = c_p_row[1]
                    pos = c_p_row[2]

                    if not chrom in duplicate_pos_dict:
                        duplicate_pos_dict[chrom] = dict();

                    if not pos in duplicate_pos_dict[chrom]:
                        duplicate_pos_dict[chrom][pos] = pos

                #------------------------
                # file name to write out result to text
                out_txt_file = distin_dict[proc]['out_path']
                log.info("out_txt_file={}.".format(out_txt_file))

                utl.save_to_tmpfile(out_txt_file)


                with open(out_txt_file, mode='a') as f:

                    header = distin_dict['formsafe']['hdr_text']
                    if (proc == "formsafe"):
                        #alint_header = ["targ_ano", "vseq_ano_str"]
                        alint_header = ["vseq_ano_str"]
                        sample_nickname_ordered_list, \
                        sample_fullname_ordered_list = \
                            utl.get_ordered_sample_list(
                                [distin_dict[0], distin_dict[1]])
                        alint_header += sample_nickname_ordered_list
                        header = "{}\t{}".format(
                            header, "\t".join(alint_header))

                    # write header
                    f.write("{}\n".format(header))

                    # each variant
                    for primer_df_row in df_distin_complete.itertuples():

                        chrom_name = primer_df_row[2]
                        pos = primer_df_row[3]

                        self._prepare_from_primer_file(
                            primer_df_row, distin_dict)

                        self._format_product(duplicate_pos_dict)

                        if (proc == "formsafe"):

                            #print("chrom_name={}, pos={}".format(
                            #    chrom_name, pos))
                            #print("{}, {}".format(chrom_name, pos))
                            #slice_one = variant_alint_dict[chrom_name][pos]
                            #print(type(slice_one))

                            #pprint.pprint(
                            #   variant_alint_dict[chrom_name][pos])

                            line = "{}\t{}".format(
                                self.line, "\t".join(
                                    map(str,
                                        variant_alint_dict[chrom_name][pos])))

                            f.write("{}\n".format(line))

                        else:
                            # 書き出す
                            f.write("{}\n".format(self.line))


    def _get_group_product_size(self):
        ''' The product size obtained by the output of primer3 is the product
        size on the REF sequence, and the products of each group need to be
        adjusted by the difference from the variant length of REF.
        indel has no digest. CAPS can get the size information to be digested.
        Obtain the digest size by CAPS from the product size of the group.

            # self.target_gno 1,0
            # rel_frag_pad_pre_end
            #                   ||rel_around_seq_pre_stt
            #                   ||        |rel_pos
            #                   ||        --- vseq_len_ref
            # <-----------------><========<=>========><----------------->
            #                500|
            #       |313
                    <--len 188--> 
            #       |313             len=421                  733|
            #       <==>......................................<==>
            #       |product_stt_pos              product_end_pos|
            #                    188+21+1=210
            #                    <========<=>========>
            #             found_pos 21|>d_pos 1
            #       |313             len=421                  733|
            #       <g0_prod_stt   vseq_g0<=>         g0_prod_end>
            #       <g1_prod_stt   vseq_g1<=>         g1_prod_end>
            #       |313             len=421                  733|
            #       <=== L_digested ===><====== R_digested ======>
            #
            #  <-- len 188 -->1234567890123456789012
            #                                     |found
            #                                      >digest
            #                 9012345678901234567890 => 210 digest pos


        '''

#        print("self.seq_template_ref_rel_pos={}".format(
#            self.seq_template_ref_rel_pos))

        # Get information on the relative position of seq_template
        rel_frag_pad_pre_stt, \
        rel_frag_pad_pre_end, \
        rel_around_seq_pre_stt, \
        rel_around_seq_pre_end, \
        rel_pos, \
        rel_around_seq_aft_stt, \
        rel_around_seq_aft_end, \
        rel_frag_pad_aft_stt, \
        rel_frag_pad_aft_end = \
            Product.separate_seq_template_pos(
                self.seq_template_ref_rel_pos)

        #----------------------------------------------------
        # vseq information etc.
        #----------------------------------------------------
        # Correspondence table of ano to gno
        ano_gno = [self.g0_ano, self.g1_ano]
        #print("ano_gno={}".format(ano_gno))

        # List the lengths of variants in alelle number order.
        vseq_lens_ano = [int(x) for x in self.vseq_lens_ano_str.split(',')]
        #print("vseq_lens_ano={}".format(vseq_lens_ano))

        # Variant length of REF
        vseq_len_ref = vseq_lens_ano[0]
#        print("vseq_len_ref={}".format(vseq_len_ref))

        #----------------------------------------------------
        # positions and lengths of left and right primers reported by primer3
        #----------------------------------------------------
        # The position of PRIMER_LEFT is the starting point of the product.
        ref_product_stt, left_len = map(int, self.PRIMER_LEFT_0.split(','))
#        print("ref_product_stt={}, left_len={}".format(
#            ref_product_stt, left_len))

        # The position of PRIMER_RIGHT is the end point of the product
        # because it is the start point of the reverse complement sequence.
        ref_product_end, right_len = map(int, self.PRIMER_RIGHT_0.split(','))
#        print("ref_product_end={}, right_len={}".format(
#            ref_product_end, right_len))

        #----------------------------------------------------
        # Calculation of product size considering the difference between the
        # variant length of ref and that of group 0 and group 1.
        l_product_size = list()
        l_product_end_pos = list()
        len_padding_left = rel_frag_pad_pre_end - ref_product_stt + 1
#        print("len_padding_left={}, {} - {} + 1".format(
#            len_padding_left, rel_frag_pad_pre_end, ref_product_stt))

        # For each group
        for gno in range(2):
            # Calculate the length difference from REF for vseq

#            print("gno={}".format(gno))
            vseq_len_gno = vseq_lens_ano[ano_gno[gno]]
#            print("vseq_len_gno={}".format(vseq_len_gno))
            diff_vseq_len = vseq_len_ref - vseq_len_gno
#            print("diff_vseq_len={}".format(diff_vseq_len))

            # Adjusting the product size and end point for each group
            my_product_end = ref_product_end - diff_vseq_len
#            print("my_product_end={}".format(my_product_end))

            my_product_size = my_product_end - ref_product_stt + 1
#            print("my_product_size={}".format(my_product_size))

            l_product_end_pos.append(my_product_end)
            l_product_size.append(my_product_size)

        # Difference in product size between groups
        diff_product_size = abs(l_product_size[0] - l_product_size[1])
#        print("diff_product_size={}".format(diff_product_size))


#        # ----------------------------------------------------------------
#        # Digest size calculation

        # for indel
        # longer_group,
        # longer_length,
        # shorter_length,
        # diff_length,
        # digested_pos

        # If !indel, For those who are digested, use'/' to enter two numbers
        # (large or small). If not, enter the product size directly.
        digested_size = [[0 for j in range(0)] for i in range(2)]

        # Group number to be digested or longer indel_len
        digested_gno = self.target_gno
        not_digested_gno = utl.flip_gno(digested_gno)

        digested_ano = ano_gno[digested_gno]
        not_digested_ano = ano_gno[not_digested_gno]

#        print("digested_gno={}".format(digested_gno))
#        print("non_digested_gno={}".format(non_digested_gno))

        # if glv.INDEL, this is indel length diff. !=INDEL, it is digest pos
        digest_diff = diff_product_size

        if self.var_type == glv.INDEL:
            # If indel, put two product sizes directly
            digested_size_str_g0 = str(l_product_size[0])
            digested_size_str_g1 = str(l_product_size[1])
            # indel is completed
            digested_size_str_g0 = "-"
            digested_size_str_g1 = "-"

        else:
            # for ! INDEL
            enzyme_name, \
            d_gno, \
            found_pos, \
            digest_pattern, \
            digested_pos \
                = EvalVariant.split_marker_info(self.marker_info)

            found_pos = int(found_pos)
            digested_pos = int(digested_pos)

#            print()
#            print("self.marker_info={}".format(self.marker_info))
#            print("enzyme_name={}".format(enzyme_name))
#            print("digested_gno={}".format(digested_gno))
#            print("digested_ano={}".format(digested_ano))
#            print("found_pos={}".format(found_pos))
#            print("digest_pattern={}".format(digest_pattern))
#            print("digested_pos={}".format(digested_pos))

            dig_pos = len_padding_left + found_pos + digested_pos
            #print("dig_pos={}, {}+{}+{}".format(
            #    dig_pos, len_padding_left, found_pos, digested_pos))

            L_digested_len = dig_pos
            R_digested_len = l_product_size[digested_gno] - L_digested_len

            #                    |ref_pos
            # <----------><======<.....>=====><----------->
            # <----------><======|=====><----------->

            # dig_gno=0
            # <----------><======<.....>=====><----------->
            #                       <-^-->
            # <----------><======|=====><----------->

            # dig_gno=1
            # <----------><======<.....>=====><----------->
            # <----------><======|=====><----------->
            #                   <-^-->

            #                    |ref_pos
            # <----------><======|=====><----------->
            # <----------><======<.....>=====><----------->

            # dig_gno=0
            #                    |ref_pos
            # <----------><======|=====><----------->
            #                   <-^-->
            # <----------><======<.....>=====><----------->

            # dig_gno=1
            #                    |ref_pos
            # <----------><======|=====><----------->
            # <----------><======<.....>=====><----------->
            #               <-^-->


            if str(digested_ano) in self.g0_gt:
                digested_size[0].append(L_digested_len)
                digested_size[0].append(R_digested_len)

            if str(not_digested_ano) in self.g0_gt:
                digested_size[0].append(
                    l_product_size[not_digested_gno])
                
            if str(digested_ano) in self.g1_gt:
                digested_size[1].append(L_digested_len)
                digested_size[1].append(R_digested_len)

            if str(not_digested_ano) in self.g1_gt:
                digested_size[1].append(
                    l_product_size[not_digested_gno])
                
#            print("targ_ano={}".format(self.targ_ano))
#            print("digested_gno={}, digested_ano={}".format(
#                digested_gno, digested_ano))
#            print("g0_gt={}, {}".format(self.g0_gt, digested_size[0]))
#            print("g1_gt={}, {}".format(self.g1_gt, digested_size[1]))

#            print("digested_ano={}".format(digested_ano))
#            print("not_digested_ano={}".format(not_digested_ano))

            digested_size_str_g0 = \
                "/".join(
                    map(str, sorted(set(digested_size[0]), reverse=True)))
            digested_size_str_g1 = \
                "/".join(
                    map(str, sorted(set(digested_size[1]), reverse=True)))
#
#            print("digested_size_str_g0={}".format(digested_size_str_g0))
#            print("digested_size_str_g1={}".format(digested_size_str_g1))

        return \
            l_product_size[0], \
            l_product_size[1], \
            digested_gno, \
            digested_ano, \
            digested_size_str_g0, \
            digested_size_str_g1, \
            diff_product_size, \
            digest_diff


    def _add_comment(self, duplicate_pos_dict):

        comment_list = list()

        #log.debug("{}".format(duplicate_pos_dict))

        if len(duplicate_pos_dict) == 0:
            pass
        elif self.pos in duplicate_pos_dict[self.chrom]:
            comment_list = [glv.COMMENT_dup]

        if len(comment_list) != 0:
            comment_list += [self.set_enz_cnt]

        if len(comment_list) == 0:
            comment_list = '-'

        return ','.join(comment_list)


    def _format_product(self, duplicate_pos_dict):

        product_size_0, \
        product_size_1, \
        digested_gno, \
        digested_ano, \
        digested_size_0, \
        digested_size_1, \
        diff_product_size, \
        digest_diff = \
            self._get_group_product_size()

        comment = self._add_comment(duplicate_pos_dict)

        #--------------------
        line_list = list()


        line_list += [self.chrom]
        line_list += [self.pos]

        line_list += [self.g0_vseq]
        line_list += [self.g1_vseq]
        line_list += [self.g0_gt]
        line_list += [self.g1_gt]
        line_list += [self.targ_ano]

        line_list += [self.var_type]
        line_list += [comment]

        line_list += [self.enzyme_name]
        line_list += [self.g0_name]
        line_list += [self.g1_name]
        line_list += [product_size_0]
        line_list += [product_size_1]

        # **
        #line_list += [digest_diff]
        line_list += [diff_product_size]
        line_list += [digested_size_0]
        line_list += [digested_size_1]

        line_list += [digested_gno]
        line_list += [digested_ano]
        line_list += [self.try_cnt]
        line_list += [self.complete]

        line_list += [self.marker_id]
        line_list += [self.gts_segr_lens]

        line_list += [self.left_primer_id]
        line_list += [self.PRIMER_LEFT_0_SEQUENCE]

        line_list += [self.right_primer_id]
        line_list += [self.PRIMER_RIGHT_0_SEQUENCE]

        self.line = '\t'.join(map(str, line_list))


    def _prepare_from_primer_file(self, primer_df_row, distin_dict):

        hdr_dict = distin_dict['primer']['hdr_dict']

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
            utl.get_basic_primer_info(primer_df_row, hdr_dict)

        self.g0_vseq, self.g1_vseq = self.vseq_gno_str.split(",")
        #11/00,hoho_1,1.1/1.1
        self.g0_gt, self.g1_gt = self._get_genotype()

        #log.debug("self.chrom={} pos={}".format(self.chrom, self.pos))

        self.try_cnt = str(primer_df_row[hdr_dict['try_cnt']])
        self.complete = str(primer_df_row[hdr_dict['complete']])
        self.blast_check = str(primer_df_row[hdr_dict['blast_check']])

        self.g0_seq_target_len = \
            int(primer_df_row[hdr_dict['g0_seq_target_len']])
        self.g0_seq_target = \
            str(primer_df_row[hdr_dict['g0_seq_target']])
        self.g1_seq_target_len = \
            int(primer_df_row[hdr_dict['g1_seq_target_len']])
        self.g1_seq_target = \
            str(primer_df_row[hdr_dict['g1_seq_target']])

        self.seq_template_ref_len = \
            int(primer_df_row[hdr_dict['seq_template_ref_len']])
        self.seq_template_ref_abs_pos = \
            str(primer_df_row[hdr_dict['seq_template_ref_abs_pos']])
        self.seq_template_ref_rel_pos = \
            str(primer_df_row[hdr_dict['seq_template_ref_rel_pos']])

        self.PRIMER_PAIR_0_PRODUCT_SIZE = \
            int(primer_df_row[hdr_dict['PRIMER_PAIR_0_PRODUCT_SIZE']])
        self.PRIMER_LEFT_0 = \
            str(primer_df_row[hdr_dict['PRIMER_LEFT_0']])
        self.left_primer_id = \
            str(primer_df_row[hdr_dict['left_primer_id']])
        self.PRIMER_LEFT_0_SEQUENCE = \
            str(primer_df_row[hdr_dict['PRIMER_LEFT_0_SEQUENCE']])
        self.PRIMER_RIGHT_0 = \
            str(primer_df_row[hdr_dict['PRIMER_RIGHT_0']])
        self.right_primer_id = \
            str(primer_df_row[hdr_dict['right_primer_id']])
        self.PRIMER_RIGHT_0_SEQUENCE = \
            str(primer_df_row[hdr_dict['PRIMER_RIGHT_0_SEQUENCE']])
        self.SEQUENCE_TEMPLATE = \
            str(primer_df_row[hdr_dict['SEQUENCE_TEMPLATE']])

    def _get_genotype(self):
        '''
        '''

        # 11/00,hoho_1,1.1/1.1
        # self.gts_segr_lens
        #print("self.gts_segr_lens={}".format(self.gts_segr_lens))
        gts = self.gts_segr_lens.split(",")[0] 
        #print("gts={}".format(gts))
        g0_gt, g1_gt = gts.split("/")
        #print("g0_gt={}, g1_gt={}".format(g0_gt, g1_gt))
        g0_gtl = list(g0_gt)
        #print("g0_gtl={}".format(g0_gtl))
        g1_gtl = list(g1_gt)

        g0_genotype = "{}/{}".format(g0_gtl[0], g0_gtl[1])
        g1_genotype = "{}/{}".format(g1_gtl[0], g1_gtl[1])

        #print("{}, {}".format(g0_genotype, g1_genotype))

        return g0_genotype, g1_genotype





