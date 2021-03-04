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


class OutList(object):

    def __init__(self):

        self.outf_prefix = {
            'prepare'       : {'no':  0,  'fn': ''},
            'variant'       : {'no': 10,  'fn': '010_variant'},
            'marker'        : {'no': 20,  'fn': '020_marker'},
            'primer'        : {'no': 30,  'fn': '030_primer'},
            'formfail'      : {'no': 40,  'fn': '040_formatF'},
            'formsafe'      : {'no': 50,  'fn': '050_formatS'},
        }

        # list of full path of out text
        self.distin_files = list()

        # stopと、progressと。
        # ある地点の、現在地番号。
        #   現在地番号が、   <progress(ここから実行)  飛ばす
        #                   ==progress(ここから実行)  実行する
        #                    >progress(ここから実行)  実行する
        #   現在地番号が、   <stop(ここが終わったら止まる)  実行する 
        #                   ==stop(ここが終わったら止まる)  実行する
        #                    >stop(ここが終わったら止まる)  止まる
        # 
        # のprogressのstatusと stopのstatusを返して、
        #               prepare variant marker  primer  formfail  formsafe
        # progress(all) x         |->     |->     |->     |->       |-> 
        # stop(no)      x       ->|     ->|     ->|     ->|       ->|

    def open_log(self):

        global log
        log = LogConf.open_log(__name__)


    def prepare_distin_files(self):
        ''' access distinct_group information with outfile info
        '''

        for distin_dict in glv.conf.distinguish_groups_list:
            # self.distinguish_groups_list
            # [ {
            #   0: 'ref',
            #   1: 'g0',
            #   'indel_size': '50-200',
            #   'pick_mode': 'indel',
            #   'regions': ['chrom_01', 'chrom_17']
            #   },

            # self.distin_files=
            # [ {
            #   0: 'ref',
            #   1: 'g0',
            #   'indel_size': '50-200',
            #   'pick_mode': 'indel',
            #   'region': 'chrom_01',
            #   'regions_list': ['chrom_01', 'chrom_17'],
            #   .....
            #   },

            # For each region in the list
            for reg in distin_dict['regions']:

                # reconstruct cause region is list
                distin_out_dict = {
                    0: distin_dict[0],
                    1: distin_dict[1],
                    'region': reg,
                    'disting_str': distin_dict['disting_str'],
                    'regions_list': distin_dict['regions'],
                    'pick_mode': distin_dict['pick_mode'],
                    'indel_size': distin_dict['indel_size'],
                }
                #log.debug("{}".format(distin_out_dict))

                out_file_dict = dict()

                for key, no_fn in self.outf_prefix.items():

                    # file name to write out result to text
                    out_file_path, basename = self._make_distin_fname(
                        no_fn['fn'],
                        distin_dict[0],             # g0
                        distin_dict[1],             # g1
                        reg,                        # region
                        distin_dict['pick_mode'])   # pick_mode

                    hdr_text, hdr_list, hdr_dict = self.make_header(key)

                    out_file_dict[key] = {
                        'out_path': out_file_path,
                        'base_nam': basename,
                        'hdr_text': hdr_text,
                        'hdr_list': hdr_list,
                        'hdr_dict': hdr_dict,
                    }

                # add dictionary item
                distin_out_dict.update(out_file_dict)
                #print("{}\n\n".format(distin_out_dict))

                # make list of dictionary
                self.distin_files.append(distin_out_dict)

        #log.info("\nself.distin_files=\n{}\n".format(
        #    pprint.pformat(self.distin_files)))


    def _make_distin_fname(
        self, outf_pref, distin_0, distin_1, rg, pick_mode):
        """
        """

        #distin~gHitomebore~gKaluheenati~rg0~all~50-200.txt
        basename = "{}~{}~{}~{}~{}~i{}-{}~p{}-{}".format(
                outf_pref,
                distin_0,
                distin_1,
                rg,
                pick_mode,
                glv.conf.min_indel_len,
                glv.conf.max_indel_len,
                glv.conf.min_product_size,
                glv.conf.max_product_size)

        out_file_path = "{}/{}.txt".format(glv.conf.out_dir_path, basename)

        return out_file_path, basename


    def make_header(self, type):

        hd_l = list()
        hd_d = dict()
        i = 1   # 0 is index

        if type == 'variant':
            # Synchronize with allele_select.py _construct_var_line
            hd_l, hd_d, i = self._mkmap('chrom',         hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('pos',           hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('targ_grp',      hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('targ_ano',      hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('vseq_gno_str',  hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('gts_segr_lens', hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('var_type',      hd_l, hd_d, i)
            # ----------------------
            hd_l, hd_d, i = self._mkmap('set_n',         hd_l, hd_d, i)
            # ----------------------
            hd_l, hd_d, i = self._mkmap('len_g0g1_dif_long', hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('vseq_ano_str',  hd_l, hd_d, i)
            #Allele int information for each sample will be added

        elif type == 'marker':
            # Sync with eval_variant.py
            # copy_line_for_effective_restriction_enzymes
            hd_l, hd_d, i = self._mkmap('marker_id',   hd_l, hd_d, i)
            # --------------
            hd_l, hd_d, i = self._mkmap('chrom',         hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('pos',           hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('targ_grp',      hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('targ_ano',      hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('vseq_gno_str',  hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('gts_segr_lens', hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('var_type',      hd_l, hd_d, i)
            # ----------------------
            hd_l, hd_d, i = self._mkmap('set_enz_cnt',   hd_l, hd_d, i)
            # ----------------------
            hd_l, hd_d, i = self._mkmap('marker_info',   hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('vseq_lens_ano_str', hd_l, hd_d, i)
            # --------------

            # g0
            hd_l, hd_d, i = self._mkmap('g0_seq_target_len',
                 hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('g0_seq_target',
                 hd_l, hd_d, i)
            # g1
            hd_l, hd_d, i = self._mkmap('g1_seq_target_len',
                 hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('g1_seq_target',
                 hd_l, hd_d, i)

            # seq_template_ref
            hd_l, hd_d, i = self._mkmap('seq_template_ref_len',
                 hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('seq_template_ref_abs_pos',
                 hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('seq_template_ref_rel_pos',
                 hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('SEQUENCE_TARGET',
                 hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('seq_template_ref',
                 hd_l, hd_d, i)

        elif type == 'primer':  # there is no header

            # Synchronize with primer.py _primer_complete_to_line
            hd_l, hd_d, i = self._mkmap('marker_id',   hd_l, hd_d, i)
            # --------------
            hd_l, hd_d, i = self._mkmap('chrom',       hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('pos',         hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('targ_grp',    hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('targ_ano',  hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('vseq_gno_str',  hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('gts_segr_lens', hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('var_type',    hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('set_enz_cnt', hd_l, hd_d, i)
            # --------------
            hd_l, hd_d, i = self._mkmap('marker_info', hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('vseq_lens_ano_str', hd_l, hd_d, i)
            # --------------
            hd_l, hd_d, i = self._mkmap('g0_seq_target_len',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('g0_seq_target',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('g1_seq_target_len',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('g1_seq_target',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('seq_template_ref_len',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('seq_template_ref_abs_pos',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('seq_template_ref_rel_pos',
                hd_l, hd_d, i)
            # --------------

            hd_l, hd_d, i = self._mkmap('try_cnt',     hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('complete',    hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('blast_check', hd_l, hd_d, i)
            # --------------

            hd_l, hd_d, i = self._mkmap('PRIMER_PAIR_0_PRODUCT_SIZE',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('PRIMER_LEFT_0',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('left_primer_id',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('PRIMER_LEFT_0_SEQUENCE',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('PRIMER_RIGHT_0',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('right_primer_id',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('PRIMER_RIGHT_0_SEQUENCE',
                hd_l, hd_d, i)

            hd_l, hd_d, i = self._mkmap('SEQUENCE_TARGET',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('SEQUENCE_EXCLUDED_REGION',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('SEQUENCE_TEMPLATE',
                hd_l, hd_d, i)

        elif type == 'formsafe':

            # Synchronize with formtxt.py _format_product

            hd_l, hd_d, i = self._mkmap('chrom',            hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('pos',              hd_l, hd_d, i)

            hd_l, hd_d, i = self._mkmap('g0_vseq',          hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('g1_vseq',          hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('g0_gt',          hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('g1_gt',          hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('targ_ano',      hd_l, hd_d, i)

            hd_l, hd_d, i = self._mkmap('var_type',         hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('comment',          hd_l, hd_d, i)

            hd_l, hd_d, i = self._mkmap('enzyme',           hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('g0_name',          hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('g1_name',          hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('g0_product_size',  hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('g1_product_size',  hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('diff_length',      hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('g0_digested_size', hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('g1_digested_size', hd_l, hd_d, i)

            hd_l, hd_d, i = self._mkmap('digested_gno',     hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('digested_ano',     hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('try_cnt',          hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('complete',         hd_l, hd_d, i)

            hd_l, hd_d, i = self._mkmap('marker_id',            hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('gts_segr_lens',    hd_l, hd_d, i)

            hd_l, hd_d, i = self._mkmap('left_primer_id',   hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('PRIMER_LEFT_0_SEQUENCE',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('right_primer_id',  hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('PRIMER_RIGHT_0_SEQUENCE',
                hd_l, hd_d, i)


        hd_str = '\t'.join(map(str, hd_l))
        return hd_str, hd_l, hd_d


    def _mkmap(self, key, header_list, header_dict, idx):

        #log.debug("{} {} {} {}".format(
        #    key, header, header_dict, idx,
        #    ))

        header_list += [key]
        add_dict = {key: idx}
        header_dict.update(add_dict)
        idx += 1

        return header_list, header_dict, idx


    def complete_list(self):
        pass

