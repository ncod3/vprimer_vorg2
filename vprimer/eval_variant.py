import sys
import os
import errno
import pprint
import re

import logging
log = logging.getLogger(__name__)

# global configuration
import vprimer.glv as glv
import vprimer.utils as utl

from Bio import Restriction
from Bio.Seq import Seq
# error
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

from vprimer.product import Product

class EvalVariant(object):

    def __init__(self, enzyme_name_list):

        self.enzyme_name_list = enzyme_name_list

        self.chrom = ''
        self.pos = 0

        self.targ_grp = ''
        self.g0_name = ''
        self.g1_name = ''
        self.gname_gno = list()

        self.gts_segr_lens = ''
        self.set_my_no = 0
        self.set_total_no = ''
        self.set_n = ""

        self.targ_ano = ''
        self.g0_ano = -1
        self.g1_ano = -1
        self.ano_gno = list()

        self.len_g0g1_dif_long = ''
        self.g0_len = 0
        self.g1_len = 0
        self.diff_len = 0
        self.long_gno = -1
        self.short_gno = -1
        self.longer_len = 0
        self.shorter_len = 0

        self.var_type = ''
        self.vseq_ano_str = ''
        self.vseq_ano = list()
        self.vseq_lens_ano = list()

        # --------------------------------
        self.marker_available = False

        self.vseq_ref = ''
        self.vseq_ref_len = ''

        # around_seq
        self.abs_around_seq_pre_stt = 0
        self.abs_around_seq_pre_end = 0
        self.abs_around_seq_aft_stt = 0
        self.abs_around_seq_aft_end = 0

        self.around_seq_pre = ''
        self.around_seq_aft = ''

        self.around_seq_pre_len = 0
        self.around_seq_aft_len = 0

        # target ref
        self.seq_target_ref = ''
        self.seq_target_ref_len = 0
        self.seq_target_ref_rel_pos = ''

        self.SEQUENCE_TARGET = ''

        # product
        self.g0_product = Product()
        self.g1_product = Product()
        self.gr_product = [self.g0_product, self.g1_product]

        # caps
        self.caps_result = dict()
        self.caps_found = 0

        # digest_gno:g0_prod_size:g1_prod_size:128/29/20
        #self.digest_info = ''

        # seq_template
        self.fragment_pad_len = 0

        self.abs_frag_pad_pre_stt = 0
        self.abs_frag_pad_pre_end = 0
        self.abs_frag_pad_aft_stt = 0
        self.abs_frag_pad_aft_end = 0

        self.seq_template_ref_abs_pos = ''
        self.seq_template_ref_rel_pos = ''

        self.frag_pad_pre = ''
        self.frag_pad_aft = ''
        self.seq_template_ref = ''

        self.frag_pad_pre_len = 0
        self.frag_pad_aft_len = 0
        self.seq_template_ref_len = 0


    def evaluate_for_marker(self, variant_df_row, distin_dict):
        ''' 1 variant 1 exec for palallel
        '''

        # Read the data from the 010_variant file
        self._prepare_from_variant(variant_df_row, distin_dict)

        # skip if pick_mode is different
        if not utl.is_my_pick_mode(self.var_type, distin_dict['pick_mode']):
            return

        # Cut out two around_seq before and after the REF variant.
        self._pick_ref_around_seq()

        # Product()
        # Hold related sequences for each group to be compared.
        for gno in range(2):
            ano = self.ano_gno[gno]
            self.gr_product[gno].set_info(
                gno,
                ano,
                self.gname_gno[gno],
                self.chrom,
                self.pos,
                self.var_type,
                self.vseq_ano[ano],
                self.vseq_lens_ano[ano],
                self.around_seq_pre,
                self.around_seq_aft,
            )

        # If var_type is indel, it is decided to be available
        if self.var_type == glv.INDEL:
            self.marker_available = True
        else:
            # If var_type is not indel, after checking caps
            self.marker_available = self._check_caps()


    def _check_caps(self):
        '''
        '''

        marker_available = False

        caps_checked_list = list()

        # If you need a detailed report of the enzyme,
        # Specify the command parameter analyse_caps as True.
        mes_list = list()
        mes_list = ["\ntest start\n"]

        for gno in range(2):
            mes_list += self._analong_header(gno, mes_list)

            # Check the effect of enzyme on two group sequences
            # caps_check_dict
            # caps_check_dict[enzyme_string] = {
            #     'ResType': enzyme_RestrictionType,
            #     'result': caps_ResTyp_dict[enzyme_RestrictionType],
            # }
            caps_check_dict, \
            enzyme_map_txt = \
                self._check_effect_of_enzyme(
                    self.gr_product[gno].seq_target,
                    self.enzyme_name_list)

            mes_list += ["{}".format(enzyme_map_txt)]
            caps_checked_list.append(caps_check_dict)

        mes_list += ["\ntest end\n"]

        #log.debug("\n{}".format(pprint.pformat(mes_list)))

        # Compare the result of caps_check_dict to determine
        # if it is available in caps
        self.caps_found = \
            self._compare_digestion_of_enzymes(caps_checked_list)

        if self.caps_found > 0:
            marker_available = True
            mes_list += ["caps is available {}".format(self.caps_found)]

            for enzyme in self.caps_result.keys():
                enz_info = "{} {} {} {} {}".format(
                    enzyme,
                    self.caps_result[enzyme]['digested_gno'],
                    self.caps_result[enzyme]['found_pos'],
                    self.caps_result[enzyme]['digest_pattern'],
                    self.caps_result[enzyme]['digested_pos'])

            mes_list += [enz_info]

        else:    
            marker_available = False
            mes_list += ["caps is not available {}".format(self.caps_found)]

        mes_list += ["<><><><><><><><><><><><><><><><><>\n\n"]


        if glv.conf.analyse_caps == True:
            mes_line = "\n".join(mes_list)
            log.info("\n{}".format(mes_line))

        return marker_available


    def _compare_digestion_of_enzymes(self, caps_checked_list):
        '''
        '''

        # caps_checked_list: Digest status of each group
        list_0 = list(caps_checked_list[0].keys())
        list_1 = list(caps_checked_list[1].keys())


        #pprint.pprint(caps_checked_list)
        #print(list_0)
        #print(list_1)

        # xor list
        uniq_list = list(set(list_0) ^ set(list_1))

        # type(BsmFI) RestrictionType
        for enzyme_name in uniq_list:

            digested_gno = 0

            if enzyme_name in list_1:
                digested_gno = 1

            enz_res_type = \
                caps_checked_list[digested_gno][enzyme_name]['ResType']
            result_list = \
                caps_checked_list[digested_gno][enzyme_name]['res_list']

            digest_pattern = enz_res_type.elucidate()
            found_pos = int(result_list[0])

            # digested_pos
            digested_pos = 0
            pattern_5_digest = ""

            if "_" in digest_pattern:
                pattern_5_digest = re.sub(r"_", "", digest_pattern)
                pos_0 = pattern_5_digest.find("^")
                if pos_0 != -1:
                    # Last position remaining on the 5'side
                    # AAT^_ATT | AAT^ATT -> 3
                    digested_pos = pos_0

            # new
            self.caps_result[enzyme_name] = {
                'digested_gno':     digested_gno,
                'found_pos':        found_pos,
                'digest_pattern':   digest_pattern,
                'digested_pos':     digested_pos}
            #log.debug("{} | {} -> {}".format(
            #    digest_pattern, pattern_5_digest, digested_pos))

        caps_found = len(self.caps_result)
        return caps_found


    def _check_effect_of_enzyme(self, seq_target, enzyme_name_list):
        ''' http://biopython.org/DIST/docs/cookbook/Restriction.html
        biopython <= 1.76 for IUPACAmbiguousDNA()
        '''

        caps_ResTyp_dict = dict()
        caps_check_dict = dict()
        enzyme_map_txt = ""

        # 4.1 Setting up an Analysis
        # 4.2 Full restriction analysis
        multi_site_seq = Seq(seq_target, IUPACAmbiguousDNA())
        rb = Restriction.RestrictionBatch(enzyme_name_list)
        Analong = Restriction.Analysis(rb, multi_site_seq)

        # 4.5 Fancier restriction analysis
        #
        # full()
        #   all the enzymes in the RestrictionBatch
        #   {KpnI: [], EcoRV: [], EcoRI: [33]}
        # with_sites()
        #   output only the result for enzymes which have a site
        #   result_dict = {EcoRI: [33]}

        caps_ResTyp_dict = Analong.with_sites()

        # make dictionary as string enzyme name
        for enzyme_RestrictionType in caps_ResTyp_dict.keys():
            enzyme_string = str(enzyme_RestrictionType)

            # caps_check_dict
            caps_check_dict[enzyme_string] = {
                'ResType': enzyme_RestrictionType,
                'res_list': caps_ResTyp_dict[enzyme_RestrictionType],
            }

        # detail information: make a restriction map of a sequence
        if glv.conf.analyse_caps == True:
            Analong.print_as('map')
            enzyme_map_txt_all = Analong.format_output()
            enzyme_map_txt = ""

            for line in enzyme_map_txt_all.split('\n'):
                if " Enzymes which " in line:
                    break
                enzyme_map_txt += "{}\n".format(line)

            enzyme_map_txt += "caps_check_dict={}".format(
                caps_check_dict)

        return caps_check_dict, \
            enzyme_map_txt

            
    def _analong_header(self, gno, mes_list):
        '''
        '''

        m_list = list()
        m_list = ["gno={}, chrom={}, pos={}".format(gno,
            self.gr_product[gno].chrom, self.gr_product[gno].pos)]
        m_list += ["seq_target={}".format(self.gr_product[gno].seq_target)]
        m_list += ["seq_target_len={}".format(
            len(self.gr_product[gno].seq_target))]

        return m_list


    def _prepare_from_variant(self, variant_df_row, distin_dict):
        ''' Read the data from the 010_variant file, save it locally,
            and split the variables as needed
        '''

        hdr_dict = distin_dict['variant']['hdr_dict']

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
        self.set_n, \
        self.marker_info, \
        self.vseq_lens_ano_str, \
        self.enzyme_name, \
        self.digest_pattern, \
        self.target_gno, \
        self.target_len = \
            utl.get_basic_primer_info(variant_df_row, hdr_dict)

        # chrom pos
#        self.chrom = str(variant_df_row[hdr_dict['chrom']])
#        self.pos = int(variant_df_row[hdr_dict['pos']])

        # c2,c3
#        self.targ_grp = str(variant_df_row[hdr_dict['targ_grp']])
        # g0_name g1_name
#        self.g0_name, self.g1_name = self.targ_grp.split(',')
        # list
        self.gname_gno = [self.g0_name, self.g1_name]

        # targ_ano = ano_corresponding_to_g0_g1
#        self.targ_ano = str(variant_df_row[hdr_dict['targ_ano']])
#        self.g0_ano, self.g1_ano = map(int, self.targ_ano.split(','))
        # convert table
        self.ano_gno = [self.g0_ano, self.g1_ano]

        # vseq_gno_str
#        self.vseq_gno_str = str(variant_df_row[hdr_dict['vseq_gno_str']])

        # 00/01,hohe_s1,1.1/1.1
#        self.gts_segr_lens = str(variant_df_row[hdr_dict['gts_segr_lens']])

        # 1/1
#        self.set_n = str(variant_df_row[hdr_dict['set_n']])
        # set_my_no, set_total_no
        self.set_my_no, self.set_total_no = map(int, self.set_n.split('/'))

        # 1,1,0,-1
        self.len_g0g1_dif_long = str(
            variant_df_row[hdr_dict['len_g0g1_dif_long']])
        # g0_len, g1_len, diff_len, long_gno
        self.g0_len, self.g1_len, self.diff_len, self.long_gno = \
            map(int, self.len_g0g1_dif_long.split(','))

#        self.var_type = str(variant_df_row[hdr_dict['var_type']])

        # T,G access to all vseq
        self.vseq_ano_str = str(variant_df_row[hdr_dict['vseq_ano_str']])
        self.vseq_ano = self.vseq_ano_str.split(',')
        self.vseq_lens_ano = [len(vseq_ano) for vseq_ano in self.vseq_ano]

        # short_gno
        if self.long_gno == glv.SAME_LENGTH:
            self.long_gno = 0
            self.short_gno = 1
        else:
            self.short_gno = 0 if self.long_gno == 1 else 1

        self.longer_len = self.vseq_lens_ano[0]
        self.shorter_len = self.vseq_lens_ano[1]


    def _pick_ref_around_seq(self):
        ''' Cut out two around_seq before and after the REF variant.
        '''

        #-------------------------------------
        # around_seq is based on vseq_ano[0]
        self.vseq_ref = self.vseq_ano[0]
        #-------------------------------------
        self.vseq_ref_len = len(self.vseq_ref)

        #print("self.pos={}".format(self.pos))
        #print("self.vseq_ref={}".format(self.vseq_ref))
        #print("self.vseq_ref_len={}".format(len(self.vseq_ref)))

        # abs pos決め
        # pos=60, 60-1=59
        self.abs_around_seq_pre_end = self.pos - 1
        #print("self.abs_around_seq_pre_end={}".format(
        #   self.abs_around_seq_pre_end))

        # 59-10+1 = 50
        self.abs_around_seq_pre_stt = \
            self.abs_around_seq_pre_end - glv.AROUND_SEQ_LEN + 1
        #print("glv.AROUND_SEQ_LEN={}".format(glv.AROUND_SEQ_LEN))
        #print("self.abs_around_seq_pre_stt={}".format(
        #   self.abs_around_seq_pre_stt))

        # 60+16 = 76
        # 60+1 = 61
        self.abs_around_seq_aft_stt = self.pos + self.vseq_ref_len
        #print("self.vseq_ref_len={}".format(self.vseq_ref_len))
        #print("self.abs_around_seq_aft_stt={}".format(
        #   self.abs_around_seq_aft_stt))

        # 76+10-1=85
        self.abs_around_seq_aft_end = \
            self.abs_around_seq_aft_stt + glv.AROUND_SEQ_LEN - 1
        #print("self.abs_around_seq_aft_end={}".format(
        #   self.abs_around_seq_aft_end))

        # preの切り出し
        self.around_seq_pre = glv.ref.pick_refseq(
            self.chrom,
            self.abs_around_seq_pre_stt,
            self.abs_around_seq_pre_end).upper()

        #print("self.around_seq_pre={}".format(self.around_seq_pre))

        # aftの切り出し
        self.around_seq_aft = glv.ref.pick_refseq(
            self.chrom,
            self.abs_around_seq_aft_stt,
            self.abs_around_seq_aft_end).upper()
        #print("self.around_seq_aft={}".format(self.around_seq_aft))

        self.around_seq_pre_len = len(self.around_seq_pre)
        #print("self.around_seq_pre_len={}".format(self.around_seq_pre_len))
        self.around_seq_aft_len = len(self.around_seq_aft)
        #print("self.around_seq_aft_len={}".format(self.around_seq_aft_len))

        self.seq_target_ref = "{}{}{}".format(
            self.around_seq_pre,
            self.vseq_ref,
            self.around_seq_aft)

        #print("self.seq_target_ref={}".format(self.seq_target_ref))

        self.seq_target_ref_len = len(self.seq_target_ref)

        #print("self.seq_target_ref_len={}".format(self.seq_target_ref_len))


    def make_seq_template_ref(self):
        '''
        '''

        #self.seq_template_ref = ''
        #self.seq_template_ref_len = 0

        # frag_padを切り出す
        self.fragment_pad_len = glv.conf.fragment_pad_len

        # abs_posを決める
        self.abs_frag_pad_pre_end = self.abs_around_seq_pre_stt - 1

        self.abs_frag_pad_pre_stt = \
            self.abs_frag_pad_pre_end - self.fragment_pad_len + 1

        self.abs_frag_pad_aft_stt = self.abs_around_seq_aft_end + 1
        self.abs_frag_pad_aft_end = \
            self.abs_frag_pad_aft_stt + self.fragment_pad_len - 1

        # templateの絶対pos stringを最初に作る。
        # 最初はposはすべて完成しているが、今後端を切るために。
        self.seq_template_ref_abs_pos = \
            self._make_abs_pos(
                self.abs_frag_pad_pre_stt,
                self.abs_frag_pad_aft_end)

        # これは最初の
        # templateの相対pos string
        self.seq_template_ref_rel_pos, \
        self.SEQUENCE_TARGET = \
            self._convert_to_rel_pos(self.seq_template_ref_abs_pos)

        # refのfragpadを取り出す。
        self.frag_pad_pre, self.frag_pad_aft, self.seq_template_ref = \
            self._get_seq_template_ref(
                self.chrom, self.seq_template_ref_abs_pos)

        self.frag_pad_pre_len = len(self.frag_pad_pre)
        self.frag_pad_aft_len = len(self.frag_pad_aft)
        self.seq_template_ref_len = len(self.seq_template_ref)

        # groupのproductにもセットする
        for gno in range(2):
            self.gr_product[gno].set_frag_pad(
                self.frag_pad_pre, self.frag_pad_aft)


    def _make_abs_pos(
        self,
        abs_frag_pad_pre_stt,
        abs_frag_pad_aft_end):

        return "{}/{}/{}/{}/{}/{}/{}/{}/{}".format(
            abs_frag_pad_pre_stt,
            self.abs_frag_pad_pre_end,
            self.abs_around_seq_pre_stt,
            self.abs_around_seq_pre_end,
            self.pos,
            self.abs_around_seq_aft_stt,
            self.abs_around_seq_aft_end,
            self.abs_frag_pad_aft_stt,
            abs_frag_pad_aft_end)


    def _convert_to_rel_pos(self, seq_template_ref_abs_pos):

        abs_frag_pad_pre_stt, abs_frag_pad_pre_end, \
        abs_around_seq_pre_stt, abs_around_seq_pre_end, \
        abs_pos, \
        abs_around_seq_aft_stt, abs_around_seq_aft_end, \
        abs_frag_pad_aft_stt, abs_frag_pad_aft_end = \
            self._separate_pos_str(seq_template_ref_abs_pos)

        # 11......21
        # 1.......11
        dif = abs_frag_pad_pre_stt

        rel_frag_pad_pre_stt = 1
        rel_frag_pad_pre_end = abs_frag_pad_pre_end - dif + 1
        rel_around_seq_pre_stt = abs_around_seq_pre_stt - dif + 1
        rel_around_seq_pre_end = abs_around_seq_pre_end - dif + 1
        rel_pos = abs_pos - dif + 1
        rel_around_seq_aft_stt = abs_around_seq_aft_stt - dif + 1
        rel_around_seq_aft_end = abs_around_seq_aft_end - dif + 1
        rel_frag_pad_aft_stt = abs_frag_pad_aft_stt - dif + 1
        rel_frag_pad_aft_end = abs_frag_pad_aft_end - dif + 1

        #log.debug("{} {} {} {}".format(
        #    self.SEQUENCE_TARGET,
        #    rel_around_seq_pre_stt,
        #    rel_around_seq_aft_end,
        #    rel_around_seq_pre_stt))

        # 100-1=99+1=100
        sequence_target = "{},{}".format(
            rel_around_seq_pre_stt,
            rel_around_seq_aft_end - rel_around_seq_pre_stt + 1)

        return "{}/{}/{}/{}/{}/{}/{}/{}/{}".format(
            rel_frag_pad_pre_stt,
            rel_frag_pad_pre_end,
            rel_around_seq_pre_stt,
            rel_around_seq_pre_end,
            rel_pos,
            rel_around_seq_aft_stt,
            rel_around_seq_aft_end,
            rel_frag_pad_aft_stt,
            rel_frag_pad_aft_end), \
            sequence_target


    def _separate_pos_str(self, pos_str):

        return map(int, pos_str.split('/'))


    def _get_seq_template_ref(self, chrom, seq_template_ref_abs_pos):

        abs_frag_pad_pre_stt, abs_frag_pad_pre_end, \
        abs_around_seq_pre_stt, abs_around_seq_pre_end, \
        abs_pos, \
        abs_around_seq_aft_stt, abs_around_seq_aft_end, \
        abs_frag_pad_aft_stt, abs_frag_pad_aft_end = \
            self._separate_pos_str(seq_template_ref_abs_pos)

        # update
        # pick frag_pad_pre
        frag_pad_pre = glv.ref.pick_refseq(
            chrom,
            abs_frag_pad_pre_stt,
            abs_frag_pad_pre_end).upper()

        # これは変わらない
        seq_target_ref = self.seq_target_ref

        # pick frag_pad_aft
        frag_pad_aft = glv.ref.pick_refseq(
            chrom,
            abs_frag_pad_aft_stt,
            abs_frag_pad_aft_end).upper()

        seq_template_ref = "{}{}{}".format(
            frag_pad_pre,
            seq_target_ref,
            frag_pad_aft)

        return frag_pad_pre, frag_pad_aft, seq_template_ref


    def copy_line_for_effective_restriction_enzymes(self):
        '''
        '''

        enzyme_cnt_per_variant = 0

        l_marker_id = list()
        l_marker_info = list()
        l_seq_t_abs_pos = list()
        l_seq_t_rel_pos = list()
        l_seq_t_ref = list()
        l_seq_target = list()

        if self.var_type == glv.INDEL:

            enzyme_cnt_per_variant = 1

            marker_id = self._make_marker_id()
            l_marker_id = [marker_id]

            marker_info = self._make_marker_info()
            l_marker_info = [marker_info]

            l_seq_t_abs_pos = [self.seq_template_ref_abs_pos]
            l_seq_t_rel_pos = [self.seq_template_ref_rel_pos]
            l_seq_t_ref = [self.seq_template_ref]
            l_seq_target = [self.SEQUENCE_TARGET]

        else:
            # duplicate line information by enzyme
            enzyme_cnt_per_variant = len(self.caps_result)

            # Build information as a multi-line list.
            self._divide_information_by_enzyme(
                enzyme_cnt_per_variant,
                l_marker_id,
                l_marker_info,
                l_seq_t_abs_pos,
                l_seq_t_rel_pos,
                l_seq_t_ref,
                l_seq_target)
        
        line_for_each_enzyme = list()        


        for num in range(enzyme_cnt_per_variant):
            enzyme_cnt = "{}/{}".format(num+1, enzyme_cnt_per_variant)
            set_enz_cnt = "{}-{}".format(self.set_n, enzyme_cnt)
            vseq_lens_ano_str = \
                "{}".format(','.join(map(str, self.vseq_lens_ano)))

            l_list = list()

            #print(l_marker_id)
            #print(l_marker_info)
            #print(enzyme_cnt_per_variant)
            #print("")

            # out to marker out file
            # Synchronize with eval_variant.py outlist.py
            l_list += [l_marker_id[num]]
            l_list += [self.chrom]
            l_list += [self.pos]
            l_list += [self.targ_grp]
            l_list += [self.targ_ano]
            l_list += [self.vseq_gno_str]
            l_list += [self.gts_segr_lens]
            l_list += [self.var_type]
            # ----------------------
            l_list += [set_enz_cnt]
            # ----------------------
            l_list += [l_marker_info[num]]
            # vseq_lens_ano
            l_list += [vseq_lens_ano_str]
            # ----------------------

            # 1) g0_seq_target_len
            l_list += [self.gr_product[0].seq_target_len]
            # 2) g0_seq_target
            l_list += [self.gr_product[0].seq_target]

            # 1) g1_seq_target_len
            l_list += [self.gr_product[1].seq_target_len]
            # 2) g1_seq_target
            l_list += [self.gr_product[1].seq_target]

            # 1) seq_template_ref_len
            l_list += [len(l_seq_t_ref[num])]
            # 2) seq_template_ref_abs_pos
            l_list += [l_seq_t_abs_pos[num]]
            # 3) seq_template_ref_rel_pos
            l_list += [l_seq_t_rel_pos[num]]
            # 4) SEQUENCE_TARGET
            l_list += [l_seq_target[num]]
            # 5) seq_template_ref
            l_list += [l_seq_t_ref[num]]

            line_for_each_enzyme.append('\t'.join(map(str, l_list)))

        self.line = '\n'.join(map(str, line_for_each_enzyme))


    def _divide_information_by_enzyme(
        self,
        enzyme_cnt_per_variant,
        l_marker_id,
        l_marker_info,
        l_seq_t_abs_pos,
        l_seq_t_rel_pos,
        l_seq_t_ref,
        l_seq_target):
        '''
        '''

        #print(self.caps_result)

        #=============~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for enzyme in self.caps_result.keys():

            #print("enzyme={}".format(enzyme))

            # for marker id
            marker_id = self._make_marker_id(enzyme)
            l_marker_id.append(marker_id)

            # enzymeごとに、digest_posは、grに組み込む
            marker_info = self._make_marker_info(enzyme)
            l_marker_info.append(marker_info)

            # 何もなければすでにrefとしてセットされた値を
            seq_t_abs_pos = self.seq_template_ref_abs_pos
            seq_t_rel_pos = self.seq_template_ref_rel_pos
            seq_t_ref = self.seq_template_ref
            SEQUENCE_TARGET = self.SEQUENCE_TARGET

            text = ""
            caps_check_dict, text = \
                self._check_effect_of_enzyme(seq_t_ref, [enzyme])

            caps_result_cnt = len(caps_check_dict)
            #print(caps_result_cnt)
            #pprint.pprint(caps_check_dict)

            # 自分自身を除いて処理
            if caps_result_cnt != 0:
                # 切断されているなら、seq_templateを更新し、
                # 現在、
                # ref.
                # self.gr_product[num]

                # 20200715
                # それぞれ、メソッドで使っているだけ
                seq_t_abs_pos = \
                    self._change_abs_pos_by_digest(
                        caps_check_dict[enzyme]['res_list'])

                # templateの相対pos string
                seq_t_rel_pos, SEQUENCE_TARGET = \
                    self._convert_to_rel_pos(seq_t_abs_pos)

                # まとめてpick
                frag_pad_pre, frag_pad_aft, seq_t_ref = \
                    self._get_seq_template_ref(
                        self.chrom, seq_t_abs_pos)

                # 確認用 capsを検索する
                #caps_result_dict, caps_result_dict_str = \
                #    Products.get_caps_result(
                #        seq_t_ref, [enzyme])
                #log.debug("rechecked {}".format(caps_result_dict_str))
                #log.debug("{}\n".format(seq_t_rel_pos))

            else:
                pass

            l_seq_t_abs_pos.append(seq_t_abs_pos)
            l_seq_t_rel_pos.append(seq_t_rel_pos)
            l_seq_t_ref.append(seq_t_ref)
            l_seq_target.append(SEQUENCE_TARGET)


    def _make_marker_info(self, enzyme_name=''):
        '''
        '''

        marker_info = ''

        # INDEL
        if self.var_type == glv.INDEL:
            longer_group = self.long_gno
            longer_length = self.longer_len
            shorter_length = self.shorter_len
            diff_length = self.diff_len
            digested_pos = 0

            marker_info = "{}.{}.{}.{}.{}".format(
                longer_group,
                longer_length,
                shorter_length,
                diff_length,
                digested_pos)

        # not INDEL
        else:

            marker_info = "{}.{}.{}.{}.{}".format(
                enzyme_name,
                self.caps_result[enzyme_name]['digested_gno'],
                self.caps_result[enzyme_name]['found_pos'],
                self.caps_result[enzyme_name]['digest_pattern'],
                self.caps_result[enzyme_name]['digested_pos'])

        return marker_info


    @classmethod
    def split_marker_info(cls, marker_info):
        '''
        '''

        # for indel
            # longer_group, 
            # longer_length,
            # shorter_length,
            # diff_length,
            # digested_pos

        # for !indel
            # enzyme_name,      str
            # digested_gno,
            # found_pos,
            # digest_pattern,   str
            # digested_pos

        return marker_info.split(".")


    def _make_marker_id(self, enzyme_name=''):
        '''
        '''

        marker_id = ""

        # 1.chrom
        chrom = self.chrom
        # 2.pos
        pos = self.pos
        # 3. Allele number corresponding to groups 0 and 1
        ano_corresponding_to_g0_g1 = self.targ_ano
        # 4. variant type
        var_type = self.var_type

        # 5.SEQUENCE_TARGET detail
        if var_type == glv.INDEL:
            # 5.1 The longer group
            longer_group = self.long_gno
            # 5.2 The longer one, the length
            longer_length = self.longer_len
            # 5.3 The shorter one, the length
            shorter_length = self.shorter_len

            seq_target_detail = "{}.{}.{}".format(
                longer_group,
                longer_length,
                shorter_length)

        else:
            # snp or else
            # 5.1 digested group
            digested_gno = self.caps_result[enzyme_name]['digested_gno']

            # 5.2 found point
            found_pos = self.caps_result[enzyme_name]['found_pos']

            # 5.3 restriction enzyme
            # enzyme_name

            seq_target_detail = "{}.{}.{}".format(
                digested_gno,
                found_pos,
                enzyme_name)

        marker_id = "{}.{}.{}.{}.{}".format(
            chrom,
            pos,
            ano_corresponding_to_g0_g1,
            var_type,
            seq_target_detail)

        return marker_id


    def _change_abs_pos_by_digest(self, caps_dig_pos):

        # sequence_target の外側にある digested_pointだけが
        # 変更対象である。

        fixed_pre_stt = self.abs_frag_pad_pre_stt
        fixed_pre_end = self.abs_frag_pad_pre_end
        fixed_aft_stt = self.abs_frag_pad_aft_stt
        fixed_aft_end = self.abs_frag_pad_aft_end

        five_prime_biggest_pos = fixed_pre_stt
        three_prime_smallest_pos = fixed_aft_end

        # EcoRI': {'digest_pattern': 'G^AATT_C'}
        # [17]
        # 123456789012345678901234567890
        # CTCTGTTCGGTGGAAGAATTCAGATTTCAGAGTCA
        #               G^AATT_C
        #                 /-> 17
        #                 AATTCAGATTTCAGAGTCA
        # 切断されるポイントは17。これは残る側。
        # 残る側に、切断ポイントを残していいのか。
        # 今は残している。

        # digest_positionごとに調査
        for rel_digest_pos in caps_dig_pos:
            # 絶対posに変換
            # 10001    100 -> 10001+100 - 10100
            abs_digest_pos = fixed_pre_stt + rel_digest_pos - 1

            #log.debug("rel={} abs={} pre_stt<{} pre_end>{} aftstt<{}".format(
            #    rel_digest_pos,
            #    abs_digest_pos,
            #    fixed_pre_stt,
            #    fixed_pre_end,
            #    fixed_aft_stt))

            # |fixed_pre_stt
            #                |fixed_pre_end
            #                               |fixed_aft_stt
            # <--------------><=============<-------------->

            # １つポジションをずらす。それにより認識サイトが
            # 壊れる
            if abs_digest_pos < fixed_pre_end:
                # 5'側では、一つ先で切る
                # always
                five_prime_biggest_pos = abs_digest_pos + 1

            elif fixed_aft_stt < abs_digest_pos:
                # only once
                # 3'側では、一つ手前で切る
                three_prime_smallest_pos = abs_digest_pos - 1
                break

        seq_template_ref_abs_pos = \
            self._set_seq_template_ref_abs_pos(
                five_prime_biggest_pos,
                three_prime_smallest_pos)

        #log.debug("{}".format(self.seq_template_ref_abs_pos))
        #log.debug("{}".format(seq_template_ref_abs_pos))

        #sys.exit(1)

        return seq_template_ref_abs_pos


    def _get_seq_template_ref(self, chrom, seq_template_ref_abs_pos):

        abs_frag_pad_pre_stt, abs_frag_pad_pre_end, \
        abs_around_seq_pre_stt, abs_around_seq_pre_end, \
        abs_pos, \
        abs_around_seq_aft_stt, abs_around_seq_aft_end, \
        abs_frag_pad_aft_stt, abs_frag_pad_aft_end = \
            self._separate_pos_str(seq_template_ref_abs_pos)

        # update
        # pick frag_pad_pre
        frag_pad_pre = glv.ref.pick_refseq(
            chrom,
            abs_frag_pad_pre_stt,
            abs_frag_pad_pre_end).upper()

        # これは変わらない
        seq_target_ref = self.seq_target_ref

        # pick frag_pad_aft
        frag_pad_aft = glv.ref.pick_refseq(
            chrom,
            abs_frag_pad_aft_stt,
            abs_frag_pad_aft_end).upper()

        seq_template_ref = "{}{}{}".format(
            frag_pad_pre,
            seq_target_ref,
            frag_pad_aft)

        return frag_pad_pre, frag_pad_aft, seq_template_ref


    def _set_seq_template_ref_abs_pos(
        self,
        abs_frag_pad_pre_stt,
        abs_frag_pad_aft_end):

        return "{}/{}/{}/{}/{}/{}/{}/{}/{}".format(
            abs_frag_pad_pre_stt,
            self.abs_frag_pad_pre_end,
            self.abs_around_seq_pre_stt,
            self.abs_around_seq_pre_end,
            self.pos,
            self.abs_around_seq_aft_stt,
            self.abs_around_seq_aft_end,
            self.abs_frag_pad_aft_stt,
            abs_frag_pad_aft_end)





