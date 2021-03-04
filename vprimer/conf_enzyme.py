# -*- coding: utf-8 -*-

# http://rebase.neb.com/rebase/rebase.enz.html
# https://international.neb.com/tools-and-resources/selection-charts/isoschizoomers

import sys
import os
import errno
import pprint

# global configuration
import vprimer.glv as glv
import vprimer.utils as utl

from vprimer.logging_config import LogConf
import Bio.Restriction.Restriction_Dictionary as ResDict

class ConfEnzyme(object):

    def open_log_enzyme(self):

        # in glv
        global log
        log = LogConf.open_log(__name__)

        self.default_enzyme_names = [
        ['# A'],
            ['AluI',        'AGCT',             'AG^_CT',           '4'],
            ['ApaI',        'GGGCCC',           'G_GGCC^C',         '6'],
            ['AscI',        'GGCGCGCC',         'GG^CGCG_CC',       '8'],
            ['AvrII',       'CCTAGG',           'C^CTAG_G',         '6'],
            [''],
        ['# B'],
            ['BamHI',       'GGATCC',           'G^GATC_C',         '6'],
            ['BbsI',        'GAAGAC',           'GAAGACNN^NNNN_N',  '6'],
            ['BclI',        'TGATCA',           'T^GATC_A',         '6'],
            ['BglII',       'AGATCT',           'A^GATC_T',         '6'],
            ['BsaI',        'GGTCTC',           'GGTCTCN^NNNN_N',   '6'],
            ['BsiWI',       'CGTACG',           'C^GTAC_G',         '6'],
            ['BsmFI',       'GGGAC',            'GGGACNNNNNNNNNN^NNNN_N','5'],
            ['BspHI',       'TCATGA',           'T^CATG_A',         '6'],
            ['BssHII',      'GCGCGC',           'G^CGCG_C',         '6'],
            ['Bst1107I',    'GTATAC',           'GTA^_TAC',         '6'],
            ['BstBI',       'TTCGAA',           'TT^CG_AA',         '6'],
            ['BstEII',      'GGTNACC',          'G^GTNAC_C',        '7'],
            ['BstXI',       'CCANNNNNNTGG',     'CCAN_NNNN^NTGG',   '12'],
            [''],
        ['# C'],
            ['ClaI',        'ATCGAT',           'AT^CG_AT',         '6'],
            [''],
        ['# D'],
            ['DdeI',        'CTNAG',            'C^TNA_G',          '5'],
            ['DpnI',        'GATC',             'GA^_TC',           '4'],
            ['DraI',        'TTTAAA',           'TTT^_AAA',         '6'],
            ['DraIII',      'CACNNNGTG',        'CAC_NNN^GTG',      '9'],
            [''],
        ['# E'],
            ['Eco52I',      'CGGCCG',           'C^GGCC_G',         '6'],
            ['EcoO109I',    'RGGNCCY',          'RG^GNC_CY',        '7'],
            ['EcoO65I',     'GGTNACC',          'G^GTNAC_C',        '7'],
            ['EcoRI',       'GAATTC',           'G^AATT_C',         '6'],
            ['EcoRV',       'GATATC',           'GAT^_ATC',         '6'],
            ['EcoT14I',     'CCWWGG',           'C^CWWG_G',         '6'],
            [''],
        ['# F'],
            ['FseI',        'GGCCGGCC',         'GG_CCGG^CC',       '8'],
            [''],
        ['# H'],
            ['HaeII',       'RGCGCY',           'R_GCGC^Y',         '6'],
            ['HincII',      'GTYRAC',           'GTY^_RAC',         '6'],
            ['HindIII',     'AAGCTT',           'A^AGCT_T',         '6'],
            ['HinfI',       'GANTC',            'G^ANT_C',          '5'],
            ['HpaI',        'GTTAAC',           'GTT^_AAC',         '6'],
            ['HphI',        'GGTGA',            'GGTGANNNNNNN_N^N', '5'],
            [''],
        ['# K'],
            ['KpnI',        'GGTACC',           'G_GTAC^C',         '6'],
            [''],
        ['# M'],
            ['MluI',        'ACGCGT',           'A^CGCG_T',         '6'],
            ['MseI',        'TTAA',             'T^TA_A',           '4'],
            [''],
        ['# N'],
            ['NcoI',        'CCATGG',           'C^CATG_G',         '6'],
            ['NdeI',        'CATATG',           'CA^TA_TG',         '6'],
            ['NheI',        'GCTAGC',           'G^CTAG_C',         '6'],
            ['NlaIII',      'CATG',             '_CATG^',           '4'],
            ['NotI',        'GCGGCCGC',         'GC^GGCC_GC',       '8'],
            ['NruI',        'TCGCGA',           'TCG^_CGA',         '6'],
            ['NsiI',        'ATGCAT',           'A_TGCA^T',         '6'],
            [''],
        ['# P'],
            ['PacI',        'TTAATTAA',         'TTA_AT^TAA',       '8'],
            ['PmeI',        'GTTTAAAC',         'GTTT^_AAAC',       '8'],
            ['PmlI',        'CACGTG',           'CAC^_GTG',         '6'],
            ['Psp1406I',    'AACGTT',           'AA^CG_TT',         '6'],
            ['PstI',        'CTGCAG',           'C_TGCA^G',         '6'],
            ['PvuII',       'CAGCTG',           'CAG^_CTG',         '6'],
            [''],
        ['# R'],
            ['RsaI',        'GTAC',             'GT^_AC',           '4'],
            [''],
        ['# S'],
            ['SacI',        'GAGCTC',           'G_AGCT^C',         '6'],
            ['SacII',       'CCGCGG',           'CC_GC^GG',         '6'],
            ['SalI',        'GTCGAC',           'G^TCGA_C',         '6'],
            ['SapI',        'GCTCTTC',          'GCTCTTCN^NNN_N',   '7'],
            ['SbfI',        'CCTGCAGG',         'CC_TGCA^GG',       '8'],
            ['ScaI',        'AGTACT',           'AGT^_ACT',         '6'],
            ['SfiI',        'GGCCNNNNNGGCC',    'GGCCN_NNN^NGGCC',  '13'],
            ['SmaI',        'CCCGGG',           'CCC^_GGG',         '6'],
            ['SnaBI',       'TACGTA',           'TAC^_GTA',         '6'],
            ['SpeI',        'ACTAGT',           'A^CTAG_T',         '6'],
            ['SphI',        'GCATGC',           'G_CATG^C',         '6'],
            ['SspI',        'AATATT',           'AAT^_ATT',         '6'],
            ['StuI',        'AGGCCT',           'AGG^_CCT',         '6'],
            ['SwaI',        'ATTTAAAT',         'ATTT^_AAAT',       '8'],
            [''],
        ['# T'],
            ['TaqI',        'TCGA',             'T^CG_A',           '4'],
            ['Tth111I',     'GACNNNGTC',        'GACN^N_NGTC',      '9'],
            [''],
        ['# X'],
            ['XbaI',        'TCTAGA',           'T^CTAG_A',         '6'],
            ['XhoI',        'CTCGAG',           'C^TCGA_G',         '6'],
            ['XmaI',        'CCCGGG',           'C^CCGG_G',         '6'],
            [''],
        ['# end'],
        ]


    def read_enzyme_file(self):
        '''
        '''
        enzyme_files_list = list()

        basename_enzyme = os.path.basename(self.enzyme_files_user_list[0])

        # When not specified by the user
        if basename_enzyme == "no_enzyme":
            new_enzyme_file = "{}/{}".format(
                self.ref_dir_path, "enzyme_names.txt")
            new_enzyme_file_path = utl.full_path(new_enzyme_file)

            # default enzyme names file
            if os.path.isfile(new_enzyme_file_path):
                mes = "Since the enzyme names file is not specified, "
                mes += "set {} as enzyme names file.".format(
                    new_enzyme_file_path)
                log.info(mes)

            else:
                mes = "Since the enzyme names file is not specified, "
                mes += "make new {} as enzyme names file.".format(
                    new_enzyme_file_path)
                log.info(mes)
                self._write_new_enzyme_names_file(new_enzyme_file_path)
                
            # set
            self.enzyme_files_user_list.clear()
            self.enzyme_files_user_list.append(new_enzyme_file_path)


        for enzyme_file_user in self.enzyme_files_user_list:

            # enzyme_file_user
            if os.path.isfile(enzyme_file_user):
                log.info("found {}.".format(enzyme_file_user))

                basename_user = os.path.basename(enzyme_file_user)
                enzyme_file_slink_system = "{}/{}{}".format(
                    self.ref_dir_path, "slink_", basename_user)
                enzyme_files_list.append(enzyme_file_slink_system)

                # enzyme_file_slink_system
                if os.path.isfile(enzyme_file_slink_system):
                    log.info("found {}.".format(enzyme_file_slink_system))
                else:
                    log.info("not found {}.".format(enzyme_file_slink_system))
                    utl.ln_s(enzyme_file_user, enzyme_file_slink_system)

            else:
                log.info("not found {}. exit.".format(enzyme_file_user))
                sys.exit(1)

        enzyme_name_list = list()

        for enzyme_file in enzyme_files_list:
            log.info("open {}".format(enzyme_file))
            with open(enzyme_file, "r") as f:
                # iterator
                for r_liner in f:
                    r_line = r_liner.strip()    # cr, ws

                    if r_line.startswith('#') or r_line == '':
                        continue
                    r_line = utl.strip_hash_comment(r_line)
                    # for enzyme name (line top)
                    l_list = r_line.split()
                    enzyme_name_list.append(l_list[0])

        # add hand enzyme
        if self.enzyme_str != "":
            hand_enzyme_list = self.enzyme_str.split(',')
            hand_cnt = len(hand_enzyme_list)
            enzyme_name_list += hand_enzyme_list

        org_cnt = len(enzyme_name_list)
        enzyme_name_list = list(set(enzyme_name_list))
        enzyme_name_list = sorted(enzyme_name_list)
        set_cnt = len(enzyme_name_list)

        for enzyme_name in enzyme_name_list:
            if enzyme_name not in ResDict.rest_dict:
                log.critical("your ENZYME <{}> is not in list.".format(
                    enzyme_name))
                sys.exit(1)

        if org_cnt != set_cnt:
            diff_cnt = org_cnt - set_cnt
            log.info("There were {} duplicates out of {}.".format(
                diff_cnt, org_cnt))

        log.info("{}".format(enzyme_name_list))
        log.info("total enzyme cnt={}.".format(set_cnt))

        return enzyme_files_list, enzyme_name_list


    def _write_new_enzyme_names_file(self, new_enzyme_file_path):

        header = "# name recogseq elucidate seqlen\n\n"

        with open(new_enzyme_file_path, mode='w') as f:
            f.write(header)

            for line in self.default_enzyme_names:
                f.write("{}\n".format("\t".join(line)))


