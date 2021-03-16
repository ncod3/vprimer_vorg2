# -*- coding: utf-8 -*-

# 現在の設定を書き出すのは、ini_fileを書き出すのがいいと思う。
# それを指定すると、同じことができることにすればいい。

import sys
import os
import errno
import re

import pprint

# global variants
import vprimer.glv as glv
import vprimer.utils as utl

from vprimer.logging_config import LogConf
from vprimer.conf_disting import ConfDistinG
from vprimer.conf_vcf_file import ConfVcfFile
from vprimer.conf_curr_setting import ConfCurrSet
from vprimer.conf_enzyme import ConfEnzyme

### class Conf(ConfBase, ConfDistinG, ConfVcfFile):
###     ''' Split class and join with multiple inheritance
###     '''
###     pass

class ConfBase(object):

    def __init__(self):

        self.ini = None
        self.param = None

        # ----------------------------------------------------------
        # A dictionary that associates a value with a variable name 
        # for all parameters
        self.conf_dict = dict()

        # ----------------------------------------------------------
        self.conf_dict = {
            # show_genotype  no / gt / int
            'show_genotype':{'dtype': 'str',   'default': 'no'},

            # debug
            'analyse_caps': {'dtype': 'bool', 'default': 'False'},

            # ini, param, ini, easy
            'ini_version':  {'dtype': 'str',    'default': glv.ini_version},
            'ini_file':     {'dtype': 'str',    'default': ''},
            'out_dir':      {'dtype': 'str',    'default': 'out_vprimer'},
            'thread':       {'dtype': 'int',    'default': '10'},
            'use_joblib_threading': # param
                            {'dtype': 'str',    'default': 'yes'},

            #
            'vcf':          {'dtype': 'str',    'default': ''},
            'ref':          {'dtype': 'str',    'default': ''},

            #
            'pick_mode':    {'dtype': 'str',    'default': 'all'},
            'indel_size':   {'dtype': 'str',    'default': '20-200'},
            'product_size': {'dtype': 'str',    'default': '200-500'},

            # list enzyme_file refs/enzyme_names.txt
            'enzyme_file':  {'dtype': 'str',    'default': 'no_enzyme'},
            # list enzyme
            'enzyme':       {'dtype': 'str',    'default': ''},

            # list target
            'target':       {'dtype': 'str',    'default': ''},
            # list a_samples
            'a_sample':     {'dtype': 'str',    'default': ''},
            # list b_samples
            'b_sample':     {'dtype': 'str',    'default': ''},

            # list regions
            'regions':      {'dtype': 'str',    'default': ''},
            # list distinguish_groups
            'distinguish_groups':
                            {'dtype': 'str',    'default': ''},
            # list group_members
            'group_members':
            # The default for group_member is group_members_vcf_str
            # read from vcf.
            # Do not initialize the key for safety as we will update it later.
            #                {'dtype': 'str',    'default': ''},
                            {'dtype': 'str',},
            # gen p3_params.txt
            'p3_params':    {'dtype': 'str',    'default':
                                                'no_p3_params'},
            'fragment_pad_len':
                            {'dtype': 'int',    'default': '500'},

            'blast_distance':
                            {'dtype': 'int',  'default': '10000'},
            # blast_word_size will be set later after p3 PRIMER_MIN_SIZE
            # is set.
            'blast_word_size': # PRIMER_MIN_SIZE
                            {'dtype': 'int', 'default': '23'},

            #
            'show_samples': {'dtype': 'bool',   'default': 'False'},
            'show_fasta':   {'dtype': 'bool',   'default': 'False'},
            'progress':     {'dtype': 'str',  'default': 'all'},
            'stop':         {'dtype': 'str',  'default': 'none'},

        }

        # cwd, log --------------------------------------------
        self.cwd = glv.cwd
        self.log = LogConf()

        # for debug
        self.analyse_caps = False

        # show_genotype ---------------------------------------
        self.show_genotype = ""         # ---- INI

        # ini file --------------------------------------------
        self.user_ini_file = ""         # ---- INI
        self.ini_file_path = ""
        self.ini_version_user = ""      # ---- INI
        self.ini_version_system = ""    # ---- INI

        # out_dir ---------------------------------------------
        self.user_out_dir = ""          # ---- INI
        self.out_dir_path = ""
        self.log_dir_path = ""
        self.out_bak_dir_path = ""

        # ref_dir ---------------------------------------------
        self.ref_dir_path = ""          # ---- INI

        # out_curr_setting ------------------------------------
        self.curr_setting_file = ""
        self.curr_setting_file_path = ""

        # thread ----------------------------------------------
        self.thread = 0                 # ---- INI
        self.use_joblib_threading = "yes"
        self.parallel = True
        self.parallel_full_thread = 0
        self.parallel_blast_cnt = 0
        self.blast_num_threads = 0

        # vcf -------------------------------------------------
        self.user_vcf_file = ""         # ---- INI
        self.user_vcf_file_path = ""
        self.vcf_file_slink_system = ""
        self.vcf_file_path = ""

        # sample_nickname -------------------------------------
        self.vcf_sample_name_file = ""

        self.vcf_sample_nickname_list = list()
        self.vcf_sample_basename_list = list()
        self.vcf_sample_fullname_list = list()

        self.vcf_sample_nickname_dict = dict()
        self.vcf_sample_basename_dict = dict()
        self.vcf_sample_fullname_dict = dict()

        self.group_members_vcf_str = ""

        # show_samples-----------------------------------------
        self.show_samples = False

        # show_fasta ------------------------------------------
        self.show_fasta = False

        # ref -------------------------------------------------
        self.user_ref_fasta = ""        # ---- INI
        self.user_ref_fasta_path = ""

        # It will be set in main later
        # glv.ref = glv.ref.prepare_ref()
        self.ref_fasta_slink_system = ""
        self.ref_fasta_path = ""
        self.ref_fasta_chrom_dict_list = list()
        self.ref_fasta_chrom_list = list()
        self.ref_fasta_chrom_region_list = list()

        self.ref_fasta_fai = ""
        self.ref_fasta_chrom_txt = ""
        self.ref_fasta_pickle = ""

        # pick_mode -------------------------------------------
        self.pick_mode = ""             # ---- INI

        # indel len -------------------------------------------
        self.indel_size = ""            # ---- INI
        self.min_indel_len = 0
        self.max_indel_len = 0

        # product size ----------------------------------------
        self.product_size = ""          # ---- INI
        self.min_product_size = 0
        self.max_product_size = 0

        # enzyme ----------------------------------------------
        self.enzyme_files_user_str = ""      # ---- INI
        self.enzyme_files_user_list = list()
        self.enzyme_files_list = list()

        self.enzyme_str = ""
        #self.enzyme_list = list()
        self.enzyme_name_list = list()

        # region group member string --------------------------
        # select string by priority, next make dict or list variables
        self.regions_str = ""           # ---- INI
        self.group_members_str = ""     # ---- INI
        self.distinguish_groups_str = ""    # ---- INI

        self.region_name_list = list()
        self.group_name_list = list()

        self.regions_dict = dict()
        self.group_members_dict = dict()
        self.distinguish_groups_list = list()

        # primer3 ---------------------------------------------
        self.p3_params_file = ""        # ---- INI
        self.p3_params_file_path = ""
        self.primer3_header_dict = dict()

        self.fragment_pad_len = 0       # ---- INI

        # primer3 params
        self.p3key = { \
            'PRIMER_MIN_SIZE': 23,
            'PRIMER_OPT_SIZE': 25,
            'PRIMER_MAX_SIZE': 27,
            'PRIMER_MIN_GC': 40,
            'PRIMER_OPT_GC': 50,
            'PRIMER_MAX_GC': 60,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MAX_POLY_X': 4,
            'PRIMER_PAIR_MAX_DIFF_TM': 4,
        }

        # blast -----------------------------------------------
        self.blast_distance = 0         # ---- INI

        # not set now -----------------------------------------
        self.blast_word_size = 0
        self.blastdb_title = ""
        self.blastdb = ""

        # start stop ------------------------------------------
        self.progress = ""
        self.stop = ""


    # --- before log file open
    def collect_param_ini(self, param, ini):
        ''' aggregate all parameters into one dictionary
        '''

        self.param = param
        self.ini = ini

        for vname in self.conf_dict:

            param_value = self._get_param_value(vname, param)
            ini_value = None

            if param.p['ini_file'] is not None:
                 # ini
                ini_value = self._get_ini_value(vname, ini)

            self.conf_dict[vname]['param'] = param_value
            self.conf_dict[vname]['ini'] = ini_value


    def _get_param_value(self, vname, param):
        ''' get data from parameter
        param handles values with the correct data type
        '''

        ret = None

        # a_sample="DRS_013.all.rd DRS_084.all.rd,DRS_099.all.rd ref"
        # b_sample="DRS_025.all.rd, DRS_061.all.rd DRS_101.all.rd"

        # {'a_sample': ['DRS_013.all.rd',
        #   'DRS_084.all.rd,DRS_099.all.rd', 'ref'],
        # 'b_sample': ['DRS_025.all.rd,', 'DRS_061.all.rd',
        #   'DRS_101.all.rd'],

        if vname in param.p:

            val = param.p[vname]

            if type(val) == list:
                # 'DRS_084.all.rd,DRS_099.all.rd'
                # 'DRS_025.all.rd,'
                #print("{}={}".format(vname, val))

                mod_list = list()
                for item in val:
                    sep_list = item.split(",")
                    for sep in sep_list:
                        if sep != "":
                            mod_list.append(sep)

                #ret = ','.join(val)
                ret = ','.join(mod_list)
                #print("\t{}={}".format(vname, ret))

            else:
                ret = val

        return ret


    def _get_ini_value(self, vname, ini):
        ''' get data from ini file
        ini handles values as strings
        '''

        ret = None

        if vname in ini.ini['vprimer']:

            val = ini.ini['vprimer'][vname]
            if type(val) == list:
                ret = ','.join(val)
            else:
                ret = self._cast_val(val, self.conf_dict[vname]['dtype'])

        return ret


    def _cast_val(self, value, dtype):
        ''' for ini file data, casting data to fit data type
        ''' 

        #print("_cast_val: value={}, dtype={}".format(value, dtype))
        #print("type(value={}".format(type(value)))
        #print("type(dtype={}".format(type(dtype)))

        if dtype == 'int':
            #print("int")
            return int(value)

        elif dtype == 'float':
            #print("float")
            return float(value)

        elif dtype == 'bool':

            # for param
            if type(value) == bool:
                return value
            else:
                # form ini
                if value == "True":
                    return True
                elif value == "False":
                    return False
                else:
                    return None

        elif dtype == 'str':
            #print("str")
            return str(value)

        else:
            #print("else")
            return str(value)


    def _make_out_dir_tree(self):
        '''
        '''

        # if already made at conf
        dirs = [
            self.out_dir_path,
            self.log_dir_path,
            self.out_bak_dir_path
        ]

        for dir in dirs:
            if os.path.isdir(dir):
                # prelog
                utl.prelog("exist dir {}.".format(dir), __name__)
            else:
                utl.prelog("not exist dir {}.".format(dir), __name__)
                utl.makedirs(dir)


    def out_dir_logging_start(self):
        '''
        '''

        # user defined path
        self.user_out_dir = self._value_choice('out_dir')

        # absolute path
        self.out_dir_path = utl.full_path(self.user_out_dir)
        # log_dir
        self.log_dir_path = "{}/{}".format(self.out_dir_path, "logs")
        # bak_dir
        self.out_bak_dir_path = "{}/{}".format(self.out_dir_path, "bak")

        # out_dir and log_dir
        self._make_out_dir_tree()

        # for conf
        global log
        log = self.log.logging_start(
            __name__, self.out_dir_path, self.log_dir_path)

       # for utl
        utl.open_log()


    # --- after log file open
    def choice_variables(self):
        ''' decide variable values
        '''

        # print param and ini variables
        self._print_param_ini()

        # for debug
        self.analyse_caps = self._value_choice('analyse_caps')

        # out_dir ---------------------------------------------
        self.user_out_dir = self._value_choice('out_dir')
        self.out_dir_path = utl.full_path(self.user_out_dir)

        # vcf -------------------------------------------------
        self.user_vcf_file = self._value_choice('vcf')
        self.user_vcf_file_path = utl.full_path(self.user_vcf_file)

        # ref -------------------------------------------------
        self.user_ref_fasta = self._value_choice('ref')
        self.user_ref_fasta_path = utl.full_path(self.user_ref_fasta)

        # thread ----------------------------------------------
        self.thread = self._value_choice('thread')


        if self.out_dir_path == "" or \
            self.user_vcf_file_path == "" or \
            self.user_ref_fasta_path == "":
            err_mes = "out_dir={} and vcf={} and ref={} ".format(
                self.out_dir_path,
                self.user_vcf_file_path,
                self.user_ref_fasta_path)
            err_mes += "are all required. exit."
            log.error(err_mes)
            sys.exit(1)

        log.info("thread={}".format(self.thread))


        # out_dir ---------------------------------------------
        self.out_dir_path = utl.full_path(self.user_out_dir)
        self.log_dir_path = "{}/{}".format(self.out_dir_path, "logs")
        self.out_bak_dir_path = "{}/{}".format(self.out_dir_path, "bak")

        #pprint.pprint(self.conf_dict)
        # INI show_genotype
        self.show_genotype = self._value_choice('show_genotype')
        if self.show_genotype == "":
            self.show_genotype = "gt"

        if not self.show_genotype in glv.show_genotype_list:
            err_mes = "show_genotype is selected from one of "
            err_mes += ", ".join(glv.show_genotype_list)
            log.error("{}. exit.".format(err_mes))
            log.error("show_genotype={}".format(self.show_genotype))
            sys.exit(1)


        # ini_file --------------------------------------------
        # INI
        self.ini_version_user = self._value_choice('ini_version')
        self.ini_version_system = self.conf_dict['ini_version']['default']

        self.user_ini_file = self.conf_dict['ini_file']['param']
        self.ini_file_path = utl.full_path(self.user_ini_file)

        # ref_dir ---------------------------------------------
        self.ref_dir_path = utl.full_path("refs")
        # make ref_dir
        utl.makedirs(self.ref_dir_path)

        # out_curr_setting ------------------------------------
        self.curr_setting_file = "current_setting_ini.txt"
        self.curr_setting_file_path = "{}/{}".format(
            self.out_dir_path, self.curr_setting_file)

        # thread ----------------------------------------------
        self.use_joblib_threading = self._value_choice('use_joblib_threading')

        if not self.use_joblib_threading in ['yes', 'no']:
            err_mes = "use_joblib_threading Choose from Yes or No."
            log.error("{} exit.".format(err_mes))
            log.error("use_joblib_threading={}".format(
                self.use_joblib_threading))
            sys.exit(1)

        # thread adjust
        self.parallel, \
        self.parallel_full_thread, \
        self.parallel_blast_cnt, \
        self.blast_num_threads \
            = self._thread_adjusting()

        # vcf -------------------------------------------------
        basename_user_vcf = os.path.basename(self.user_vcf_file_path)
        self.vcf_file_slink_system = "{}/{}{}".format(
            self.ref_dir_path, 'slink_', basename_user_vcf)

        # gtonly.gz
        self.vcf_file_path = "{}/{}{}".format(
            self.ref_dir_path, basename_user_vcf, "_GTonly.vcf.gz")

        # read
        self.prepare_vcf()

        # sample_nickname -------------------------------------
        basename_vcf_file = os.path.basename(self.vcf_file_path)
        self.vcf_sample_name_file = "{}/sample_name_{}.txt".format(
            self.ref_dir_path, basename_vcf_file)
        self.save_vcf_sample_name_txt()

        self.vcf_sample_nickname_list, \
        self.vcf_sample_basename_list, \
        self.vcf_sample_fullname_list, \
        self.vcf_sample_nickname_dict, \
        self.vcf_sample_basename_dict, \
        self.vcf_sample_fullname_dict, \
        self.group_members_vcf_str \
            = self.make_vcf_sample_variable()

        # illegal
        self.conf_dict['group_members']['default'] = \
            self.group_members_vcf_str

        # show_fasta-------------------------------------------
        self.show_fasta = self._value_choice('show_fasta')

        # show_samples-----------------------------------------
        self.show_samples = self._value_choice('show_samples')

        # Because it stops at show_fasta
        if self.show_fasta != True and self.show_samples == True:
            log.info("only show_samples mode, exit.")
            log.info("program finished {}\n".format(
                utl.elapsed_time(time.time(), glv.now_epochtime)))
            sys.exit(1)


        # pick_mode
        self.pick_mode = self._value_choice('pick_mode')

        # indel len
        self.indel_size = self._value_choice('indel_size')
        self.min_indel_len, self.max_indel_len = \
            [int(i) for i in self.indel_size.split('-')]

        # product size
        self.product_size = self._value_choice('product_size')
        self.min_product_size, self.max_product_size = \
            [int(i) for i in self.product_size.split('-')]

        # ref -------------------------------------------------
        # It will be set in main later
        # glv.ref = glv.ref.prepare_ref()
        self.ref_fasta_slink_system = ""
        self.ref_fasta_path = ""
        self.ref_fasta_chrom_list = []
        self.ref_fasta_fai = ""
        self.ref_fasta_chrom_txt = ""
        self.ref_fasta_pickle = ""

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # enzyme
        self.enzyme_files_user_str = self._value_choice('enzyme_file')
        #self.enzyme_files_user_list = list()
        #self.enzyme_files_list = list()
        self.enzyme_str = self._value_choice('enzyme')
        #self.enzyme_name_list = list()

        # start stop
        self.progress = self._value_choice('progress')
        self.stop = self._value_choice('stop')

        # primer3
        self.fragment_pad_len = self._value_choice('fragment_pad_len')
        self.p3_params_file = self._value_choice('p3_params')

        # blast
        self.blast_distance = self._value_choice('blast_distance')

        # not set now
        self.blast_word_size = 0
        self.blastdb_title = ""
        self.blastdb = ""

        # region group member string ---------------------------------
        # select string by priority, next make dict or list variables

        self.regions_str = self.set_regions_str()
        self.group_members_str = self.set_group_members_str()
        self.distinguish_groups_str = self.set_distinguish_groups_str()


    def setup_variables(self):
        '''
        '''

        # setup only regions and members
        if self.show_genotype != "no":
            self._setup_genotype_variables()

        # Satisfy three structural variables
        #   1) regions_dict
        #   2) distinguish_groups_list
        #   3) group_members_dict

        # 1.1) Adjusting regions_str in easy mode --------------------------
        #log.debug("org self.regions_str={}".format(self.regions_str))

        #   Easy mode lacks the region name, so make up for it here.
        if '<EASY_MODE>' in self.regions_str:
            # <EASY_MODE>chrom_01:1-200000,chrom_02,chrom_03:all
            regions_str = ""
            rg_cnt = 1
            for region in self.regions_str.split(','):
                region = re.sub(r"^<EASY_MODE>", "", region)
                # region_name = region1
                region_name = "easy_region{}".format(rg_cnt)
                regions_str += "{}:{},".format(region_name, region)
                rg_cnt += 1

            # complete
            self.regions_str = re.sub(r",$", "", regions_str)

        # ====================
        # 1.2) make regions dictionary from regions_str ....................
        self.regions_str, \
        self.regions_dict, \
        self.region_name_list \
            = self.set_regions_dict(self.regions_str)

        self._set_chosen_value('regions', self.regions_str)
        log.info("chosen self.regions_str={}".format(
            self.regions_str))

        #log.debug("glv.conf.regions_dict={}".format(self.regions_dict))
        #log.debug("glv.conf.region_name_list={}".format(
        #    self.region_name_list))

        # 2.1) group_members_str -------------------------------------------
        #log.debug("org self.group_members_str={}".format(
        #    self.group_members_str))

        # 2.2) make group_members dictionary from group_members_str ........
        self.group_members_dict, \
        self.group_name_list \
            = self.set_group_members_dict(self.group_members_str)

        self._set_chosen_value('group_members',
            self.group_members_str)
        log.info("chosen self.group_members_str={}".format(
            self.group_members_str))

        #log.debug("glv.conf.group_members_dict={}".format(
        #    self.group_members_dict))
        #log.debug("glv.conf.group_name_list={}".format(
        #    self.group_name_list))

        #print("self.group_members_str={}".format(self.group_members_str))
        #print("self.group_members_dict={}".format(self.group_members_dict))
        #print("self.group_name_list={}".format(self.group_name_list))


        # 3.1) distinguish_groups_str --------------------------------------
        #log.debug("org self.distinguish_groups_str={}".format(
        #    self.distinguish_groups_str))

        #   Easy mode lacks the region name, so make up for it here.
        if "<EASY_MODE>" in self.distinguish_groups_str:
            region_names_str = "+".join(self.region_name_list)
            self.distinguish_groups_str = re.sub(
                r"<EASY_MODE>", region_names_str,
                self.distinguish_groups_str)

        # 3.2) make distinguish_groups_list from distinguish_groups_str ....
        # Avoid checking when show_genotype
        if self.show_genotype == "no":
            self.distinguish_groups_str, \
            self.distinguish_groups_list \
                = self.set_distinguish_groups_list(
                    self.distinguish_groups_str)

        self._set_chosen_value('distinguish_groups',
            self.distinguish_groups_str)
        log.info("chosen self.distinguish_groups_str={}".format(
            self.distinguish_groups_str))

        #log.debug("distinguish_groups_list=\n{}".format(
        #    pprint.pformat(self.distinguish_groups_list)))


        # 4.1) Create a file to set the environment of primer3
        if self.p3_params_file == "no_p3_params":
            self.p3_params_file = "{}/{}".format(
                self.ref_dir_path, "p3_params.txt")

        self.p3_params_file_path = utl.full_path(self.p3_params_file)
        #log.debug("self.p3_params_file_path={}".format(
        #    self.p3_params_file_path))

        # primer3, make and read parameters
        self.primer3_header_dict = self._set_primer3_header_dict()

        # 4.2) blast_word_size = PRIMER_MAX_SIZE
        self.blast_word_size = self.primer3_header_dict['PRIMER_MIN_SIZE']

        # 5.1) File that describes the enzyme name to be handled
        # enzyme_files_list
        self.enzyme_files_user_list = list()
        for file in self.enzyme_files_user_str.split(','):
            # user file full path list
            file_path_user = utl.full_path(file)
            self.enzyme_files_user_list.append(file_path_user)


            basename_user = os.path.basename(file_path_user)
            enzyme_file_slink_system = "{}/{}{}".format(
                self.ref_dir_path, "slink_", basename_user)

            self.enzyme_files_list.append(enzyme_file_slink_system)

        #pprint.pprint(self.enzyme_files_user_list)
        #pprint.pprint(self.enzyme_files_list)

        # 5.2) enzyme
        self.enzyme_files_list, self.enzyme_name_list \
            = self.read_enzyme_file()

        log.info("enzyme_files_list={}".format(self.enzyme_files_list))
        log.info("enzyme_name_list={}".format(self.enzyme_name_list))

        # 6) progress, stop
        if not self._check_progress_stop():
            if self.progress != "all" and self.stop != "no":
                err_mes = "The progress or stop settings are incorrect."
                log.error("{} exit.".format(err_mes))
                log.error("progress={}, stop={}".format(
                    self.progress, self.stop))
                log.error("{}".format(
                    ", ".join(glv.outlist.outf_prefix.keys())))
                sys.exit(1)


    def _setup_genotype_variables(self):
        '''
        '''

        if self.regions_str == "all":
            self.regions_str = "{}".format(
                ",".join(self.ref_fasta_chrom_list))

        if self.group_members_str == "all":
            self.group_members_str = "all:{}".format(
                ",".join(self.vcf_sample_nickname_list))

        '''
        print("regions_str={}".format(
            self.regions_str))
        print("group_members_str={}".format(
            self.group_members_str))
        print("distinguish_groups_str={}".format(
            self.distinguish_groups_str))


        print("region_name_list={}".format(
            self.region_name_list))
        print("group_name_list={}".format(
            self.group_name_list))

        print("regions_dict={}".format(
            self.regions_dict))
        print("group_members_dict={}".format(
            self.group_members_dict))
        print("distinguish_groups_list={}".format(
            self.distinguish_groups_list))
        #sys.exit(1)
        '''

        
    def _check_progress_stop(self):
        '''
        '''

        ret = True

        if self.progress == "all" and self.stop == "no":
            ret = True

        elif self.progress not in glv.outlist.outf_prefix.keys():
            ret = False

        elif self.stop not in glv.outlist.outf_prefix.keys():
            ret = False

        return ret


    def _thread_adjusting(self):
        ''' in Parallel, if there are 10 threads blast cmd will use at least
            2 cores so par 1 parallel.
            main python always use 1,
            parallel use 1 thread, blast use 2 threads
        '''

        parallel = True

        parallel_full_thread = 0
        parallel_blast_cnt = 0
        blast_num_threads = 0

        if self.thread < 6 or self.use_joblib_threading != "yes":
            parallel = False

        if parallel == True:

            # unit is 4+1=5
            parallel_base = self.thread
            parallel_full_thread = parallel_base

            # blast = 4
            parallel_blast_cnt = int(parallel_base / 5)
            blast_num_threads = 4

        else:
            # 6 = 5
            full_thread = self.thread - 1
            parallel_full_thread = full_thread
            blast_num_threads = full_thread

        return parallel, parallel_full_thread,\
            parallel_blast_cnt, blast_num_threads


    def _value_choice(self, vname):
        '''
        '''

        chosen_value = ""

        # If neither param nor ini has a key, use the default value
        if self.conf_dict[vname]['param'] is None and \
            self.conf_dict[vname]['ini'] is None:
            #print("cv 1")
            chosen_value = self.conf_dict[vname]['default']

        elif self.conf_dict[vname]['param'] is not None:
            #print("cv 2")
            chosen_value = self.conf_dict[vname]['param']

        elif self.conf_dict[vname]['ini'] is not None:
            #print("cv 3")
            chosen_value = self.conf_dict[vname]['ini']

        elif self.conf_dict[vname]['default'] is not None:
            #print("cv 4")
            chosen_value = self.conf_dict[vname]['default']

        else:
            utl.prelog("not found value of key {}.".format(vname), __name__)
            sys.exit(1)

        #print("vname={}".format(vname))
        #print("dtype={}".format(self.conf_dict[vname]['dtype']))
        #print(chosen_value)
        #print(type(chosen_value))

        # cast by dtype
        chosen_value = self._cast_val(
            chosen_value, self.conf_dict[vname]['dtype'])

        #print(type(chosen_value))

        # chosen value
        self._set_chosen_value(vname, chosen_value)
        #print("{}={}".format(vname, self.conf_dict[vname]['chosen']))

        return chosen_value


    def _set_chosen_value(self, vname, chosen_value):

        # Assuming the chosen key does not exist
        self.conf_dict[vname]['chosen'] = chosen_value


    def _is_chrom_name(self, chrom_name):
        '''
        '''
        ret = False
        #print("_is_chrom_name, {}, {}".format(
        #    chrom_name, self.ref_fasta_chrom_list))
        if chrom_name in self.ref_fasta_chrom_list:
            ret = True
        #print("{}".format(ret))
        return ret


    def _is_region_name(self, region_name):
        '''
        '''
        ret = False
        #print("_is_region_name, {}, {}".format(
        #    region_name, self.region_name_list))
        if region_name in self.region_name_list:
            ret = True
        #print("{}".format(ret))
        return ret


    def _is_group_name(self, group_name):
        '''
        '''
        ret = False
        #print("_is_group_name, {}, {}".format(
        #    group_name, self.group_name_list))
        if group_name in self.group_name_list:
            ret = True
        #print("{}".format(ret))
        return ret


    def _is_valid_int_range(self, range_str):
        '''
        '''
        ret = True

        if "-" not in range_str:
            ret = False
        else:
            min_size, max_size = range_str.split("-")

            if not min_size.isdecimal() or \
                not max_size.isdecimal():
                ret = False
            elif int(min_size) > int(max_size):
                ret = False

        return ret


    def _is_valid_chrom_range(self, range_str):

        ret = True

        #print("_is_valid_chrom_range, range_str={}".format(range_str))
        chrom_name, rg_str = range_str.split(':')
        #print("_is_valid_chrom_range, chrom_name={}, rg_str={}".format(
        #    chrom_name, rg_str))

        min_pos, max_pos = [int(i) for i in rg_str.split('-')]
        #print("_is_valid_chrom_range, min_pos={}, max_pos={}".format(
        #    min_pos, max_pos))

        region_str, start, end, length = self._get_chrom_info(chrom_name)
        #print("_is_valid_chrom_range, _get_chrom_info={}, {}, {}, {}".format(
        #    region_str, start, end, length))

        if min_pos < start or end < max_pos:
            ret = False
        #print("_is_valid_chrom_range, ret={}".format(ret))

        return ret


    def _get_chrom_info(self, chrom_name):
        '''
        '''

        # glv.conf.ref_fasta_chrom_dict_list
        start = 0
        end = 0
        length = 0

        for d in self.ref_fasta_chrom_dict_list:
            #end': 30583384, 'length': 30583384, 'start': 1
            if chrom_name == d.get('chrom'):
                start = d.get('start')
                end = d.get('end')
                length = d.get('length')

        region_str = chrom_name

        if start is not None:
            region_str = "{}:{}-{}".format(region_str, start, end)

        #pprint.pprint(glv.conf.ref_fasta_chrom_dict_list)
        #print("{}, {}".format(region_str, length))

        return region_str, start, end, length


    def _is_easy_mode(self):
        ''' On the command line, determine if we are currently in easy mode
        '''

        easy_mode = False
        a_sample = False
        b_sample = False

        #print("self.conf_dict['a_sample']['param']={}".format(
        #    self.conf_dict['a_sample']['param']))
        #print("self.conf_dict['b_sample']['param']={}".format(
        #    self.conf_dict['b_sample']['param']))

        if self.conf_dict['a_sample']['param'] is not None and \
            self.conf_dict['a_sample']['param'] != "":
            a_sample = True

        if self.conf_dict['b_sample']['param'] is not None and \
            self.conf_dict['b_sample']['param'] != "":
            b_sample = True

        #print("a_sample={}".format(a_sample))
        #print("b_sample={}".format(b_sample))

        if a_sample == True or b_sample == True:
            easy_mode = True

        #print("easy_mode={}".format(easy_mode))

        return easy_mode
    

    def _print_param_ini(self):
        ''' for debug
        '''
        # parameter and ini_file variables

        #utl.prelog
        log.info("\n======== param.p ====================")
        log.info("self.param.p=\n\n{}\n".format(pprint.pformat(self.param.p)))

        log.info("\n======== ini.ini['vprimer'] =========")

        if self.param.p['ini_file'] is not None:
            log.info("\nini_file={}\n\n{}\n".format(
                self.param.p['ini_file'],
                pprint.pformat(dict(self.ini.ini['vprimer']))))
        else:
            log.info("ini_file not specified.")
        log.info("\n=====================================\n")


    def _set_primer3_header_dict(self):
        '''
        '''

        primer3_header_dict = dict()

        if os.path.isfile(self.p3_params_file_path):
            log.info("found {}.".format(self.p3_params_file_path))
            # This file may have been edited by the user, so copy it 
            utl.save_to_tmpfile(self.p3_params_file_path, True, True)

        else:
            log.info("not found {}.".format(self.p3_params_file_path))
            with open(self.p3_params_file_path, mode='w') as f:
                f.write("{}={}\n".format('#PARAM', 'VALUE'))

                for key, value in list(self.p3key.items()):
                    f.write("{}={}\n".format(key, value))

        # 1.1) open and read parameters
        with open(self.p3_params_file_path, mode='r') as f:
            # iterator
            for r_liner in f:
                r_line = r_liner.strip()    # cr, ws

                if r_line.startswith('#') or r_line == '':
                    continue

                r_line = utl.strip_hash_comment(r_line)
                vname, value = r_line.split('=')
                if vname == 'PRIMER_PRODUCT_SIZE_RANGE' or \
                    vname == 'PRIMER_NUM_RETURN':
                    continue

                primer3_header_dict[vname] = value


        # constant value for primer3

        # PRIMER_FIRST_BASE_INDEX=1
        primer3_header_dict['PRIMER_FIRST_BASE_INDEX'] = str(1)
        # PRIMER_PRODUCT_SIZE_RANGE=???-???
        primer3_header_dict['PRIMER_PRODUCT_SIZE_RANGE'] = \
            "{}-{}".format(self.min_product_size, self.max_product_size)
        # PRIMER_NUM_RETURN=1
        primer3_header_dict['PRIMER_NUM_RETURN'] = str(1)

        return primer3_header_dict


class Conf(ConfBase, ConfDistinG, ConfVcfFile, ConfCurrSet, ConfEnzyme):
    ''' Split class and join with multiple inheritance
    '''
    pass


