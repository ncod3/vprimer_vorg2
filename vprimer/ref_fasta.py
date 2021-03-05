# -*- coding: utf-8 -*-

import sys
import os
import errno
import time
import datetime

# global configuration
import vprimer.glv as glv

import subprocess as sbp
import pandas as pd
import pickle
import pprint

import vprimer.utils as utl

from vprimer.blast import Blast
from vprimer.logging_config import LogConf


class RefFasta(object):

    def __init__(self):

        # as glv.ref.refseq
        self.refseq = dict()

        # as glv.ref.ref_fasta_pickle
        self.ref_fasta_pickle = ''

    def open_log(self):

        global log
        log = LogConf.open_log(__name__)


    def prepare_ref(self):

        #AJLKJHLKJHLKAJSDLKABSDLKJAHSLKDJ
        #AJLKJHLKJHLKAJSDLKABSDLKJAHSLKDJ
        #AJLKJHLKJHLKAJSDLKABSDLKJAHSLKDJ

        # user_ref_fasta_path: existence confirmation
        if os.path.isfile(glv.conf.user_ref_fasta_path):
            log.info("found {}.".format(glv.conf.user_ref_fasta_path))
        else:
            log.info("not found {}. exit.".format(
                glv.conf.user_ref_fasta_path))
            sys.exit(1)

        # ext, basename, without_ext
        # https://note.nkmk.me/python-os-basename-dirname-split-splitext/
        basename_user = os.path.basename(glv.conf.user_ref_fasta_path)
        root_ext_pair = os.path.splitext(glv.conf.user_ref_fasta_path)
        without_ext = root_ext_pair[0]
        basename_without_ext = os.path.basename(without_ext)
        ext = root_ext_pair[1]

        # ref_fasta_slink_system
        # symbolic link user's fasta to sys_ref_dir as .org(.gz)
        if ext == '.gz':
            glv.conf.ref_fasta_slink_system = "{}/{}{}".format(
                glv.conf.ref_dir_path, "slink_", basename_user)
            # for blast
            glv.conf.blastdb_title = basename_without_ext
            
        else:
            glv.conf.ref_fasta_slink_system = "{}/{}{}".format(
                glv.conf.ref_dir_path, "slink_", basename_user)
            # for blast
            glv.conf.blastdb_title = basename_user

        glv.conf.blastdb = "{}/{}{}".format(
            glv.conf.ref_dir_path, glv.conf.blastdb_title, '.blastdb')


        log.info("glv.conf.blastdb={}".format(glv.conf.blastdb))

        if os.path.isfile(glv.conf.ref_fasta_slink_system):
            log.info("found {}.".format(glv.conf.ref_fasta_slink_system))
        else:
            utl.ln_s(
                glv.conf.user_ref_fasta_path, glv.conf.ref_fasta_slink_system)

        log.info("ext ({}).".format(ext))

        # convert to bgz if ext is .gz and set to ref_fasta
        if ext == '.gz':

            # it should be convert to bgz in ref_dir_path
            glv.conf.ref_fasta_path = "{}/{}".format(
                glv.conf.ref_dir_path, basename_user)

            log.info("ext {}, glv.conf.ref_fasta_path={}.".format(
                ext, glv.conf.ref_fasta_path))

            # half of thread?
            cmd1 = 'bgzip -cd -@ {} {} | bgzip -@ {} > {}'.format(
                        glv.conf.parallel_full_thread,
                        glv.conf.ref_fasta_slink_system,
                        glv.conf.parallel_full_thread,
                        glv.conf.ref_fasta_path)

        else:
            # it should be convert to bgz in ref_dir_path
            glv.conf.ref_fasta_path = "{}/{}{}".format(
                glv.conf.ref_dir_path,
                basename_user,
                '.gz')

            cmd1 = 'bgzip -c -@ {} {} > {}'.format(
                        glv.conf.parallel_full_thread,
                        glv.conf.ref_fasta_slink_system,
                        glv.conf.ref_fasta_path)

        # execute
        if os.path.isfile(glv.conf.ref_fasta_path):
            log.info("found {}.".format(glv.conf.ref_fasta_path))

        else:
            log.info("not found {}. do cmd={}".format(
                glv.conf.ref_fasta_path, cmd1))

            utl.try_exec(cmd1)
 
        # make fai file
        cmd2 = 'samtools faidx {}'.format(
            glv.conf.ref_fasta_path, glv.conf.log_dir_path)

        glv.conf.ref_fasta_fai = \
            "{}{}".format(glv.conf.ref_fasta_path, '.fai')
        glv.conf.ref_fasta_chrom_txt = \
            "{}{}".format(glv.conf.ref_fasta_path, '.chrom.txt')

        if os.path.isfile(glv.conf.ref_fasta_fai):
            log.info("found {}.".format(glv.conf.ref_fasta_fai))

        else:
            utl.try_exec(cmd2)

        # read fasta to dict vprimer.cnf.refseq
        glv.conf.ref_fasta_pickle = "{}{}".format(
            glv.conf.ref_fasta_path, '.pickle')
        self._read_fasta()

        # ref to makeblastdb
        Blast.makeblastdb()

        return self

    def pick_refseq(self, chrom, start_coordinate, end_coordinate):
        ''' for refseq substr, etc...
        '''

        slice_stt = start_coordinate - 1
        slice_end = end_coordinate

        #   1   2   3   4   5   6   coordinate   3-5 tho
        # +---+---+---+---+---+---+
        # | 0 | 1 | 2 | 3 | 4 | 5 | idx
        # | P | y | t | h | o | n |
        # +---+---+---+---+---+---+
        # 0   1   2   3   4   5   6 slice        2-5
        #

        return self.refseq[chrom][slice_stt:slice_end]


    def _read_fasta(self):
        '''
        '''
        if os.path.isfile(glv.conf.ref_fasta_pickle):
            log.info("found {}.".format(glv.conf.ref_fasta_pickle))
            self._read_fasta_pickle()
        else:
            log.info("not found {}.".format(glv.conf.ref_fasta_pickle))
            self._read_fasta_first()


    def _read_fasta_pickle(self):
        '''
        '''
        with open(glv.conf.ref_fasta_pickle, 'rb') as f:
            glv.ref.refseq = pickle.load(f)

        # dictionary
        glv.conf.ref_fasta_chrom_dict_list, \
        glv.conf.ref_fasta_chrom_list, \
        glv.conf.ref_fasta_chrom_region_list = \
            self._get_fai_info()


    def _get_fai_info(self):
        '''
        '''
# glv.conf.ref_fasta_fai
#chr01   43270923    7   60  61
#chr02   35937250    43992120    60  61
#chr03   36413819    80528332    60  61
#chr04   35502694    117549055   60  61
#chr05   29958434    153643468   60  61
#chr06   31248787    184101217   60  61
#chr07   29697621    215870825   60  61
#chr08   28443022    246063414   60  61
#chr09   23012720    274980494   60  61
#chr10   23207287    298376767   60  61
#chr11   29021106    321970850   60  61
#chr12   27531856    351475649   60  61

        # get chrom list from fai text
        df_fai = pd.read_csv(
            glv.conf.ref_fasta_fai, sep = '\t',
            header = None, index_col = None)

        ref_fasta_chrom_dict_list = list()
        ref_fasta_chrom_list = list()
        ref_fasta_chrom_region_list = list()

        for row in df_fai.itertuples():
            chrom_dict = dict()
            chrom_dict = {
                'chrom': row[1],
                'start': 1,
                'end': row[2],
                'length': row[2]
            }
            ref_fasta_chrom_dict_list.append(chrom_dict)
            ref_fasta_chrom_list.append(row[1])
            region = "{}:{}-{}".format(row[1], 1, row[2])
            ref_fasta_chrom_region_list.append(region)


        log.info("ref_fasta_chrom_dict_list={}.".format(
            ref_fasta_chrom_dict_list))

        log.info("ref_fasta_chrom_list={}.".format(
            ref_fasta_chrom_list))

        log.info("fai {}, chrom cnt={}".format(
            glv.conf.ref_fasta_fai, len(ref_fasta_chrom_list)))

        pd_json = pd.json_normalize(ref_fasta_chrom_dict_list)
        log.info("\n{}\n".format(pd_json))

        # exist or not, vcf_sample_name_file
        if os.path.isfile(glv.conf.ref_fasta_chrom_txt):
            log.info("found. {}".format(glv.conf.ref_fasta_chrom_txt))
        else:
            # write to vcf_sample_name_file
            with open(glv.conf.ref_fasta_chrom_txt, mode='w') as f:
                f.write("{}\n".format(pd_json))
            log.info("save. {}".format(glv.conf.ref_fasta_chrom_txt))

        '''
        Creating dataframe by converting dict to list of items
        '''
        if glv.conf.show_fasta == True:
            log.info("only show_fasta mode, exit.")
            log.info("program finished {}\n".format(
                utl.elapsed_time(time.time(), glv.now_epochtime)))
            sys.exit(1)

        return ref_fasta_chrom_dict_list, \
            ref_fasta_chrom_list, \
            ref_fasta_chrom_region_list
            

    def _read_fasta_first(self):
        '''
        '''

        # read fai and set dictionary
        glv.conf.ref_fasta_chrom_dict_list, \
        glv.conf.ref_fasta_chrom_list, \
        glv.conf.ref_fasta_chrom_region_list = \
            self._get_fai_info()

        # for each chrom name
        log.info("read refseq by samtools faidx from fasta {}".format(
            glv.conf.ref_fasta_path))
        start = time.time()

        chrom_seq_list = []
        last_chrom = ''
        for chrom in glv.conf.ref_fasta_chrom_list:

            # get sequence from samtools command
            cmd1 = "samtools faidx {} {}".format(
                glv.conf.ref_fasta_path, chrom)
            cmd_list = cmd1.split(' ')

            # log.info("{}".format(cmd_list))

            # using command output by pipe, get sequence into python
            proc = sbp.Popen(
                cmd_list, stdout = sbp.PIPE, stderr = sbp.PIPE)

            # got bytes (b'')
            for byte_line in proc.stdout:
                # bytes to str, strip \n
                b_line = byte_line.decode().strip()

                #print("{}={}(top)".format(chrom, len(chrom_seq_list)))
                # fasta header
                if b_line.startswith('>'):
                    # not the first time
                    if len(chrom_seq_list) != 0:
                        # dictionary
                        glv.ref.refseq[last_chrom] = ''.join(chrom_seq_list)
                        chrom_seq_list = []
                        continue

                else:
                    # append to list
                    chrom_seq_list.append(b_line)
                    #print("{}={}(append)".format(chrom, len(chrom_seq_list)))

            last_chrom = chrom

        if len(chrom_seq_list) != 0:
            glv.ref.refseq[last_chrom] = ''.join(chrom_seq_list)
            #print("{},last_len={}".format(
            #    chrom, len(glv.ref.refseq[last_chrom])))

        log.info("read refseq done {}\n".format(
            utl.elapsed_time(time.time(), start)))

        # pickle
        with open(glv.conf.ref_fasta_pickle, 'wb') as f:
            pickle.dump(glv.ref.refseq, f)

        log.info('dumped glv.ref.refseq->{}'.format(
            glv.conf.ref_fasta_pickle))


