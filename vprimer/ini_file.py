# -*- coding: utf-8 -*-

import sys
import os
import errno
import time

import re
import pprint

# global configuration
import vprimer.glv as glv
import vprimer.utils as utl

import configparser

class IniFile(object):

    def __init__(self):

        self.user_ini_file = ''
        self.ini_file_path = ''

        # for default
        self.ini = dict()

    def read_ini_file(self, user_ini_file):
        '''
        '''

        # don't permit empty_lines_in_values
        self.ini = configparser.ConfigParser(
            empty_lines_in_values=False
        )
        # don't convert to lower case
        self.ini.optionxform = str

        # ini_file full_path
        self.user_ini_file = user_ini_file
        self.ini_file_path = utl.full_path(user_ini_file)

        utl.prelog(
             "ini_file_path = {}".format(self.ini_file_path), __name__)
 
        if os.path.exists(self.ini_file_path):
            with open(self.ini_file_path, encoding='utf-8') as fp:
                self.ini.read_file(fp)
                #  adjustment of variable format
                self.ini = self._format_ini_variable(self.ini)
        else:
             utl.prelog(
                 "not found {}. exit.".format(self.ini_file_path), __name__)
             sys.exit(1)


    def _format_ini_variable(self, ini):

        for section in ini.sections():
            for key in ini[section]:

                val = ini[section][key]
                # remove hash comment
                val = utl.strip_hash_comment(val)
                # remove \n at the beginning of value, not necessary \n
                val = val.lstrip()
                # replace internal \n to semicolons
                val = val.replace('\n', ',')

                # replace white space to one space
                if key == 'group_members':
                    val = re.sub(r"\s+", " ", val)
                    val = re.sub(r"\s*:\s*", ":", val)
                    val = re.sub(r" ", ",", val)
                    val = re.sub(r",+", ",", val)
                    val = re.sub(r",;", ";", val)
                    val = re.sub(r";", ",", val)

                else:
                    # remove white space
                    val = re.sub(r"\s+", "", val)

                # reset
                ini[section][key] = val

        return ini

