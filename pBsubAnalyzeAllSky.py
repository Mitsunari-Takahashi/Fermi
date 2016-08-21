#!/usr/bin/env python
"""Throw bjobs for all sky analysis.
Execute in a directory which has a file of hBDT_CutValues_{bdtname}.root.
"""
import sys
import commands
import click
import subprocess

@click.command()
@click.argument('pathfiles', type=str)
@click.option('--bdtname', '-b', default='S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15_catZDIR060_BDTG500D06_catZDIR060')
def main(pathfiles, bdtname):
    cmd_ls = "ls {0}".format(pathfiles)
    ret_ls = commands.getoutput(cmd_ls)
    li_path_files = ret_ls.split("\n")
    str_date = commands.getoutput("\date +%Y%m%d%H%M")
    for path_file in li_path_files:
        print '----------'
        print path_file
        str_year = path_file[-17:-13]
        str_week = path_file[-8:-5]
        li_cmd = ['bsub', '-ologs/{0}anaLPA{1}_{2}_{3}.log'.format(str_date, str_year, str_week, bdtname), '-Jana{0}_{1}'.format(str_year, str_week), '-W200', 'python', '/u/gl/mtakahas/work/PythonModuleMine/Fermi/pAnalyzeAllSky.py', bdtname, '{0}_{1}'.format(str_year, str_week), path_file]
        print li_cmd
        subprocess.call(li_cmd)


if __name__ == '__main__':
    main()
