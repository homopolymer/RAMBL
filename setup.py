#!/usr/bin/env python

from setuptools import setup
from setuptools.command.install import install
from distutils.command.build import build
from subprocess import call
import glob,os

BASEPATH = os.path.dirname(os.path.abspath(__file__))
STRAINCALL_PATH = os.path.join(BASEPATH,'StrainCall')
STRAINCALL_CPP = glob.glob(os.path.join(STRAINCALL_PATH,'*.cpp'))
BLASTOUT2ABUNDANCE_PATH = os.path.join(BASEPATH,'scripts')
BLASTOUT2ABUNDANCE_CPP = [os.path.join(BASEPATH,'scripts/blastout2abundance.cpp')]
EXTRACTMAPREAD_PATH = os.path.join(BASEPATH,'scripts')
EXTRACTMAPREAD_CPP=[os.path.join(BASEPATH,'scripts/extract_mapped_reads.cpp')]

class ExBuild(build):
    def run(self):
        # run the original build code
        build.run(self)
        build_path = os.path.abspath(self.build_temp)
        
        # build StrainCall
        cmd = [os.environ['CXX'],'-std=c++11','-o',os.path.join(STRAINCALL_PATH,'StrainCall')]
        cmd.extend(STRAINCALL_CPP)

        def compile_straincall():
            call(cmd)

        self.execute(compile_straincall, [], 'Compiling StrainCall')

        # build blastout2abundance
        cmd = [os.environ['CXX'],'-std=c++11','-o',os.path.join(BLASTOUT2ABUNDANCE_PATH,'blastout2abundance')]
        cmd.extend(BLASTOUT2ABUNDANCE_CPP) 

        def compile_blastout2abundance():
            call(cmd)

        self.execute(compile_blastout2abundance, [], 'Compiling blastout2abundance')

        # build extract_mapped_reads
        cmd = [os.environ['CXX'],'-std=c++11','-o',os.path.join(EXTRACTMAPREAD_PATH,'extract_mapped_reads')]
        cmd.extend(EXTRACTMAPREAD_CPP)

        def compile_extract_mapped_reads():
            call(cmd)

        self.execute(compile_extract_mapped_reads, [], 'Compiling extract_mapped_reads')

        # copy resulting tool to library build folder
        self.mkpath(self.build_lib)

        if not self.dry_run:
            self.copy_file(os.path.join(STRAINCALL_PATH,'StrainCall'),self.build_lib)
            self.copy_file(os.path.join(BLASTOUT2ABUNDANCE_PATH,'blastout2abundance'),self.build_lib)
            self.copy_file(os.path.join(EXTRACTMAPREAD_PATH,'extract_mapped_reads'),self.build_lib)


class ExInstall(install):
    def initialize_options(self):
        install.initialize_options(self)
        self.build_scripts = None

    def finalize_options(self):
        install.finalize_options(self)
        self.set_undefined_options('build', ('build_scripts', 'build_scripts'))

    def run(self):
        # run original install code
        install.run(self)

        # install StrainCall executables
        self.copy_tree(self.build_lib, os.path.join(self.exec_prefix,'bin'))
            
        # install copy number file
        self.copy_file(os.path.join(BASEPATH,'scripts/rrnDB_RDP.tsv'), os.path.join(self.exec_prefix,'bin/rrnDB_RDP.tsv'))

install_requires = ['scipy >= 0.15.1','numpy >= 1.9.2','ete2 >= 2.3.1','mpi4py >= 1.3.1']
scripts = glob.glob('scripts/*.py')

setup(
    name = 'rbra',
    description = 'a pipeline to assemble strain-level full-length 16S rRNA genes',
    author = 'Feng Zeng',
    author_email = 'zengfeng@xmu.edu.cn',
    url = '',
    scripts = scripts,
    install_requires = install_requires,   
    cmdclass={
        'build': ExBuild,
        'install': ExInstall,
    }
)

