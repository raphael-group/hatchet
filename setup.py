import os
import re
import sys
import platform
import subprocess
import shutil

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_INSTALL_PREFIX=' + extdir,
                      '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DGUROBI_HOME=' + os.environ.get('GUROBI_HOME', ''),
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        # TODO: "Release" generates a solve binary that coredumps on Linux - investigate.
        cfg = 'Debug' if self.debug else ''  # else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if os.path.exists(self.build_temp):
            shutil.rmtree(self.build_temp)
        os.makedirs(self.build_temp)

        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.', '--target', 'install'] + build_args, cwd=self.build_temp)


setup(
    name='hatchet',
    version='0.0.1',
    packages=['hatchet', 'hatchet.utils', 'hatchet.bin'],
    package_dir={'': 'src'},
    package_data={'hatchet': ['hatchet.ini']},
    ext_modules=[CMakeExtension('hatchet.solve')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,

    python_requires='~=2.7',

    install_requires=[
        'cython',
        'biopython==1.76',
        'bnpy',
        'futures',
        'importlib_resources>=1.0.2',
        'matplotlib',
        'matplotlib-venn',
        'munkres<=1.0.12',
        'opencv-python<=4.3.0.36',
        'pandas',
        'psutil',
        'pysam',
        'seaborn',
        'scikit-learn',
        'scipy'
    ],

    extras_require = {
        'dev': ['pytest', 'mock']
    }

)
