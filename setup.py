# The first version of setuptools capable of dealing with a pyproject.toml
# See https://pip.pypa.io/en/stable/reference/pip/#pep-517-and-518-support

# We *could* in addition, specify the first version of pip capable of recognizing `python_requires` (9.0.1)
# See https://packaging.python.org/guides/distributing-packages-using-setuptools/#python-requires
# But that fails under a conda environment (as does any approach that does `import pip`)
# so this is best we can do for now.

import pkg_resources
pkg_resources.require(['setuptools >= 40.8.0'])

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
    version='0.3.3',
    packages=['hatchet', 'hatchet.utils', 'hatchet.bin'],
    package_dir={'': 'src'},
    package_data={'hatchet': ['hatchet.ini']},
    ext_modules=[CMakeExtension('hatchet.solve')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,

    python_requires='>=3.7',

    entry_points={
        'console_scripts': [
            'hatchet = hatchet.__main__:main',
        ],
    },

    install_requires=[
        'biopython',
        'matplotlib',
        'pandas',
        'psutil',
        'pysam',
        'requests',
        'seaborn',
        'scikit-learn',
        'scipy'
    ],

    extras_require={
        'dev': ['pytest', 'mock', 'numpydoc', 'sphinx', 'sphinxcontrib-bibtex<2.0.0', 'sphinx-rtd-theme', 'recommonmark', 'sphinx-markdown-tables']
    }

)
