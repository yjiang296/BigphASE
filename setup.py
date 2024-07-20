from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install
import sys
import os
import subprocess
import shutil
class BuildExt(build_ext):
    def run(self):
        for ext in self.extensions:
            self.build_extension(ext)
        self.copy_executable()

    def build_extension(self, ext):
        sources = ext.sources or []
        for source in sources:
            if source.endswith('.cpp'):
                self.compile_cpp(source)

    def compile_cpp(self, source):
        # Ensure g++ version is >= 4.8 (which supports C++11)
        min_gcc_version = (4, 8)

        # Check g++ version
        try:
            output = subprocess.check_output(['g++', '--version']).decode()
            version_line = output.splitlines()[0]
            version_str = version_line.split()[-1]
            version_parts = version_str.split('.')
            gcc_version = tuple(map(int, version_parts[:2]))
        except (subprocess.CalledProcessError, IndexError, ValueError):
            raise RuntimeError('Failed to determine g++ version')

        if gcc_version < min_gcc_version:
            raise RuntimeError(f'g++ version {gcc_version[0]}.{gcc_version[1]} is too old. '
                               f'Please update to g++ >= {min_gcc_version[0]}.{min_gcc_version[1]}.')

        # Compile using g++ with C++11 standard
        object_file = source.replace('.cpp', '.o')
        include_dirs = ['src/BigphASE/delta2hap/include']  
        compile_cmd = ['g++', '-std=c++11', '-c', source, '-o', object_file] + \
                      [f'-I{dir}' for dir in include_dirs]
        print(f'Compiling {source}...')
        sys.stdout.flush()
        if os.system(' '.join(compile_cmd)) != 0:
            raise RuntimeError('Compilation failed')

    def copy_executable(self):
        # Link the object files and create the executable
        object_files = [
            'src/BigphASE/delta2hap/src/delta2hap.o',
            'src/BigphASE/delta2hap/src/delta.o',
            'src/BigphASE/delta2hap/src/tigrinc.o',
            'src/BigphASE/delta2hap/src/translate.o'
        ]
        
        output_executable = 'src/BigphASE/delta2hap/bin/delta2hap'
        # Ensure the bin directory exists
        os.makedirs(os.path.dirname(output_executable), exist_ok=True)
        
        link_cmd = ['g++', '-o', output_executable] + object_files
        print(f'Linking {output_executable}...')
        sys.stdout.flush()
        if os.system(' '.join(link_cmd)) != 0:
            raise RuntimeError('Linking failed')
        
        # Move the executable to the bin directory
        if not os.path.exists(output_executable):
            raise RuntimeError('Executable not found after compilation')

class InstallCommand(install):
    def run(self):
        self.run_command('build_ext')
        super().run()
        self.cleanup()
        
    def cleanup(self):
        # Define the directories and files to be removed
        cleanup_paths = [
            'src/BigphASE/delta2hap/include',
            'src/BigphASE/delta2hap/src',
            'src/BigphASE/delta2hap/CMakeLists.txt'
        ]

        for path in cleanup_paths:
            full_path = os.path.abspath(path)
            if os.path.isdir(full_path):
                shutil.rmtree(full_path)
                print(f'Removed directory {full_path}')
            elif os.path.isfile(full_path):
                os.remove(full_path)
                print(f'Removed file {full_path}')



cpp_extension = Extension('delta2hap',
                          sources=['src/BigphASE/delta2hap/src/delta2hap.cpp',
                                   'src/BigphASE/delta2hap/src/delta.cpp',
                                   'src/BigphASE/delta2hap/src/tigrinc.cpp',
                                   'src/BigphASE/delta2hap/src/translate.cpp'])

setup(
    name='BigphASE',
    version='0.1.1',
    description='This package was devised to analyze hybrid RNA-seq using bi-parental graph strategy.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Yi Jiang, Yi Mu Ran , Zhule Liu',
    author_email='jiangy296@mail2.sysu.edu.cn',
    url='https://github.com/yjiang296/BigphASE',
    license='GNU GPL v3.0 License',
    platforms='any',
    ext_modules=[cpp_extension],
    cmdclass={'build_ext': BuildExt, 'install': InstallCommand},
    packages=['BigphASE'],
    include_package_data=True,
    package_dir={'': 'src'},
    install_requires=[],
    classifiers=[
        'Natural Language :: English',
        'Programming Language :: Python :: 3.12',
    ],
)
