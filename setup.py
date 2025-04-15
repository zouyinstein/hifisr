from setuptools import setup, find_packages

setup(
    name='hifisr',
    version='0.5.0',
    author="Yi Zou",
    author_email="zouyi.nju@gmail.com",
    url='https://github.com/zouyinstein/hifisr',
    packages=find_packages(),
    py_modules=['hifisr_functions'],
    install_requires=[
        'biopython',
        'pysam',
        'pandas',
        'numpy',
        'openpyxl',
        'xlsxwriter',
        'matplotlib',
        'polars',
        'fastexcel',
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)
