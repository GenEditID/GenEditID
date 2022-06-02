from setuptools import setup, find_packages

version = '0.3'

print("""
-----------------------------------
 Installing GeneditID version {}
-----------------------------------
""".format(version))

install_requires = [
    # configuration
    'pyyaml==5.4',
    # database
    # 'psycopg2',  # for production only
    'sqlalchemy==1.3.20',
    # test
    'pytest==6.2.0',
    # data
    'xlrd==1.2.0',
    'pandas==1.1.5',
    # plots
    'plotly==4.6.0',
    'colorlover==0.3.0',
    # geneditidapp
    'pyramid==1.10.5',
    'pyramid_chameleon==0.3',
    'pyramid_jinja2==2.8',
    'pyramid_debugtoolbar==4.9',
    'pyramid_tm==2.4',
    'deform==2.0.15',
    'transaction==3.0.1',
    'zope.sqlalchemy==1.3',
    'waitress==2.1.2',
    # amplicount and variantid analysis
    'pyfaidx==0.5.8',
    'biopython==1.78',
    'varcode==1.0.3'
]

setup(
    name='geneditid',
    version=version,
    author='GeneditID team',
    description='GeneditID database, WebApp and tools to identify variants and create plate-base visualisation to identify clones with desirable gene edits',
    url='https://geneditid.github.io/',
    license='MIT License',
    packages=find_packages(),
    python_requires='>=3.0, <4',
    include_package_data=True,
    install_requires=install_requires,
    entry_points={
        'paste.app_factory': [
            'main = geneditidapp:main',
        ],
        'console_scripts': [
            'geneditid_run_ampli_find = geneditidtools.run_ampli_find:main',
            'geneditid_run_amplicount = geneditidtools.run_ampli_count:main',
            'geneditid_run_variantid = geneditidtools.run_ampli_plots:main',
        ],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 3',
        'Framework :: Pyramid',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],

)

print("""
--------------------------------
GeneditID installation complete!
--------------------------------
For help in running GeneditID, please see the documentation available
at https://geneditid.github.io/
""")
