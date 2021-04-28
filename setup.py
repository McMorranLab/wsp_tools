from setuptools import setup, find_packages

setup(
	name = 'wsp_tools-test-0',
	packages = find_packages(),
	version = '1.2.02',
	author = 'William S. Parker',
	author_email = 'wparker4@uoregon.edu',
	description = 'Utilities for TEM data analysis and simulation.',
	url = 'https://github.com/McMorranLab/wsp_tools',
	project_urls={
        "Bug Tracker": "https://github.com/pypa/sampleproject/issues",
    },
	classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
	long_description = open('README.md').read(),
	long_description_content_type = "text/markdown",
	python_requires='>=3.6',
	install_requires=['numpy','matplotlib','scipy','pdoc3','deprecated']
)
