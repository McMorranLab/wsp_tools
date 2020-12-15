from setuptools import setup, find_packages

setup(
    name = 'wsp_tools',
    packages = find_packages(),
	version = '1.0.4',
    author = 'William Parker',
    author_email = 'wparker4@uoregon.edu',
    description = 'McMorran Lab tools developed by WSP. ',
    url = 'https://github.com/McMorranLab/wsp_tools',
    long_description = open('README.md').read(),
    long_description_content_type = "text/markdown",
    python_requires='>=3.6',
    install_requires=['numpy','matplotlib','scipy']
)
