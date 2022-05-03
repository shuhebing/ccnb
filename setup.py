from setuptools import setup, find_packages


def readme():
    with open('README.en.md', encoding='utf-8') as f:
        content = f.read()
    return content


install_requires = [
    "ase==3.15.0",
    "pymatgen==2022.0.10",
    "cavd==0.5.4",
]

setup(name="ccnb",
      version="0.2.2",
      author='He Bing',
      author_email='bhe@t.shu.edu.cn',
      long_description=readme(),
      long_description_content_type='text/markdown',
      url='https://gitee.com/shuhebing/CCNB',
      packages=find_packages(),
      include_package_data=True,
      install_requires=install_requires,
      python_requires='>=3.7,<3.10',
      entry_points={
          'console_scripts': [
              'ccnb = ccnb:main.main',
          ],
      })
