from platform import python_version
from setuptools import setup, find_packages


def readme():
    with open('README.md', encoding='utf-8') as f:
        content = f.read()
    return content


setup(name="CCNB",
      version="0.1",
      author='Bing He',
      author_email='bhe@t.shu.edu.cn',
      long_description=readme(),
      long_description_content_type='text/markdown',
      url='https://gitee.com/shuhebing/CCNB',
      packages=find_packages(),
      include_package_data=True,
      entry_points={
          'console_scripts': [
              'ccnb = ccnb:main.main',
          ],
      })
