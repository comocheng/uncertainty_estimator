# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

import re

my_comment = 'Gas phase thermo for [CH]CC from Thermo group additivity estimation: group(Cs-CsCsHH)'
print(my_comment)

new_start = my_comment.find('from ', len('Gas phase thermo for ')) + len('from ')

my_comment[new_start:]

my_comment.re

regex = r'Gase phase thermo for . from '

re.split(regex, my_comment)


