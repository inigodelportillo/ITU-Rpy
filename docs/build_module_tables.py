# -*- coding: utf-8 -*-
import re
import requests
import importlib
from functools import reduce


def make_rst_table(grid):
    cell_width = 2 + max(reduce(lambda x,y: x+y, [[len(item) for item in row] for row in grid], []))
    num_cols = len(grid[0])
    rst = table_div(num_cols, cell_width, 0)
    header_flag = 1
    for row in grid:
        new_row = '| ' + '| '.join([normalize_cell(x, cell_width-1) for x in row])
        rst = rst + new_row
        rst = rst + ' ' * ((cell_width + 1) * num_cols - len(new_row)) + '|\n'
        rst = rst + table_div(num_cols, cell_width, header_flag)
        header_flag = 0
    return rst


def table_div(num_cols, col_width, header_flag):
    if header_flag == 1:
        return num_cols*('+' + (col_width)*'=') + '+\n'
    else:
        return num_cols*('+' + (col_width)*'-') + '+\n'


def normalize_cell(string, length):
    return string + ((length - len(string)) * ' ')


URL = 'https://www.itu.int/rec/R-REC-P.{0}/en'
URL_PDF = 'https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.{0}!!PDF-E.pdf'
RECs = [618, 676, 453, 530, 835, 836, 837, 838, 839, 840, 1510, 1511]
# RECs = [840]

for rec in RECs:
    rec = str(rec)
    url = URL.format(rec)
    resp = requests.get(url)
    txt = resp.text.replace(u'\xa0', u' ')
    header = re.findall('<strong>(.+)</strong>', txt)[0]
    title = header.split(':')[0].strip().capitalize()
    desc = header.split(':')[-1].strip()
    approved_in = re.findall('Approved in (.+)</p>', txt)[0]

    recs_and_links = re.findall('<a href="(.+)"><strong>(.+)</strong></a>', txt)

    # Latest version
    link_latest, name_latest = recs_and_links[0]
    name_latest, date_latest = name_latest.strip().split(' ')
    date_latest = date_latest.replace('(', '').replace(')', '')
    code_latest = rec + link_latest.split(rec)[-1]
    link_latest = URL.format(code_latest)
    pdf_latest = URL_PDF.format(code_latest)

    table = [['Title', 'PDF', 'Latest approved in'],
             ['`Recommendation {0} <{1}>`_'.format('ITU-R ' + title, url),
              '`[PDF] <{0}>`_'.format(pdf_latest), approved_in],
             ['{0}'.format(desc)],
             ['**Current recommendation version (In force)**', '', '**Date**'],
             ['`Recommendation {0} <{1}>`_'.format('ITU-R ' + name_latest, link_latest),
              '`[PDF] <{0}>`_'.format(pdf_latest), date_latest]]

    rec_implemented = [
             ['**Recommendations implemented in ITU-Rpy**', '', '**Date**']]

    rec_not_implemented = [
             ['**Recommendations not implemented in ITU-Rpy**', '', '**Date**']]

    module = importlib.import_module('itur.models.itu{0}'.format(rec))

    recs_and_links = sorted(recs_and_links, key=lambda x: (len(x[1]), x))[::-1]
    for rec_link, rec_name in recs_and_links:
        name, date = rec_name.split(' ')
        date = date.replace('.', '/').replace('(', '').replace(')', '')
        code_rec = rec + rec_link.split(rec)[-1]

        link = URL.format(code_rec)
        pdf_link = URL_PDF.format(code_rec)
        version = int(name.split('-')[-1].strip())

        if len(date) == 5 and date[3] == '9':
            date = date[:3] + '19' + date[3:]
        elif len(date) == 5 and date[3] in ['0', '1', '2']:
            date = date[:3] + '20' + date[3:]

        # Check if version is available in ITU-Rpy
        try:
            module.change_version(version)
            rec_implemented.append([
              '`Recommendation {0} <{1}>`_'.format('ITU-R ' + name, link),
              '`[PDF] <{0}>`_'.format(pdf_link),
              date])
        except Exception:
            rec_not_implemented.append([
              '`Recommendation {0} <{1}>`_'.format('ITU-R ' + name, link),
              '`[PDF] <{0}>`_'.format(pdf_link),
              date])

    # Join tables that have rows
    if len(rec_implemented) > 1:
        table += rec_implemented

    if len(rec_not_implemented) > 1:
        table += rec_not_implemented

    # Write table to a file
    table_str = make_rst_table(table)
    with open('./apidoc/itu{0}_table.rst'.format(rec), 'w') as fd:
        fd.write(table_str)
