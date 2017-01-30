"""
Read a file and extract the splicer blocks.
"""
from __future__ import print_function
from __future__ import absolute_import

import os


def get_splicers(fname, out):
    """
    fname - input file name
    out - dictionary to update

    tags of the form
    begin  value ...
    """
    str_begin = 'splicer begin'
    str_end = 'splicer end'
    str_push = 'splicer push'
    str_pop = 'splicer pop'

    state_look = 1
    state_collect = 2
    state = state_look

    top = out
    stack = [top]

    begin_tag = ''

    with open(fname, 'r') as fp:
        for line in fp.readlines():
            if state == state_look:
                i = line.find(str_begin)
                if i > 0:
                    fields = line[i+len(str_begin):].split()
                    tag = fields[0]
                    begin_tag = tag
                    subtags = tag.split('.')
                    begin_subtag = subtags[-1]
                    for subtag in subtags[:-1]:
                        top = top.setdefault(subtag, {})
                        stack.append(top)
#                    print("BEGIN", begin_tag)
                    save = []
                    state = state_collect
                    continue

            elif state == state_collect:
                i = line.find(str_end)
                if i > 0:
                    fields = line[i+len(str_end):].split()
                    end_tag = fields[0]
#                    print("END", end_tag)
                    if begin_tag != end_tag:
                        raise RuntimeError(
                            "Mismatched tags  '%s' '%s'", (begin_tag, end_tag))
                    if end_tag in top:
                        raise RuntimeError(
                            "Tag already exists - '%s'" % begin_tag)
                    top[begin_subtag] = save
                    top = out
                    stack = [top]
                    state = state_look
                else:
                    save.append(line.rstrip())


def get_splicer_based_on_suffix(name, out):
        fileName, fileExtension = os.path.splitext(name)
        if fileExtension in ['.f', '.f90']:
            d = out.setdefault('f', {})
            get_splicers(name, d)
        elif fileExtension in [
                '.c', '.cpp', '.hpp', '.cxx', '.hxx', '.cc', '.C']:
            d = out.setdefault('c', {})
            get_splicers(name, d)
        elif fileExtension in ['.py']:
            d = out.setdefault('py', {})
            get_splicers(config, name, d)
        elif fileExtension in ['.lua']:
            d = out.setdefault('lua', {})
            get_splicers(config, name, d)


# def print_tree(out):
#     import json
#     print(json.dumps(out, indent=4, sort_keys=True))

if __name__ == '__main__':
    import glob
    import json

    out = {}
    for name in glob.glob('../tests/example/*'):
        get_splicer_based_on_suffix(name, out)
    print(json.dumps(out, indent=4, sort_keys=True))

    for name in ['fsplicer.f', 'csplicer.c', 'pysplicer.c']:
        out = {}
        get_splicers(os.path.join('..', 'tests', name), out)
        print(json.dumps(out, indent=4, sort_keys=True))
