#!/bin/env python3

import collections
import string

fmt = string.Formatter()

def wformat(template, d):
    # shorthand, wrap format
    return fmt.vformat(template, None, d)

def eval_templates(templates, node, fmt_dict):
    """ If variable is not already set in node, then
    set based on template.
    Update value in fmt_dict.
    """
    for tname in ( templates ) :
        if tname not in node:
            node[tname] = fmt.vformat(
                node['options'][tname + '_template'],
                None, fmt_dict)
        fmt_dict[tname] = node[tname]


# http://stackoverflow.com/questions/1175208/elegant-python-function-to-convert-camelcase-to-camel-case
def un_camel(text):
    """ Converts a CamelCase name into an under_score name. 

        >>> un_camel('CamelCase')
        'camel_case'
        >>> un_camel('getHTTPResponseCode')
        'get_http_response_code'
    """
    result = []
    pos = 0
    while pos < len(text):
        if text[pos].isupper():
            if pos-1 > 0 and text[pos-1].islower() or pos-1 > 0 and \
            pos+1 < len(text) and text[pos+1].islower():
                result.append("_%s" % text[pos].lower())
            else:
                result.append(text[pos].lower())
        else:
            result.append(text[pos])
        pos += 1
    return "".join(result)

# http://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth
def update(d, u):
    for k, v in u.items():
        if isinstance(v, collections.Mapping):
            r = update(d.get(k, {}), v)
            d[k] = r
        else:
            d[k] = u[k]
    return d



if __name__ == '__main__':
    print(un_camel('incrementCount'))
    print(un_camel('local_function1'))
    print(un_camel('getHTTPResponseCode'))


    a = dict(name='aa', type='int', attrs=dict(ptr=True))
    b = dict(attrs=dict(len='somefunction'))

    print("AAAA", a)
    print("BBBB", b)
#    b.update(a)
    update(b, a)
    print("CCCC", b)
