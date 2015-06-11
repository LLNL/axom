#!/bin/env python
from __future__ import print_function


import collections
import string
import json

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


class Argument(object):
    """ Collect fields for an argument.
    This used to be a dict but a class has better access semantics: i.attr vs d['attr']
    It also initializes default values to avoid  d.get('attr', default)
    """
    # valid fields
    fields = ['name', 'type', 'other']
    defaults = dict(
        type='dtype',
        other='dother'
        )

    def __init__(self, name, **kw):
        self.name = name
#        for key, defvalue in self.defaults.items():
#            setattr(self, key, defvalue)
        self.__dict__.update(self.defaults) # set all default values
        for key in kw:
            if key in self.defaults:
                setattr(self, key, kw[key])
            else:
                raise RuntimeError("Unknown key for Argument %s", key)

    def _to_dict(self):
        """Convert instance to a dictionary for json.
        """
        a = {}
        for key in self.fields:
            a[key] = getattr(self, key)
        return a


class ExpandedEncoder(json.JSONEncoder):
    def default(self, obj):
        if hasattr(obj, '_to_dict'):
            return obj._to_dict()
         # Let the base class default method raise the TypeError
        return json.JSONEncoder.default(self, obj)


class Option(object):
    """
    If attribute is not found, look in parent's.
    A replacement for a dictionary to allow obj.name syntax.
    It will automatically look in __parent to attribute if not found to allow
    A nesting of options.
    Use __attr to avoid xporting to json
    """
    def __init__(self, parent, **kw):
        self.__parent = parent
        self.__hidden = 43
        for key, value in kw.items():
            setattr(self, key, value)

    def __getattr__(self, name):
        # we get here if the attribute does not exist in current instance
        if self.__parent:
            return getattr(self.__parent, name)
        else:
            raise AttributeError("%r object has no attribute %r" %
                                 (self.__class__.__name__, name))

    def _to_dict(self):
        d = {}
        skip = '_' + self.__class__.__name__ + '__'
        for key, value in self.__dict__.items():
            if not key.startswith(skip):
                d[key] = value
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


    # Argument
    print("Test Argument")
    a = Argument('abc', type='def') #, bird='abcd')
    print(json.dumps(a, cls=ExpandedEncoder, sort_keys=True))

    # Option
    print("Test Option")
    lev0 = Option(None, a=1, b=2, c=3)
    lev1 = Option(lev0, x=100, y=01, z=102)
    lev0.c2 = 32
    lev1.z2 = 103

    print(lev0.a)
    print(lev0.c2)
    #print(lev0.z)
    print(lev1.a)
    print(lev1.z)
    print(lev1.c2)
    print(lev1.z2)
    try:
        print(lev1.nosuch)
    except AttributeError:
        print("Passed nosuch attribute")

    print(json.dumps(lev0, cls=ExpandedEncoder, sort_keys=True))
    print(json.dumps(lev1, cls=ExpandedEncoder, sort_keys=True))
