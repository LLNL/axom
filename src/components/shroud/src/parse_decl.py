"""
Parse a declaration and return a named-tuple tree.

decl = ( type, name, attrs, [ args ] )

arg = ( type, name, attrs )

attrs = { key : value }


This module just parses syntax.  Semantics, like valid type names
and duplicate argument names are check later.
"""
from __future__ import print_function
from __future__ import absolute_import

import parsley


def add_to_dict(d, key, value):
    d[key] = value
    return d


x = parsley.makeGrammar("""
name = < (letter | '_') (letter | digit | '_' | ':')* >

#type = name:t ?( t in types ) ^(C-type) -> t
type = name:t

digits = <digit*>
floatPart :sign :ds = <('.' digits exponent?) | exponent>:tail
                     -> float(sign + ds + tail)
exponent = ('e' | 'E') ('+' | '-')? digits

number = spaces ('-' | -> ''):sign (digits:ds (floatPart(sign ds)
                                               | -> int(sign + ds)))

string = (('"' | "'"):q <(~exactly(q) anything)*>:xs exactly(q))
                     -> xs

parens = <'('  (~')' anything)* ')'>

value = name | string | number

attr = '+' ws name:n ( ('=' value) | parens | -> True ):v
        -> (n,v)

qualifier = 'const' -> [('const', True)]
            |       -> []

pointer = '*' -> [('ptr', True)]
        | '&' -> [('reference', True)]
        |     -> []

default = ws '=' ws value:default -> [('default', default)]
                 |                -> []

declarator = qualifier:qu ws type:t ws pointer:pp ws name:n  ( ws attr )*:at default:df
        -> dict(type=t, name=n, attrs=dict(qu+pp+at+df))

parameter_list = declarator:first ( ws ',' ws declarator)*:rest -> [first] + rest
                 | -> []

argument_list = ( '(' ws parameter_list:l ws ')' ) -> l
                | -> []

decl = declarator:dd ws argument_list:args ws qualifier:qual (ws attr)*:at
        -> dict( result=dd, args=args, attrs=dict(qual + at))
""", {})


def check_decl(expr, parser=x):
    """ parse expr as a declaration, return list/dict result.
    """
    return parser(expr).decl()


if __name__ == '__main__':
    import sys
    import json
    for test in [
        "const",
        "",
    ]:
        r = x(test).qualifier()
        print('qualifier: "{0}"'.format(test))
        json.dump(r, sys.stdout, sort_keys=True, indent=4)
        print('\n')

    for test in [
        "*",
        "&",
        "",
    ]:
        r = x(test).pointer()
        print('pointer: "{0}"'.format(test))
        print(r)
        json.dump(r, sys.stdout, sort_keys=True, indent=4)
        print('\n')

    for test in [
        "+intent",
        "+intent=in",
        '+name="abcd"',
        "+name='def'",
        "+ii=12",
        "+d1=-12.0",
        "+d2=11.3e-10",
        "+d3=11e10",
        "+intent()",
        "+intent(in)",
        "+dimension",
        "+dimension(*)",
        "+dimension(len)",
    ]:
        r = x(test).attr()
        print('attr: "{0}"'.format(test))
        json.dump(r, sys.stdout, sort_keys=True, indent=4)
        print('\n')

    for test in [
        "int arg",
        "const int arg",
        "badtype arg",
    ]:
        r = x(test).declarator()
        print('declarator: "{0}"'.format(test))
        json.dump(r, sys.stdout, sort_keys=True, indent=4)
        print('\n')

    for test in [
        'int arg',
        'int *arg',
        'int arg1, double arg2',
        'int arg +in',
        'int arg +in +value',
        'const string& getName',
        'std::string getName',
    ]:
        r = x(test).parameter_list()
        print('parameter_list: "{0}"'.format(test))
        json.dump(r, sys.stdout, sort_keys=True, indent=4)
        print('\n')

    for test in [
        "()",
        "(int arg1)",
        "(int arg1, double arg2)",
        "(int arg1, double arg2 = 0.0)",
    ]:
        r = x(test).argument_list()
        print('argument_list: "{0}"'.format(test))
        json.dump(r, sys.stdout, sort_keys=True, indent=4)
        print('\n')

    for test in [
        "void foo",
        "void foo +alias=junk",
        "void foo()",
        "void foo() const",
        "void foo(int arg1)",
        "void foo(int arg1, double arg2)",
        "const std::string& getName() const",
        "const void foo(int arg1+in, double arg2+out)",
        "void new() + constructor",
    ]:
        r = x(test).decl()
        print('decl: "{0}"'.format(test))
        json.dump(r, sys.stdout, sort_keys=True, indent=4)
        print('\n')
