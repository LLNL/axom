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
