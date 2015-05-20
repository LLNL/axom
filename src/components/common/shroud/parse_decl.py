"""
Parse a declaration and return a named-tuple tree.

decl = ( type, name, attrs, [ args ] )

arg = ( type, name, attrs )

attrs = { key : value }


This module just parses syntax.  Semantics, like valid type names
and duplicate argument names are check later.
"""


import parsley

def add_to_dict(d, key, value):
    d[key] = value
    return d

x = parsley.makeGrammar("""
name = < (letter | '_') (letter | digit | '_')* >

#type = name:t ?( t in types ) ^(C-type) -> t
type = name:t

string = (('"' | "'"):q <(~exactly(q) anything)*>:xs exactly(q))
                     -> xs

integer = <digit*>:i -> int(i)

value = name | string | integer

attr = '+' name:n ( '=' value | -> True ):v
        -> (n,v)

qualifier = 'const' -> [('const', True)]
            |       -> []

pointer = '*' -> [('ptr', True)]
        | '&' -> [('reference', True)]
        |     -> []

declarator = qualifier:qu ws type:t ws pointer:pp ws name:n  ( ws attr )*:at
        -> dict(type=t, name=n, attrs=dict(qu+pp+at))

parameter_list = declarator:first ( ws ',' ws declarator)*:rest -> [first] + rest
                 | -> []

argument_list = ( '(' ws parameter_list:l ws ')' ) -> l
                | -> []

decl = declarator:dd ws argument_list:args
        -> dict( result=dd, args=args)
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
        ]:
        r = x(test).parameter_list()
        print('parameter_list: "{0}"'.format(test))
        json.dump(r, sys.stdout, sort_keys=True, indent=4)
        print('\n')

    for test in [
        "()",
        "(int arg1)",
        "(int arg1, double arg2)",
        ] :
        r = x(test).argument_list()
        print('argument_list: "{0}"'.format(test))
        json.dump(r, sys.stdout, sort_keys=True, indent=4)
        print('\n')

    for test in [
        "void foo",
        "void foo +alias=junk",
        "void foo()",
        "void foo(int arg1)",
        "void foo(int arg1, double arg2)",
        "const void foo(int arg1+in, double arg2+out)",
        ]:
        r = x(test).decl()
        print('decl: "{0}"'.format(test))
        json.dump(r, sys.stdout, sort_keys=True, indent=4)
        print('\n')


