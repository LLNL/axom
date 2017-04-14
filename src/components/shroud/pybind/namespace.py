mod = Module('MyModule')
outer = mod.add_cpp_namespace('Outer')
outer.add_class('MyClass')
outer.add_function('Do', None, [])
inner = outer.add_cpp_namespace('Inner')
inner.add_class('MyClass')
inner.add_function('Do', None, [])

