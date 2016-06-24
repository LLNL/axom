#
# Routines to generate splicers for wrappers.
# Used to generate several variations of a routine for Fortran.
# Similar to templates in C++.
#
from __future__ import print_function
import sys

# types to use for generic routines
types = (
    ( 'int',    'integer(C_INT)',  'SIDRE_INT_ID'),
    ( 'long',   'integer(C_LONG)', 'SIDRE_LONG_ID'),
    ( 'float',  'real(C_FLOAT)',   'SIDRE_FLOAT_ID'),
    ( 'double', 'real(C_DOUBLE)',  'SIDRE_DOUBLE_ID'),
)

# maximum number of dimensions of generic routines
maxdims = 4

def XXnum_metabuffers():
    return len(types) * (maxdims + 1) # include scalars
######################################################################

def group_get_scalar(d):
    """Create methods on DataGroup to get a scalar.
    """
    return """
! Generated by genfsidresplicer.py
subroutine datagroup_get_scalar_{typename}(group, name, value)
    use iso_c_binding
    class(datagroup), intent(IN) :: group
    character(*), intent(IN) :: name
    {f_type}, intent(OUT) :: value
    integer(C_INT) :: lname
    type(C_PTR) view

    lname = len_trim(name)
    view = c_datagroup_get_view_from_name_bufferify(group%voidptr, name, lname)
    value = c_dataview_get_data_{typename}(view)
end subroutine datagroup_get_scalar_{typename}""".format(**d)

def group_set_scalar(d):
    """Create methods on DataGroup to set a scalar.
    """
    return """
! Generated by genfsidresplicer.py
subroutine datagroup_set_scalar_{typename}(group, name, value)
    use iso_c_binding
    class(datagroup), intent(IN) :: group
    character(*), intent(IN) :: name
    {f_type}, intent(OUT) :: value
    integer(C_INT) :: lname
    type(C_PTR) view

    lname = len_trim(name)
    view = c_datagroup_get_view_from_name_bufferify(group%voidptr, name, lname)
    call c_dataview_set_scalar_{typename}(view, value)
end subroutine datagroup_set_scalar_{typename}""".format(**d)

def group_create_array_view(d):
    # typename - part of function name
    # nd       - number of dimensions
    # f_type   - fortran type
    # shape     - :,:, to match nd
    if d['rank'] == 0:
        extents_decl = 'extents(1)'
        extents_asgn = 'extents(1) = 1_SIDRE_LENGTH'
    else:
        extents_decl = 'extents(%d)' % d['rank']
        extents_asgn = 'extents = shape(value, kind=SIDRE_LENGTH)'

    return """
! Generated by genfsidresplicer.py
function datagroup_create_array_view_{typename}{nd}(group, name, value) result(rv)
    use iso_c_binding
    implicit none

    class(datagroup), intent(IN) :: group
    character(*), intent(IN) :: name
    {f_type}, target, intent(IN) :: value{shape}
    integer(C_INT) :: lname
    type(dataview) :: rv
    integer(SIDRE_LENGTH) :: {extents_decl}
    integer(C_INT), parameter :: type = {sidre_type}
    type(C_PTR) addr

    lname = len_trim(name)
    call SHROUD_C_LOC(value, addr)
    if (c_associated(addr)) then
      {extents_asgn}
      rv%voidptr = c_datagroup_create_view_external_bufferify( &
          group%voidptr, name, lname, addr)
      call c_dataview_apply_type_shape(rv%voidptr, type, {rank}, extents)
    else
      rv%voidptr = c_datagroup_create_view_from_type_bufferify( &
          group%voidptr, name, lname, type, 0_C_LONG)
    endif
end function datagroup_create_array_view_{typename}{nd}""".format(
        extents_decl=extents_decl,
        extents_asgn=extents_asgn, **d)

def group_set_array_data_ptr(d):
    """
    call view%set_external_data_ptr
    hide c_loc call and add target attribute
    """
    # XXX - should this check the type/shape of value against the view?
    # typename - part of function name
    # nd       - number of dimensions
    # f_type   - fortran type
    # shape     - :,:, to match nd
    if d['rank'] == 0:
        extents_decl = 'extents(1)'
        extents_asgn = 'extents(1) = 1_SIDRE_LENGTH'
    else:
        extents_decl = 'extents(%d)' % d['rank']
        extents_asgn = 'extents = shape(value, kind=SIDRE_LENGTH)'

    return """
! Generated by genfsidresplicer.py
! This function does nothing if view name does not exist in group.
subroutine datagroup_set_array_data_ptr_{typename}{nd}(group, name, value)
    use iso_c_binding
    implicit none

    class(datagroup), intent(IN) :: group
    character(len=*), intent(IN) :: name
    {f_type}, target, intent(IN) :: value{shape}
    integer(C_INT) :: lname
    type(C_ptr) view
!    integer(SIDRE_LENGTH) :: {extents_decl}
!    integer(C_INT), parameter :: type = {sidre_type}
    type(C_PTR) addr

    lname = len_trim(name)
!    {extents_asgn}
    view = c_datagroup_get_view_from_name_bufferify(group%voidptr, name, lname)
    if (c_associated(view)) then
        call SHROUD_C_LOC(value, addr)
        call c_dataview_set_external_data_ptr_only(view, addr)
!        call c_dataview_apply_type_shape(rv%voidptr, type, {rank}, extents)
    endif
end subroutine datagroup_set_array_data_ptr_{typename}{nd}""".format(
        extents_decl=extents_decl,
        extents_asgn=extents_asgn, **d)

def view_set_array_data_ptr(d):
    """
    call view%set_external_data_ptr
    hide c_loc call and add target attribute
    """
    # XXX - should this check the type/shape of value against the view?
    # typename - part of function name
    # nd       - number of dimensions
    # f_type   - fortran type
    # shape     - :,:, to match nd
    if d['rank'] == 0:
        extents_decl = 'extents(1)'
        extents_asgn = 'extents(1) = 1_SIDRE_LENGTH'
    else:
        extents_decl = 'extents(%d)' % d['rank']
        extents_asgn = 'extents = shape(value, kind=SIDRE_LENGTH)'

    return """
! Generated by genfsidresplicer.py
subroutine dataview_set_array_data_ptr_{typename}{nd}(view, value)
    use iso_c_binding
    implicit none

    class(dataview), intent(IN) :: view
    {f_type}, target, intent(IN) :: value{shape}
!    integer(SIDRE_LENGTH) :: {extents_decl}
!    integer(C_INT), parameter :: type = {sidre_type}
    type(C_PTR) addr

!    lname = len_trim(name)
!    {extents_asgn}
    call SHROUD_C_LOC(value, addr)
    call c_dataview_set_external_data_ptr_only(view%voidptr, addr)
!    call c_dataview_apply_type_shape(rv%voidptr, type, {rank}, extents)
end subroutine dataview_set_array_data_ptr_{typename}{nd}""".format(
        extents_decl=extents_decl,
        extents_asgn=extents_asgn, **d)

def print_get_data(d):
    # typename - part of function name
    # nd       - number of dimensions
    # f_type   - fortran type
    # shape     - :,:, to match nd
    if d['rank'] == 0:
        return """
! Generated by genfsidresplicer.py
subroutine dataview_get_data_{typename}{nd}{suffix}(view, value)
    use iso_c_binding
    implicit none
    class(dataview), intent(IN) :: view
    {f_type}, pointer, intent(OUT) :: value{shape}
    type(C_PTR) cptr

    cptr = view%get_void_ptr()
    if (c_associated(cptr)) then
      call c_f_pointer(cptr, value)
    else
      nullify(value)
    endif
end subroutine dataview_get_data_{typename}{nd}{suffix}""".format(**d)

    else:
        return """
! Generated by genfsidresplicer.py
subroutine dataview_get_data_{typename}{nd}{suffix}(view, value)
    use iso_c_binding
    implicit none
    class(dataview), intent(IN) :: view
    {f_type}, pointer, intent(OUT) :: value{shape}
    type(C_PTR) cptr
    integer rank
    integer(SIDRE_LENGTH) extents({rank})

    cptr = view%get_void_ptr()
    if (c_associated(cptr)) then
      rank = view%get_shape({rank}, extents)
      call c_f_pointer(cptr, value, extents)
    else
      nullify(value)
    endif
end subroutine dataview_get_data_{typename}{nd}{suffix}""".format(**d)


class AddMethods(object):
    """Create lines necessary to add generic methods to a derived type.
    Loops over types and rank.

    procedure :: {stem}_{typename}{nd}{suffix} => {wrap_class}_{stem}_{typename}{nd}{suffix}
    generic :: {stem} => &
         gen1, &
         genn
    """
    def __init__(self, wrap_class):
        self.wrap_class = wrap_class
        self.lines = []
        self.methods = []

    @staticmethod
    def type_bound_procedure_part(d):
        return 'procedure :: {stem}_{typename}{nd}{suffix} => {wrap_class}_{stem}_{typename}{nd}{suffix}'.format(**d)

    @staticmethod
    def type_bound_procedure_generic(d):
        return '{stem}_{typename}{nd}{suffix}'.format(**d)

    def add_method(self, stem, fcn, scalar=False, **kwargs):
        self.methods.append((stem, fcn, scalar, kwargs))

    def gen_type_bound(self):
        lines = []
        for stem, fcn, scalar, kwargs in self.methods:
            generics = []
            extra = dict(
                wrap_class=self.wrap_class,
                stem=stem,
                )
            extra.update(kwargs)
            foreach_type(lines, AddMethods.type_bound_procedure_part, scalar=scalar, **extra)
            foreach_type(generics, AddMethods.type_bound_procedure_generic, scalar=scalar, **extra)
            
            lines.append('generic :: {stem} => &'.format(stem=stem))
            for gen in generics[:-1]:
                lines.append('    ' + gen + ',  &')
            lines.append('    ' + generics[-1])
        return lines

    def gen_body(self):
        lines = []
        for stem, fcn, scalar, kwargs in self.methods:
            foreach_type(lines, fcn, scalar=scalar, **kwargs)
        return lines


def foreach_type(lines, fcn, scalar=False, **kwargs):
    """ Call fcn once for each type and rank, appending to lines.
    kwargs - additional values for format dictionary.
    """
    shape = []
    lbound = []
    for nd in range(maxdims + 1):
        shape.append(':')
        lbound.append('lbound(value,%d)' % (nd+1))
    d = dict(
        suffix=''     # suffix of function name
    )
    d.update(kwargs)
    indx = 0
    for typetuple in types:
        d['typename'], d['f_type'], d['sidre_type'] = typetuple

        # scalar values
        # XXX - generic does not distinguish between pointer and non-pointer
#        d['rank'] = -1
#        d['nd'] = 'scalar'
#        d['shape'] = ''
#        lines.append(fcn(d))

        # scalar pointers
        d['index'] = indx
        indx += 1
        d['rank'] = 0
        d['shape'] = ''
        d['lower_bound'] = ''

        if scalar:
            d['nd'] = ''
            lines.append(fcn(d))
        else:
            d['nd'] = '_scalar'
            lines.append(fcn(d))
            for nd in range(1,maxdims+1):
                d['index'] = indx
                indx += 1
                d['rank'] = nd
                d['nd'] = '_%dd' % nd
                d['shape'] = '(' + ','.join(shape[:nd]) + ')'
                d['lower_bound'] = '(' + ','.join(lbound[:nd]) + ')'
                lines.append(fcn(d))

#----------------------------------------------------------------------

def group_string():
    """Text for functions with get and set strings for a group.

    get_string  =>   grp->getView(name)->getString()
    set_string  =>   grp->getView(name)->setString()
    """
    return """
subroutine datagroup_get_string(group, name, value)
    use iso_c_binding
    class(datagroup), intent(IN) :: group
    character(*), intent(IN) :: name
    character(*), intent(OUT) :: value
    integer(C_INT) :: lname
    type(C_PTR) view

    lname = len_trim(name)
    view = c_datagroup_get_view_from_name_bufferify(group%voidptr, name, lname)
    call c_dataview_get_string_bufferify(view, value, len(value, kind=C_INT))
end subroutine datagroup_get_string

subroutine datagroup_set_string(group, name, value)
    use iso_c_binding
    class(datagroup), intent(IN) :: group
    character(*), intent(IN) :: name
    character(*), intent(IN) :: value
    integer(C_INT) :: lname
    type(C_PTR) view

    lname = len_trim(name)
    view = c_datagroup_get_view_from_name_bufferify(group%voidptr, name, lname)
    call c_dataview_set_string_bufferify(view, value, len_trim(value, kind=C_INT))
end subroutine datagroup_set_string
"""


#----------------------------------------------------------------------

def gen_fortran():
    """Generate splicers used by Shroud.
    """
    print('! Generated by genfsidresplicer.py')

    # DataGroup
    t = AddMethods('datagroup')
    t.add_method('get_scalar', group_get_scalar, True)
    t.add_method('set_scalar', group_set_scalar, True)
    t.add_method('create_array_view', group_create_array_view)
    t.add_method('set_array_data_ptr', group_set_array_data_ptr)

    print('! splicer begin class.DataGroup.type_bound_procedure_part')
    for line in t.gen_type_bound():
        print(line)
    print('procedure :: get_string => datagroup_get_string')
    print('procedure :: set_string => datagroup_set_string')
    print('! splicer end class.DataGroup.type_bound_procedure_part')

    print()
    print('------------------------------------------------------------')
    print()

    print('! splicer begin class.DataGroup.additional_functions')
    for line in t.gen_body():
        print(line)
    print(group_string())
    print('! splicer end class.DataGroup.additional_functions')


    # DataView
    t = AddMethods('dataview')
    t.add_method('get_data', print_get_data, suffix='_ptr')
    t.add_method('set_array_data_ptr', view_set_array_data_ptr)

    print('! splicer begin class.DataView.type_bound_procedure_part')
    for line in t.gen_type_bound():
        print(line)
    print('! splicer end class.DataView.type_bound_procedure_part')

    print()
    print('------------------------------------------------------------')
    print()

    print('! splicer begin class.DataView.additional_functions')
    for line in t.gen_body():
        print(line)
    print('! splicer end class.DataView.additional_functions')

######################################################################

if __name__ == '__main__':
    try:
        cmd = sys.argv[1]
    except IndexError:
        raise RuntimeError("Missing command line argument")

    if cmd == 'fortran':
        # fortran splicers
        gen_fortran()
    elif cmd == 'test':
        AllocateAllocatable(print)
    else:
        raise RuntimeError("Unknown command")
