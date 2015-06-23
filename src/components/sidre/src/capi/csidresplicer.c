! C code that will be inserted into code via shroud splicer blocks



// DataGroup

// splicer begin class.DataGroup.C_definition
extern const char * ATK_datagroup_get_group_name_with_error_check(const ATK_datagroup * self, ATK_IndexType idx);

// splicer end class.DataGroup.C_definition


// splicer begin class.DataGroup.additional_functions
// identical to ATK_datagroup_get_group_name except return a pointer to an empty string.
static const char * empty_string  = "";

const char * ATK_datagroup_get_group_name_with_error_check(const ATK_datagroup * self, ATK_IndexType idx)
{
    const DataGroup *selfobj = static_cast<const DataGroup *>(self);
    const std::string & rv = selfobj->getGroupName(idx);
    return isNameValid(rv) ? rv.c_str() : empty_string;
}
// splicer end class.DataGroup.additional_functions
