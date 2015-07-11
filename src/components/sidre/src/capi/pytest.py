import sidre
ds = sidre.DataStore()

print ds
print dir(ds)

grp = ds.getRoot()
print "Group", grp
#print "Dir of Group", dir(grp)

print "Name of Group", grp.getName()

print "Name of non-existent view", grp.getViewName(100)



#xx = sidre.DataGroup()
#print xx

