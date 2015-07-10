import sidre
ds = sidre.DataStore()

print ds
print dir(ds)

grp = ds.getRoot()
print "Group", grp
#print "Dir of Group", dir(grp)

print "Name of Group", grp.getName()



#xx = sidre.DataGroup()
#print xx

