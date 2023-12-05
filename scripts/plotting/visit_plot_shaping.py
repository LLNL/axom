# This is a script for VisIt that can plot a shaping output from Axom's
# shaping driver or quest_intersection_shaper_test. The script will save
# a tiled image of the various volume fractions, their sum, and the materials
# shown as a set of Isovolume'd plots.
#
# visit -v 3.3 -cli -nowin -s visit_plot_shaping.py [filename]
#
import sys
matcolors = {"a":(100,100,100,255),
             "b":(240,240,240,255),
             "mat0":(255,0,0,255),
             "mat1":(0,255,0,255),
             "mat2":(80,80,220,255),
             "mat3":(255,255,0,255),
             "mat4":(180,0,220,255),
             "mat5":(255,153,0,255),
             "mat6":(0,0,150,255)
            }

# Plot everything as Isovolume'd plots.
def PlotMaterials(mats, fmt):
    for matname in mats:
        if fmt == "MFEM":
            AddPlot("Subset", "main", 0, 0)
        else:
            AddPlot("Subset", "shaping_mesh", 0, 0)
        atts = SubsetAttributes(1)
        atts.colorType = atts.ColorBySingleColor
        atts.singleColor = matcolors[matname]
        atts.legendFlag = 0
        SetPlotOptions(atts)
        AddOperator("Isovolume", 0)
        iso = IsovolumeAttributes(1)
        iso.lbound = 0.5
        if fmt == "MFEM":
            iso.variable = "vol_frac_"+matname
        else:
            iso.variable = "shaping_mesh/vol_frac_"+matname
        SetOperatorOptions(iso)
    DrawPlots()
    
# Plot everything as volume fractions.
def PlotVolumeFractions(mats, fmt):
    first = True
    dx = 0
    dy = 0
    expr = ""
    index = 0
    for matname in mats:
        row = int(index / 4)
        col = index % 4
        if fmt == "MFEM":
            var = "vol_frac_"+matname
        else:
            var = "shaping_mesh/vol_frac_"+matname
        AddPlot("Pseudocolor", var, 0, 0)
        if first:
            DrawPlots()
            Query("SpatialExtents")
            obj = GetQueryOutputObject()
            ext = obj["extents"]
            dx = 1.1 * ext[1]
            dy = 1.1 * -ext[3]
            first = False

            expr = "<" + var + ">"
        else:
            atts = PseudocolorAttributes()
            atts.legendFlag = 0
            SetPlotOptions(atts)
            #print(f"index={index}, row={row}, col={col}")
            AddOperator("Transform", 0)
            trans = TransformAttributes()
            trans.doTranslate = 1
            trans.translateX = dx * col
            trans.translateY = dy * row
            SetOperatorOptions(trans)

            expr = expr + "+<" + var + ">"
        index = index + 1

    # Plot the volume fraction sum. Filter it so it does not look speckly when
    # it's all really close to 1.
    DefineScalarExpression("sum", expr)
    DefineScalarExpression("filter_sum", "if(gt(sum,0.99999), 1., sum)")
    AddPlot("Pseudocolor", "filter_sum", 0, 0)
    atts = PseudocolorAttributes()
    atts.legendFlag = 0
    SetPlotOptions(atts)
    AddOperator("Transform", 0)
    row = int(index / 4)
    col = index % 4
    trans = TransformAttributes()
    trans.doTranslate = 1
    trans.translateX = dx * col
    trans.translateY = dy * row
    SetOperatorOptions(trans)

    DrawPlots()

def main():
    SetWindowLayout(2)

    fmt = "Blueprint"
    db = "shaping.root"
    for a in Argv():
        if ".mfem_root" in a:
            db = a
            fmt = "MFEM"
            print(f"Using {db}")
            break
        elif ".root" in a:
            db = a
            fmt = "Blueprint"
            print(f"Using {db}")
            break

    OpenDatabase(db)
    md = GetMetaData(db)

    # Look at the scalars to get the materials. Note that we sort the
    # material names in the end because the order of scalars in the file
    # could differ depending on the order in which the VFs were added to
    # the dataset.
    mats = []
    for i in range(md.GetNumScalars()):
        name = md.GetScalars(i).name
        if fmt == "MFEM":
            key = "vol_frac_"
            pos = name.find(key)
            if pos == 0:
                matname=name[pos+len(key):]
                #print(f"name={name}, pos={pos}, matname={matname}")
                mats.append(matname)
        else:
            key = "/vol_frac_"
            pos = name.find(key)
            if pos != -1:
                matname=name[pos+len(key):]
                #print(f"name={name}, pos={pos}, matname={matname}")
                mats.append(matname)
    mats = sorted(mats)
    SetActiveWindow(1)
    PlotVolumeFractions(mats, fmt)

    SetActiveWindow(2)
    DeleteAllPlots()
    PlotMaterials(mats, fmt)

    SetActiveWindow(1)
    swa = GetSaveWindowAttributes()
    swa.saveTiled = 1
    swa.width = 2048
    SetSaveWindowAttributes(swa)
    SaveWindow()
    sys.exit()

if __name__ == "__main__":
    main()
