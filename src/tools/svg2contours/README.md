## Running the svg2contours Python script

The `svg2contours` script converts [SVG](https://developer.mozilla.org/en-US/docs/Web/SVG) images to [MFEM NURBS meshes](https://mfem.org/mesh-format-v1.0/#nurbs-meshes) using the [svgpathtools](https://github.com/mathandy/svgpathtools) library. 

The latter can be used with axom's `quest_winding_number` example application to 
sample the winding number field over the generated curves.

Full SVG support requires a (slightly) patched copy of [svgpathtools@1.6.1](https://github.com/mathandy/svgpathtools/releases/tag/v1.6.1), as described in this document.

### Create a virtual environment

```shell
> python3 -m venv venv

# linux
> source venv/bin/activate
# windows bash
> source ./venv/Scripts/activate

> pip3 install -r requirements.txt
```

### Apply patches to svgpathtools for proper treatment of rotated ellipses and rounded rectangles

[svgpathtools@1.6.1](https://github.com/mathandy/svgpathtools/releases/tag/v1.6.1) has a bug in applying the correct rotation angle when transforming ellipses and elliptical arcs. 
This can be resolved by applying the following patch:
```shell
> patch  -p1 venv/lib/python3.9/site-packages/svgpathtools/path.py -i svgpathtools-1.6.1-ellipse-rotation.patch --verbose 
```
See: https://github.com/mathandy/svgpathtools/pull/221

It has another bug related to rounded rectangles, which can be resolved by applying the following pathc:
```shell
> patch  -p1 venv/lib/python3.9/site-packages/svgpathtools/svg_to_paths.py  -i svgpathtools-1.6.1-rounded-rect.patch --verbose
```
See: https://github.com/mathandy/svgpathtools/pull/222

#### Developer's note:
These patches were generated from a git commit in each of the above pull requests using the following commands:
```shell
 > git format-patch -1 260a44ed2c0d114f77d57016d6d143a50729aca9 --stdout > svgpathtools-1.6.1-ellipse-rotation.patch
 > git format-patch -1 ec1e1101037fcd66967caa40bc2b038c928bae4f --stdout > svgpathtools-1.6.1-rounded-rect.patch
```

### Run the script on an input SVG mesh

To convert an SVG to an mfem NURBS mesh, run the following command:
```shell
> cd <axom_root>/<build_dir>
> ../src/tools/svg2contours/svg2contours.py -i ../data/contours/svg/shapes.svg 

SVG dimensions: width='210mm' height='297mm' viewBox='0 0 210 297'
Wrote 'drawing.mesh' with 54 vertices and NURBS 27 elements
```
Note: This assumes your axom clone has the `data` submodule located
at `<axom_root>/data`

### Run the quest winding number example
Now that we have an MFEM NURBS mesh, we can run our winding number application

```shell
> cd <axom_root>/<build_dir>/
> ./examples/quest_winding_number_ex      \
    -i ./drawing.mesh  \
    query_mesh --min 0 0 --max 250 250 --res 500 500 
```
