
Generating figures
------------------

1. Introductory example: 
   - VUE mind-mapping software (http://vue.tufts.edu/) was used to draw
     the initial diagram sidre_datastore_example.vue.  This was then saved
     as a PDF.
   - Inkscape vector graphics program (https://inkscape.org/en/) was used
     to import the PDF and perform manual editing to produce the vector
     graphics version sidre_datastore_example.svg and raster version
     sidre_datastore_example.png.
2. Sidre View state table:
   - MS Powerpoint was used to generate the table (file: sidre-view-states.pptx)
   - Table was saved as a *.png file.
3. Sidre interaction with Conduit:
   - The example cartoon tiny_mesh.png was drawn by hand and scanned.
   - VUE was used to draw the initial diagram sidreconduit.vue, which was
     saved as a PDF.
   - Inkscape was then used to
     - import the PDF,
     - rearrange the logical grouping of the graphical elements into layers,
     - save raster versions of the graphic as successive layers were set to
       "visible":
       - ds.png,
       - cds.png,
       - cdscoords.png,
       - cdstopo.png,
       - cdsfields.png;      
     - save the overall vector graphic file sidreconduit.svg.
     - Isolated shapes were edited and saved in raster format for the
       legend figures:
       - hexagon.png
       - rectangle.png
       - roundrectangle.png
   - Running the Sidre example sidre_createdatastore_ex produced two
     data files: tinymesh.root and tinymesh.json.  VisIt was used to open
     tinymesh.root and produce the rendering, tiny_mesh_rendered.png.
