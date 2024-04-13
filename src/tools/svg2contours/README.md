### Create a virtual environment

```shell
> python3 -m venv venv

# linux
> source venv/bin/activate
# windows bash
> source ./venv/Scripts/activate

> pip3 install -r requirements.txt
```

### Open jupyter notebook in vs code

Need to run the notebook using a kernel from this venv.

#### This might work

* Open the ipynb file in vs code with jupyter extension
* Point jupyter to the ipython3 interpretter from your environment
```
CTRL+SHIFT+P -> "python select interpretter"
<select path to venv/bin/ipython3>
```
* Note: Might need to restart vs code to get this working

< based on https://stackoverflow.com/a/74410917 >

#### If not, try to create a kernel

```shell
> python3 -m ipykernel install --user --name=svg2contour
```
Then point jupyter to this kernel via "select kernel" button


### Run the notebook on your SVG file
This generates an mfem `*.mesh` file

### Run the quest winding number example

```shell
>./examples/quest_winding_number_ex -f ../data/contours/svg/drawing.mesh -v inline_mesh --min 0 0 --max 250 250 --res 500 500 
```
