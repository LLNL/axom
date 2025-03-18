import os, sys

runs = {
  "build-rzansel-blueos_3_ppc64le_ib_p9-clang@10.0.1.2_cuda-release" : ["seq", "omp", "cuda"],
  "build-rzwhippet-toss_4_x86_64_ib-clang@14.0.6-release" : ["seq", "omp"],
  "build-rzwhippet-toss_4_x86_64_ib-gcc@10.3.1-release" : ["seq", "omp"],
  "build-rzwhippet-toss_4_x86_64_ib-intel@2022.1.0-release" : ["seq", "omp"],
  "build-rzvernal-toss_4_x86_64_ib_cray-rocmcc@6.2.1_hip-release" : ["seq", "hip"],
  "build-rzvernal-toss_4_x86_64_ib_cray-clang@17.0.0_hip-release" : ["seq", "hip"],
  "build-rzadams-toss_4_x86_64_ib_cray-cce@18.0.0_hip-release" :["seq", "hip"],
  "build-rzadams-toss_4_x86_64_ib_cray-rocmcc@6.2.1_hip-release" :["seq", "hip"]
}

sizes = (50, 100, 200, 500, 1000, 1500, 2000, 4000, 8000)

def generate(method):
  for r in runs:
    if not os.path.exists(r):
      print(f"Skipping {r}.")
      continue
    filename = os.path.join(r, "run_concentric_circles.bash")
    f = open(filename, "wt")
    f.write("#!/bin/bash\n\n")
    f.write("CONCENTRIC_CIRCLES=./examples/mir_concentric_circles\n\n")

    for s in sizes:
      f.write(f"# Size {s}\n")
      for policy in runs[r]:
        f.write(f'echo "Running --gridsize {s} --numcircles 5 --policy {policy} --method {method}"\n')
        f.write(f'$CONCENTRIC_CIRCLES --gridsize {s} --numcircles 5 --policy {policy} --method {method} > result_{policy}_{s}.txt\n\n')

    f.close()
    os.chmod(filename, 0o700)
    print(f"Wrote {filename}")

def read_timings(filename, searchKey):
  retval = "" # no data
  try:
    lines = open(filename, "rt").readlines()
    for line in lines:
      pos = line.find(searchKey)
      if pos != -1:
        toks = [x for x in line.split() if x != '']
        retval = float(toks[2]) # timings (I)
        break
  except:
    pass
  return retval

def make_csv(method, outputfile):
  # I was measuring runMIR before, which includes device transfers and HDF5 output.
  # I don't want that stuff in the timings because it is not relevant to timing MIR
  # and it can vary wildly.

  # Measure just the MIR algorithm
  searchKey = "EquiZAlgorithm"
  if method == "elvira":
    searchKey = "ElviraAlgorithm"

  columns = []
  # Add sizes column
  sc = ["Sizes"]
  for s in sizes:
    sc.append(s)
  columns.append(sc)

  # Gather data.
  for r in sorted(runs.keys()):
    if not os.path.exists(r):
      continue
    for policy in runs[r]:
      buildname = r[6:-8]
      name = f"{buildname} {policy.upper()}"
      data = [name]
      for s in sizes:
        filename = os.path.join(r, f"result_{policy}_{s}.txt")
        value = read_timings(filename, searchKey)
        data.append(value)
      columns.append(data)

  # Write data
  f = open(outputfile, "wt")
  nrows = len(columns[0])
  for i in range(nrows):
    rowdata = [str(c[i]) for c in columns]
    line = ",".join(rowdata)
    f.write(f"{line}\n")
  f.close()

def main():
  method = "equiz"
  if "--equiz" in sys.argv:
    method = "equiz"
  if "--elvira" in sys.argv:
    method = "elvira"

  if "--generate" in sys.argv:
    generate(method)
  else:
    make_csv(method, "concentric_circle_timings.csv")

main()
