# Items a benchmark deliberately leaves out

One file per benchmark (`<name>.txt`), one item key per line, a tab, and why. These items are
dropped when the benchmark is read (`Benchmark.__init__`) and are not counted as failures by the
download scripts, so a benchmark's contents are defined here on purpose instead of by whichever
downloads happened to succeed. A pinned homology split (see `Benchmark.split_map`) is only valid
for exactly the data it was built from, so adding or removing a line here changes the dataset
and needs a new split.
