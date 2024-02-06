# CircumferenceChecker
This program contains various tools for checking graph properties related to paths and cycles. In particular, it can be used to determine the circumference of a graph, the length of a longest induced cycle and the length of a longest induced path. The length of a path is always the number of edges it contains.

The latest version of this program can be obtained from <https://github.com/JarneRenders/CircumferenceChecker>.

### Installation

This requires a working shell and `make`. On Windows an easy way to simulate this is by using Windows Subsystem for Linux (WSL).

- Compile using: 
  * `make` to create a binary for the 64 bit version
  * `make 128bit` to create a binary for the 128 bit version
  * `make 192bit` to create a binary for the 192 bit version
  * `make 256bit` to create a binary for the 256 bit version
  * `make all` to create all of the above

The 64 bit version can handle graphs up to 64 vertices, the 128 bit version up to 128 vertices, etc.
Lower bit versions are always faster than the higher bit ones, hence it is recommended to use the version which higher, but closest to the order of the graphs you want to inspect.


### Usage of circumferenceChecker

This helptext can be found by executing `./circumferenceChecker -h`.

Usage: `./circumferenceChecker [-cf#|-pf#] [-Cdo#] [-h]`

Count and or filter graphs depending on their circumference, induced cycles
or paths.

Graphs are read from stdin in graph6 format. Graphs are sent to stdout in
graph6 format. If the input graph had a graph6 header, so will the
output graph (if it passes through the filter).

```
    -c, --induced-cycle\n\
            count the longest induced cycle of each graph and print in a table.
    -C, --complement
            print all the graphs that would not have been sent to stdout and
            vice versa.
    -d, --difference
            count difference with order of the graph.
    -f#, --forbidden=#
            send all graphs to stdout that contain an induced path or cycle
            of length # depending on the presence of -c or -p. Length of path
            is the number of edges.
    -h, --help
            print help message
    -o#, --output=#
            send all graphs with value # in the table to stdout.
    -p, --induced-path
            count the longest induced path of each graph and print in a table.
```