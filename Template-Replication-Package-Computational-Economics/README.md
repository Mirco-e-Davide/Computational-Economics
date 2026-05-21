# CompEcon Project Template 

This is the replication project template. 

## Step Number 1

* Click top right on `use this template`, then `create a new repository`.
* Choose a suitable name for your replication project.
* Clone your new repository to your computer.
* Open this locally in VScode.

(you get a nice preview of this document in VScode by executing command `Markdown: Open preview` - type this into the command 😉 )

This page should contain concise information for how to

1. install your julia package
2. run the unit tests
3. run the replication output (there needs to be a single entry point that runs everything, ideally called `run_all()`)
4. gives instructions for how to compile and where to see your replication [report](report.qmd). See [https://quarto.org](https://quarto.org) for instructions for how to use quarto.

## First steps

1. Copy the original replication package into the folder `replication-package`
2. (optional) follow instructions in `replication-package` to replicate desired exhibits
3. set up your julia package in this repo. To do so, start the julia REPL in this location (alt-J alt-O, otherwise look in command palette)
4. Choose a name for your replication package, for example `XYZ.jl` (choose something better and replace `XYZ` throughout with your choice!)
5. Generate your julia package in this directory with
   
    ```julia
    using Pkg
    Pkg.generate("XYZ.jl")
    ```
6. In the VScode file browser on the left, investigate the newly created directory
7. Activate your package:
   
    ```julia
    julia> Pkg.activate("XYZ.jl")

    # load it
    julia> using XYZ

    # run function
    julia> XYZ.greet()
    Hello World!
    ```
8. Edit the code in `XYZ.jl/src/XYZ.jl`, maybe adding a word to the greet function. Save the file. Call the function again.
   
    ```julia
    julia> XYZ.greet()
    Hello World, Earthlings!
    ```


