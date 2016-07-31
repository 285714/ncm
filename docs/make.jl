using Documenter
include("../master/_core.jl")

makedocs()

run(`pandoc build/index.md --latex-engine=xelatex -o build/index.pdf`)

deploydocs(repo = "github.com/285714/ncm.git")
