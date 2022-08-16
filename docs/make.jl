using Documenter
using TaylorFourierTransform

makedocs(
    modules     = [TaylorFourierTransform],
    format      = Documenter.HTML(mathengine = Documenter.MathJax()),
    sitename    = "TaylorFourierTransform.jl",
    authors     = "Tom Van Acker",
    pages       = [ "Home"              => "index.md",
                    "DSL Manual"        =>
                         ["DTFT"                       => "dtft.md",
                          "Utilities"                  => "util.md"],
                    "Developer"         =>
                         ["Developer's Guide"          => "guide.md",
                          "Mathematical Background"    => "math.md"]
                  ]
)

deploydocs(
     repo = "github.com/timmyfaraday/TFT.jl.git"
)