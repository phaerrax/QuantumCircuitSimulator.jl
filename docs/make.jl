using Documenter, DocumenterCitations
using QuantumCircuitSimulator

# doctest dependencies
using ITensors, ITensorMPS

bib = CitationBibliography(joinpath(@__DIR__, "src", "bibliography.bib"); style=:numeric)

makedocs(;
    modules=[QuantumCircuitSimulator],
    sitename="QuantumCircuitSimulator.jl",
    repo=Remotes.GitHub("phaerrax", "QuantumCircuitSimulator.jl"),
    checkdocs=:exported,
    authors="Davide Ferracin <davide.ferracin@icloud.com> and contributors",
    pages=[
        "Home" => "index.md",
        "Building and running circuits" => "circuits.md",
        "Gates" => "gates.md",
        "Sparse Pauli-Lindblad noise model" => "noise.md",
        "Pauli-string sampling" => "paulistrings.md",
        "Library reference" => "api.md",
    ],
    plugins=[bib],
    format=Documenter.HTML(;
        mathengine=Documenter.MathJax3(
            Dict(
                :tex => Dict(
                    :inlineMath => [["\$", "\$"], [raw"\(", raw"\)"]],
                    :tags => "ams",
                    :packages =>
                        ["base", "ams", "autoload", "configmacros", "mathtools"],
                    :macros => Dict(
                        :C => [raw"\mathbb{C}"],
                        :Im => [raw"\mathrm{Im}"],
                        :N => [raw"\mathbb{N}"],
                        :Num => [raw"\mathscr{N}"],
                        :R => [raw"\mathbb{R}"],
                        :Re => [raw"\mathrm{Re}"],
                        :abs => [raw"\lvert #1 \rvert", 1],
                        :adj => [raw"#1^\dagger", 1],
                        :avg => [raw"\langle #1\rangle", 1],
                        :blank => [raw"{-}"],
                        :bra => [raw"\langle #1 \rvert", 1],
                        :braket => [raw"\langle #1 \vert #2 \rangle", 2],
                        :conj => [raw"\bar{#1}", 1],
                        :defeq => [raw"\mathrel{\mathop:}="],
                        :diag => [raw"\operatorname{diag}"],
                        :dissipator => [raw"\mathscr{D}"],
                        :eu => [raw"\mathrm{e}"],
                        :fcomm => [raw"\{#1\}", 1],
                        :id => [raw"1"],
                        :imat => [raw"I_{#1}", 1, ""],
                        :innp => [raw"\langle #1, #2\rangle", 2],
                        :iu => [raw"\mathrm{i}"],
                        :ket => [raw"\lvert #1 \rangle", 1],
                        :lindblad => [raw"\mathscr{L}"],
                        :norm => [raw"\lVert #1 \rVert", 1],
                        :outp => [raw"\lvert #1\rangle \langle #2\rvert", 2],
                        :paulim => [raw"\sigma_-^{(#1)}", 1, ""],
                        :paulip => [raw"\sigma_+^{(#1)}", 1, ""],
                        :paulix => [raw"\sigma_x^{(#1)}", 1, ""],
                        :pauliy => [raw"\sigma_y^{(#1)}", 1, ""],
                        :pauliz => [raw"\sigma_z^{(#1)}", 1, ""],
                        :phantomadj => [raw"^{\vphantom{\dagger}}"],
                        :proj => [raw"\outp{#1}{#1}", 1],
                        :sb => [raw"_{#1}", 1],
                        :set => [raw"\{\, #1 \;\vert\; #2\,\}", 2],
                        :spin => [raw"S_{#1}", 1],
                        :spinup => [raw"\uparrow"],
                        :spindown => [raw"\downarrow"],
                        :tr => [raw"\operatorname{tr}"],
                        :transpose => [raw"#1^{\mathrm{T}}", 1],
                        :tsp => [raw"^{\otimes #1}", 1],
                    ),
                ),
            ),
        ),
    ),
)

deploydocs(; repo="github.com/phaerrax/QuantumCircuitSimulator.jl")
