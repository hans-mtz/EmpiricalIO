# using Pkg; Pkg.add.(["Plots", "DSP"])

using Weave

filename = normpath(Weave.EXAMPLE_FOLDER, "FIR_design.jmd")
filename = "/Volumes/SSD Hans/Github/Emp IO/FIR_design.jmd"
path = "/Volumes/SSD Hans/Github/Emp IO"
# Julia markdown to HTML
weave(filename; doctype = "md2html", out_path = :pwd)

# Julia markdown to PDF
weave(filename; doctype = "md2pdf", out_path = path)

# Julia markdown to Pandoc markdown
weave(filename; doctype = "pandoc", out_path = :pwd)
