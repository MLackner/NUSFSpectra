using Revise, Pkg

Pkg.activate("/Users/lackner/Documents/Projects/Northwestern/")

using NUSFSpectra

folder = "/Users/lackner/Documents/Projects/Northwestern/Michael/2018-11-15/"

spectra = batchprocess(folder)
