from Bio import Phylo

tree = Phylo.read(snakemake.input.tree, "newick")
tree.root_at_midpoint()
Phylo.write(tree, snakemake.output.rooted, "newick")
