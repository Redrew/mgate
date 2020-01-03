# This script creates a global list `MGATE_AA` describing the different amino acids and their properties.
# Used to parse input data and plot heatmaps.
amino_acid <- function(code, term, name, property, order) {
  amino_acid <- code
  attributes(amino_acid) <- list(term=term, name=name, property=property, order=order)
  return(amino_acid)
}

MGATE_AA <- list(amino_acid("H", "His", "Histidine", "+", 1),
                 amino_acid("K", "Lys", "Lysine", "+", 2), 
                 amino_acid("R", "Arg", "Arginine", "+", 3),
                 amino_acid("D", "Asp", "Aspartic Acid", "-", 4),
                 amino_acid("E", "Glu", "Glutamic Acid", "-", 5),
                 amino_acid("N", "Asn", "Asparagine", "Polar/neutral", 6),
                 amino_acid("Q", "Gln", "Glutamine", "Polar/neutral", 7),
                 amino_acid("S", "Ser", "Serine", "Polar/neutral", 8),
                 amino_acid("T", "Thr", "Threonine", "Polar/neutral", 9),
                 amino_acid("C", "Cys", "Cysteine", "Polar/neutral", 10),
                 amino_acid("M", "Met", "Methionine", "Polar/neutral", 11),
                 amino_acid("A", "Ala", "Alanine", "Non-polar", 12),
                 amino_acid("I", "Ile", "Isoleucine", "Non-polar", 13),
                 amino_acid("L", "Leu", "Leucine", "Non-polar", 14),
                 amino_acid("V", "Val", "Valine", "Non-polar", 15),
                 amino_acid("F", "Phe", "Phenylalanine", "Non-polar", 16),
                 amino_acid("W", "Trp", "Tryptophan", "Non-polar", 17),
                 amino_acid("Y", "Tyr", "Tyrosine", "Non-polar", 18),
                 amino_acid("G", "Gly", "Glycine", "Unique", 19),
                 amino_acid("P", "Pro", "Proline", "Unique", 20),
                 amino_acid("%", "Ter", "Termination", "Terminated", 21))
                 
names(MGATE_AA) <- MGATE_AA
property_levels = sapply(MGATE_AA, function(aa) {attr(aa, "property")}) %>% unique()
for (aa in MGATE_AA) {
  attr(MGATE_AA[[aa]], "property") <- factor(attr(aa, "property"), levels=property_levels)
}