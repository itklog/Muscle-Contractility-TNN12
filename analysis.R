Step 1. Set up your R project
install.packages(c("rentrez", "ggplot2"))
library(rentrez)
library(ggplot2)

Step 2. Retrieve your human NCBI Gene record using R
gene_symbol <- "TNNI2"          
transcript_accession <- "NM_003282.4"
clinvar_id <- "12436"
variant_protein <- "p.Arg156Ter" 

for (x in c("data", "results", "figures")) {
  dir.create(x, showWarnings = FALSE)
}
hits <- entrez_search(
  db = "gene", 
  term = paste0(gene_symbol, "[Gene Name] AND txid9606[Organism]"), 
  retmax = 20
)

if (length(hits$ids) == 0) stop("Gene not found")

gene_id <- NA_character_
gene_info <- NULL

for (id in hits$ids) {
  rec <- entrez_summary(db = "gene", id = id)
  if (toupper(rec$name) == toupper(gene_symbol)) {
    gene_id <- id
    gene_info <- rec
    break
  }
}

if (is.na(gene_id)) stop("Exact human gene not found")

print(gene_info)

writeLines(
  c(
    paste("Gene symbol:", gene_info$name), 
    paste("NCBI Gene ID:", gene_id), 
    paste("Official name:", gene_info$description), 
    paste("Chromosome:", gene_info$chromosome)
  ), 
  "results/gene_summary.txt"
)

Step 3. Retrieve the matching protein sequence from the assigned transcript
gb <- entrez_fetch(
  db = "nuccore", 
  id = transcript_accession, 
  rettype = "gb", 
  retmode = "text"
)

pat <- '/protein_id="NP_[0-9]+[.][0-9]+"'
fragment <- regmatches(gb, regexpr(pat, gb))

if (length(fragment) == 0) {
  stop("NP protein ID not found; inspect the transcript record")
}

protein_accession <- gsub('^/protein_id="|"$', '', fragment)

protein_fasta <- entrez_fetch(
  db = "protein", 
  id = protein_accession, 
  rettype = "fasta", 
  retmode = "text"
)

if (!startsWith(protein_fasta, ">")) {
  stop("Invalid FASTA format retrieved")
}
directory
fasta_file <- file.path("data", paste0(gene_symbol, "_protein.fasta"))
writeLines(protein_fasta, fasta_file)

cat("Matched protein:", protein_accession, "\n")
writeLines(paste("NP accession:", protein_accession), "results/protein_accession.txt")


Step 4. Analyze protein length and amino-acid composition in R

fasta_lines <- readLines(fasta_file)
sequence <- paste(fasta_lines[!grepl("^>", fasta_lines)], collapse = "")
aa <- strsplit(toupper(sequence), "")[[1]]
protein_length <- length(aa)

standard <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
counts <- table(factor(aa, levels = standard))
results <- data.frame(AminoAcid = standard, Count = as.integer(counts))

write.csv(results, "results/amino_acid_composition.csv", row.names = FALSE)

# Generate Amino Acid Composition Bar Chart
fig1 <- ggplot(results, aes(x = AminoAcid, y = Count)) + 
  geom_col(fill = "grey45") + 
  theme_bw() + 
  labs(title = paste(gene_symbol, "Amino-acid composition"), x = "Amino acid", y = "Number of residues")

ggsave("figures/01_amino_acid_composition.png", fig1, width = 8, height = 5, dpi = 300)

Step 6. Produce Option A: Linear protein map showing total reference-length and affected residue

if (nzchar(variant_protein)) {
  # Extract the numeric amino-acid position from the HGVS string (e.g., "156" from "p.Arg156Ter")
  matched <- regmatches(variant_protein, regexpr("[0-9]+", variant_protein))
  position <- as.integer(matched)
  
  if (!is.na(position) && position >= 1 && position <= protein_length) {
    
    # Build the linear protein map plot
    fig2 <- ggplot() +
      geom_segment(aes(x = 1, xend = protein_length, y = 0, yend = 0), linewidth = 2, color = "grey40") +
      geom_point(aes(x = position, y = 0), size = 4, color = "firebrick") +
      annotate("text", x = position, y = 0.25, label = variant_protein, size = 4, fontface = "bold") +
      xlim(0, protein_length + 10) + 
      ylim(-0.6, 0.6) +
      theme_void() +
      labs(
        title = paste(gene_symbol, "Reference Protein and Variant Position"), 
        x = "Amino-acid position"
      ) +
      theme(
        axis.title.x = element_text(margin = margin(t = 10), size = 11, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold", margin = margin(b = 10))
      )
    
    print(fig2)
    
    # Save the plot as a high-resolution PNG file
    ggsave("figures/02_variant_position.png", fig2, width = 8, height = 2.5, dpi = 300)
    cat("Option A plot successfully generated and saved to figures/02_variant_position.png\n")
    
  } else {
    stop("Variant position is out of bounds or invalid based on protein length.")
  }
} else {
  stop("variant_protein is empty; Option A requires a defined protein change.")
}
