
# ==========================================
# BIOINFORMATICS PROJECT
# Comparative DNA Sequence Analysis
# ==========================================

# Step 1: Load required packages

library(Biostrings)
library(ggplot2)


# Step 2: Create required project folders

dir.create("data", showWarnings = FALSE)
dir.create("results/graphs", recursive = TRUE, showWarnings = FALSE)
dir.create("report", showWarnings = FALSE)


# ==========================================
# Step 3: Download DNA Sequences from NCBI
# ==========================================

# E. coli

if (!file.exists("data/ecoli.fasta")) {
  download.file(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=J01859.1&rettype=fasta&retmode=text",
    "data/ecoli.fasta",
    mode = "wb"
  )
}

# Salmonella enterica

if (!file.exists("data/salmonella.fasta")) {
  download.file(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=MG869135.1&rettype=fasta&retmode=text",
    "data/salmonella.fasta",
    mode = "wb"
  )
}

# Klebsiella pneumoniae

if (!file.exists("data/klebsiella.fasta")) {
  download.file(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=X93214.1&rettype=fasta&retmode=text",
    "data/klebsiella.fasta",
    mode = "wb"
  )
}

# Check downloaded files

list.files("data")


# ==========================================
# Step 4: Import DNA Sequences
# ==========================================

ecoli <- readDNAStringSet("data/ecoli.fasta")

salmonella <- readDNAStringSet("data/salmonella.fasta")

klebsiella <- readDNAStringSet("data/klebsiella.fasta")


# Display imported sequences

ecoli
salmonella
klebsiella


# ==========================================
# Step 5: Calculate DNA Sequence Length
# ==========================================

ecoli_length <- as.numeric(width(ecoli))

salmonella_length <- as.numeric(width(salmonella))

klebsiella_length <- as.numeric(width(klebsiella))

# Display sequence lengths

ecoli_length
salmonella_length
klebsiella_length


# ==========================================
# Step 6: Nucleotide Composition Analysis
# ==========================================

# Calculate nucleotide frequencies

ecoli_counts <- letterFrequency(
  ecoli,
  letters = c("A", "T", "G", "C")
)

salmonella_counts <- letterFrequency(
  salmonella,
  letters = c("A", "T", "G", "C")
)

klebsiella_counts <- letterFrequency(
  klebsiella,
  letters = c("A", "T", "G", "C")
)

# Display nucleotide counts

ecoli_counts
salmonella_counts
klebsiella_counts


# ==========================================
# Step 7: Create Comparison Table
# ==========================================

dna_results <- data.frame(
  
  Species = c(
    "Escherichia coli",
    "Salmonella enterica",
    "Klebsiella pneumoniae"
  ),
  
  Sequence_Length = c(
    ecoli_length,
    salmonella_length,
    klebsiella_length
  ),
  
  A = c(
    ecoli_counts[1, "A"],
    salmonella_counts[1, "A"],
    klebsiella_counts[1, "A"]
  ),
  
  T = c(
    ecoli_counts[1, "T"],
    salmonella_counts[1, "T"],
    klebsiella_counts[1, "T"]
  ),
  
  G = c(
    ecoli_counts[1, "G"],
    salmonella_counts[1, "G"],
    klebsiella_counts[1, "G"]
  ),
  
  C = c(
    ecoli_counts[1, "C"],
    salmonella_counts[1, "C"],
    klebsiella_counts[1, "C"]
  )
)


# ==========================================
# Step 8: Calculate GC Content
# ==========================================

# Count confirmed DNA bases (A, T, G, C)

confirmed_bases <- c(
  sum(ecoli_counts),
  sum(salmonella_counts),
  sum(klebsiella_counts)
)

# Calculate GC percentage
# Ambiguous bases are excluded from the denominator

gc_percentage <- round(
  c(
    (ecoli_counts[1, "G"] + ecoli_counts[1, "C"]) /
      confirmed_bases[1] * 100,
    
    (salmonella_counts[1, "G"] + salmonella_counts[1, "C"]) /
      confirmed_bases[2] * 100,
    
    (klebsiella_counts[1, "G"] + klebsiella_counts[1, "C"]) /
      confirmed_bases[3] * 100
  ),
  2
)

# Add GC percentage to results table

dna_results$GC_Percentage <- gc_percentage


# Calculate ambiguous bases

dna_results$Ambiguous_Bases <- c(
  ecoli_length - confirmed_bases[1],
  salmonella_length - confirmed_bases[2],
  klebsiella_length - confirmed_bases[3]
)


# Display final results

print(dna_results)


# Save results as CSV

write.csv(
  dna_results,
  "results/dna_analysis_results.csv",
  row.names = FALSE
)


# ==========================================
# Step 9: Data Visualization
# ==========================================


# Graph 1: GC Content Comparison

gc_plot <- ggplot(
  dna_results,
  aes(x = Species, y = GC_Percentage, fill = Species)
) +
  
  geom_col(width = 0.6) +
  
  geom_text(
    aes(label = paste0(GC_Percentage, "%")),
    vjust = -0.5,
    size = 4
  ) +
  
  labs(
    title = "GC Content Comparison",
    x = "Bacterial Species",
    y = "GC Content (%)"
  ) +
  
  theme_minimal() +
  
  theme(legend.position = "none")


# Display graph

print(gc_plot)


# Save graph

ggsave(
  "results/graphs/gc_content.png",
  plot = gc_plot,
  width = 8,
  height = 5
)



# Graph 2: Nucleotide Composition

nucleotide_data <- data.frame(
  
  Species = rep(dna_results$Species, 4),
  
  Nucleotide = rep(
    c("A", "T", "G", "C"),
    each = 3
  ),
  
  Count = c(
    dna_results$A,
    dna_results$T,
    dna_results$G,
    dna_results$C
  )
)


nucleotide_plot <- ggplot(
  nucleotide_data,
  aes(x = Species, y = Count, fill = Nucleotide)
) +
  
  geom_col(position = "dodge") +
  
  labs(
    title = "Nucleotide Composition of Bacterial DNA",
    x = "Bacterial Species",
    y = "Number of Nucleotides",
    fill = "Nucleotide"
  ) +
  
  theme_minimal()


# Display graph

print(nucleotide_plot)


# Save graph

ggsave(
  "results/graphs/nucleotide_composition.png",
  plot = nucleotide_plot,
  width = 9,
  height = 5
)



# Graph 3: DNA Sequence Length Comparison

length_plot <- ggplot(
  dna_results,
  aes(x = Species, y = Sequence_Length, fill = Species)
) +
  
  geom_col(width = 0.6) +
  
  geom_text(
    aes(label = Sequence_Length),
    vjust = -0.5,
    size = 4
  ) +
  
  labs(
    title = "DNA Sequence Length Comparison",
    x = "Bacterial Species",
    y = "Sequence Length (bp)"
  ) +
  
  theme_minimal() +
  
  theme(legend.position = "none")


# Display graph

print(length_plot)


# Save graph

ggsave(
  "results/graphs/sequence_length.png",
  plot = length_plot,
  width = 8,
  height = 5
)


# ==========================================
# Step 10: Biological Interpretation
# ==========================================

findings <- c(
  
  "BIOINFORMATICS PROJECT - KEY FINDINGS",
  
  "",
  
  "1. GC Content:",
  "Klebsiella pneumoniae showed the highest GC content (55.06%).",
  "Escherichia coli showed the lowest GC content (54.45%).",
  
  "",
  
  "2. DNA Sequence Length:",
  "Escherichia coli had the longest sequence (1541 bp).",
  "Klebsiella pneumoniae had the shortest sequence (1458 bp).",
  
  "",
  
  "3. Nucleotide Composition:",
  "Guanine (G) was the most frequent nucleotide in all three sequences.",
  
  "",
  
  "4. Ambiguous Bases:",
  "Klebsiella pneumoniae contained five ambiguous nucleotide positions.",
  "The other two sequences contained no ambiguous positions.",
  
  "",
  
  "5. Overall Observation:",
  "The three bacterial 16S rRNA sequences showed relatively similar GC percentages."
)


# Display findings

cat(findings, sep = "\n")


# Save findings

writeLines(
  findings,
  "report/findings.txt"
)


# ==========================================
# PROJECT COMPLETED
# ==========================================

print("DNA Sequence Analysis Completed Successfully!")
