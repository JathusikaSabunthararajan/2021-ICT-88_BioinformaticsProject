# Bioinformatics Project
# Bacterial DNA Sequence Analysis
# Project: Comparative DNA Sequence Analysis


# Create required project folders

dir.create("data", showWarnings = FALSE)
dir.create("results/graphs", recursive = TRUE, showWarnings = FALSE)
dir.create("report", showWarnings = FALSE)



# Step 5: Download DNA sequences from NCBI

# E. coli
download.file(
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=J01859.1&rettype=fasta&retmode=text",
  "data/ecoli.fasta",
  mode = "wb"
)

# Salmonella enterica
download.file(
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=MG869135.1&rettype=fasta&retmode=text",
  "data/salmonella.fasta",
  mode = "wb"
)

# Klebsiella pneumoniae
download.file(
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=X93214.1&rettype=fasta&retmode=text",
  "data/klebsiella.fasta",
  mode = "wb"
)

# Check downloaded files
list.files("data")


# ==========================================
# Step 6: Import DNA sequences
# ==========================================

# Load Biostrings package
library(Biostrings)

# Import E. coli sequence
ecoli <- readDNAStringSet("data/ecoli.fasta")

# Import Salmonella sequence
salmonella <- readDNAStringSet("data/salmonella.fasta")

# Import Klebsiella sequence
klebsiella <- readDNAStringSet("data/klebsiella.fasta")

# Display imported sequences
ecoli
salmonella
klebsiella

# Calculate DNA sequence lengths

width(ecoli)
width(salmonella)
width(klebsiella)


# ==========================================
# Step 7: Nucleotide Composition Analysis
# ==========================================

# Calculate nucleotide frequencies

ecoli_counts <- letterFrequency(ecoli,
                                letters = c("A", "T", "G", "C"))

salmonella_counts <- letterFrequency(salmonella,
                                     letters = c("A", "T", "G", "C"))

klebsiella_counts <- letterFrequency(klebsiella,
                                     letters = c("A", "T", "G", "C"))

# Display nucleotide counts

ecoli_counts
salmonella_counts
klebsiella_counts



# ==========================================
# Calculate GC Content
# ==========================================

ecoli_gc <- letterFrequency(ecoli, "GC", as.prob = TRUE) * 100

salmonella_gc <- letterFrequency(salmonella, "GC", as.prob = TRUE) * 100

klebsiella_gc <- letterFrequency(klebsiella, "GC", as.prob = TRUE) * 100

# Display GC percentages

ecoli_gc
salmonella_gc
klebsiella_gc



# ==========================================
# Create Comparison Table
# ==========================================

dna_results <- data.frame(
  
  Species = c("Escherichia coli",
              "Salmonella enterica",
              "Klebsiella pneumoniae"),
  
  Sequence_Length = c(
    as.numeric(width(ecoli)),
    as.numeric(width(salmonella)),
    as.numeric(width(klebsiella))
  ),
  
  A = c(ecoli_counts[1, "A"],
        salmonella_counts[1, "A"],
        klebsiella_counts[1, "A"]),
  
  T = c(ecoli_counts[1, "T"],
        salmonella_counts[1, "T"],
        klebsiella_counts[1, "T"]),
  
  G = c(ecoli_counts[1, "G"],
        salmonella_counts[1, "G"],
        klebsiella_counts[1, "G"]),
  
  C = c(ecoli_counts[1, "C"],
        salmonella_counts[1, "C"],
        klebsiella_counts[1, "C"]),
  
  GC_Percentage = c(
    as.numeric(ecoli_gc),
    as.numeric(salmonella_gc),
    as.numeric(klebsiella_gc)
  )
)

# Display final results

print(dna_results)

# Save results as CSV

write.csv(dna_results,
          "results/dna_analysis_results.csv",
          row.names = FALSE)


# Step 7.1: Correct GC content calculation

# Count confirmed DNA bases
confirmed_bases <- c(
  sum(ecoli_counts),
  sum(salmonella_counts),
  sum(klebsiella_counts)
)

# Calculate GC percentage using confirmed bases
dna_results$GC_Percentage <- round(
  c(
    (ecoli_counts[1,"G"] + ecoli_counts[1,"C"]) / confirmed_bases[1] * 100,
    (salmonella_counts[1,"G"] + salmonella_counts[1,"C"]) / confirmed_bases[2] * 100,
    (klebsiella_counts[1,"G"] + klebsiella_counts[1,"C"]) / confirmed_bases[3] * 100
  ), 2
)

# Calculate ambiguous bases
dna_results$Ambiguous_Bases <- c(
  width(ecoli) - confirmed_bases[1],
  width(salmonella) - confirmed_bases[2],
  width(klebsiella) - confirmed_bases[3]
)

# Display updated table
print(dna_results)

# Save updated CSV
write.csv(dna_results, "results/dna_analysis_results.csv",
          row.names = FALSE)




# Step 8: Data Visualization

library(ggplot2)

# Check package






# Graph 1: GC Content Comparison

gc_plot <- ggplot(dna_results,
                  aes(x = Species, y = GC_Percentage,
                      fill = Species)) +
  
  geom_col(width = 0.6) +
  
  geom_text(aes(label = paste0(GC_Percentage, "%")),
            vjust = -0.5, size = 4) +
  
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
ggsave("results/graphs/gc_content.png",
       plot = gc_plot,
       width = 8,
       height = 5)





# Graph 2: Nucleotide Composition

nucleotide_data <- data.frame(
  
  Species = rep(dna_results$Species, 4),
  
  Nucleotide = rep(c("A", "T", "G", "C"), each = 3),
  
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
ggsave("results/graphs/nucleotide_composition.png",
       plot = nucleotide_plot,
       width = 9,
       height = 5)





# Graph 3: DNA Sequence Length Comparison

length_plot <- ggplot(
  dna_results,
  aes(x = Species, y = Sequence_Length, fill = Species)
) +
  
  geom_col(width = 0.6) +
  
  geom_text(aes(label = Sequence_Length),
            vjust = -0.5, size = 4) +
  
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
ggsave("results/graphs/sequence_length.png",
       plot = length_plot,
       width = 8,
       height = 5)




# ==========================================
# Step 9: Biological Interpretation
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
writeLines(findings, "report/findings.txt")

