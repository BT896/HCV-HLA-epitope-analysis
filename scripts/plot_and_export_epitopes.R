# HCV HLA-associated epitope analysis
#
# Purpose:
#   1. Visualise epitope-hit presence across serial HCV isolates.
#   2. Export observed epitope amino-acid sequences to CSV.
#
# Inputs:
#   data/epitope_hits_with_coordinates_P102.txt
#   data/epitope_hits_with_coordinates_P110.txt
#
# Outputs:
#   results/P102_P110_HLA_epitope_map.png
#   results/HLA_variable_epitopes_P102_P110.csv
#
# Run from the repository root directory.

# ---- Libraries ----

required_packages <- c("ggplot2", "dplyr", "tidyr", "forcats", "patchwork", "stringr", "ggtext")

missing_packages <- setdiff(required_packages, rownames(installed.packages()))

if (length(missing_packages) > 0) {
  stop(
    "Please install these packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

dir.create("results", showWarnings = FALSE)

# ---- Parameters ----
target_serotypes <- c("HLA-A1", "HLA-A2","HLA-A24", "HLA-A26", "HLA-B15", "HLA-B27", "HLA-B40", "HLA-B44")

# host ??? alleles (for the coloured header)
host_alleles <- list("P102" = c("A1","A24","B15","B44"),
                     "P110" = c("A2", "A26", "B27", "B40")
                    )

# ---- Helpers ----
normalize_hla_classI_serotype <- function(x) {
  x <- as.character(x)
  na_idx <- is.na(x)
  x <- toupper(trimws(x))
  x <- gsub("\\s+", "", x)
  x <- ifelse(grepl("^HLA-", x), x, paste0("HLA-", x))
  x <- gsub("^HLA-CW", "HLA-C", x)                       # Cw -> C
  x <- sub("^HLA-([ABC])\\*0?([0-9]{2}).*$", "HLA-\\1\\2", x)  # A/B/C*02:01 -> HLA-A02
  x <- sub("^HLA-([ABC])0?([0-9]{2})$",      "HLA-\\1\\2", x)  # HLA-A02     -> HLA-A2
  x <- sub("^(HLA-[ABC])0([0-9])$",          "\\1\\2", x)      # HLA-A02     -> HLA-A2
  x[na_idx] <- NA_character_
  x
}
parse_ddmmyy <- function(code) as.Date(strptime(code, "%d%m%y"))
prefix_rank <- function(hla) ifelse(grepl("^HLA-A", hla), 1, ifelse(grepl("^HLA-B", hla), 2, 3))
num_part    <- function(hla) suppressWarnings(as.numeric(sub("^HLA-[ABC]([0-9]+)$", "\\1", hla)))

# ---- Genome regions (for the bar under x-axis) ----
color_map <- data.frame(
  start=c(1,192,384,747,810,1027,1658,1712,1973,2421),
  end  =c(191,383,746,809,1026,1657,1711,1972,2420,3011),
  region=c("core","E1","E2","p7","NS2","NS3","NS4A","NS4B","NS5A","NS5B")
)

# ---- Read one host ----
read_host <- function(host_id) {
  cand <- c(file.path("data", paste0("epitope_hits_with_coordinates_", host_id, ".txt")),
            file.path("data", paste0("epitope_hits_with_coordinates_", host_id))
            )
  cand <- cand[file.exists(cand)]
  if (!length(cand)) stop("Missing file for ", host_id,
                          " (epitope_hits_with_coordinates_", host_id, "[.txt]).")
  df <- read.delim(cand[1], stringsAsFactors = FALSE)
  
  req <- c("Coordinate1",
           "Coordinate2",
           "HLA_Association",
           "Fasta_Header",
           "Epitope_Hit")
  if (!all(req %in% names(df))) stop("File for ", host_id, " must contain: ", paste(req, collapse=", "))
  
  df <- df %>%
    filter(!is.na(Coordinate1), !is.na(Coordinate2),
           Coordinate1 >= 1, Coordinate2 <= 3011)
  
  # extract ddmmyy from header
  rx_start <- paste0("^(?i)", host_id, "[_-]([0-9]{6})")
  rx_any   <- paste0("(?i)\\b", host_id, "[_-]([0-9]{6})\\b")
  fh <- as.character(df$Fasta_Header)
  dc_start <- stringr::str_match(fh, rx_start)[,2]
  dc_any   <- stringr::str_match(fh, rx_any)[,2]
  DateCode <- ifelse(!is.na(dc_start), dc_start, dc_any)
  
  df <- df %>%
    mutate(
      Host = host_id,
      DateCode = DateCode,
      HLA_Serotype = normalize_hla_classI_serotype(HLA_Association)
    ) %>%
    filter(!is.na(DateCode),
           !is.na(HLA_Serotype),
           HLA_Serotype %in% target_serotypes)
  
  if (!nrow(df)) stop("No rows for requested serotypes in ", host_id, ".")
  
  present_dates <- sort(unique(df$DateCode))
  order_df <- data.frame(DateCode = present_dates,
                         Date = parse_ddmmyy(present_dates),
                         stringsAsFactors = FALSE) %>%
    arrange(Date) %>%
    mutate(
      DateLabel   = format(Date, "%d/%m/%Y"),               # for y-axis
      SampleLabel = paste0(host_id, " - ", DateLabel)       # Host - Date (for labels & grouping)
    )
  
  df %>%
    left_join(order_df, by = "DateCode") %>%
    mutate(
      SampleLabel = factor(SampleLabel, levels = order_df$SampleLabel),
      DateLabel   = factor(DateLabel,   levels = order_df$DateLabel)
    )
}

# ---- Read & combine ----
P102 <- read_host("P102")
P110 <- read_host("P110")
both <- bind_rows(P102, P110)

# Interleave by actual date (across hosts)
sample_order_df <- both %>% distinct(SampleLabel, Host, Date, DateLabel) %>% arrange(Date, Host)
label_levels <- sample_order_df$SampleLabel  # Host - Date, chronological

both <- both %>%
  mutate(
    SampleLabel = factor(as.character(SampleLabel), levels = label_levels),
    DateLabel   = factor(as.character(DateLabel),   levels = unique(sample_order_df$DateLabel))
  )

# HLA order (A???B then numeric), only requested serotypes
hla_levels <- both %>%
  distinct(HLA_Serotype) %>%
  arrange(prefix_rank(HLA_Serotype), num_part(HLA_Serotype)) %>%
  pull(HLA_Serotype)

# Build y factor: HLA block ? (Host - Date) for ordering
y_levels <- as.vector(unlist(lapply(hla_levels, function(h) paste0(h, " - ", label_levels))))
both <- both %>%
  mutate(
    HLA_Serotype = factor(HLA_Serotype, levels = hla_levels),
    YComb = factor(paste0(HLA_Serotype, " - ", SampleLabel), levels = y_levels)
  )

# ---- Grey out epitopes present in ALL isolates (constant across series) ----
total_isolates <- length(unique(both$SampleLabel))  # across both hosts & all timepoints

presence_by_epitope <- both %>%
  distinct(HLA_Serotype, Coordinate1, Coordinate2, SampleLabel) %>%  # de-dup per isolate
  count(HLA_Serotype, Coordinate1, Coordinate2, name = "n_present")

both <- both %>%
  left_join(presence_by_epitope,
            by = c("HLA_Serotype","Coordinate1","Coordinate2")) %>%
  mutate(
    ColorKey = ifelse(n_present == total_isolates, "Constant", as.character(Host))
  )

# ---- Colours ----
host_cols <- c("P102" = "#1f77b4",  # blue
               "P110" = "#d62728",  # red
               "Constant" = "#bdbdbd")  # grey for epitopes present in all isolates

# ---- Main plot with on-the-fly coloured y tick labels ----
p_main <- ggplot(both,
                 aes(x = Coordinate1, xend = Coordinate2,
                     y = YComb, yend = YComb, color = ColorKey)) +
  geom_segment(lineend = "butt", linewidth = 1) +
  scale_x_continuous(name = "HCV amino acid position",
                     limits = c(1, 3011),
                     breaks = seq(0, 3000, 200),
                     expand = c(0,0)) +
  # Compute coloured labels from the tick text x itself: "HLA - Host - dd/mm/yyyy"
  scale_y_discrete(
    name = "Class I HLA serotype - sample date",
    drop = TRUE,
    labels = function(x) {
      host <- sub("^.*-\\s*(P[0-9]+)\\s*-.*$", "\\1", x)
      date <- sub("^.*-\\s*P[0-9]+\\s*-\\s*", "", x)
      
      paste(host, "\u2212", date)
    }
  ) +
  scale_color_manual(values = host_cols) +
  guides(color = "none") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7),   
        axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(
          face = "bold",
          size = 14,
          margin = margin(r = 90),
          hjust = 0.5,
          vjust = 1),
    plot.title   = element_text(face = "bold", size = 18, hjust = 0.5, vjust = 4),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 10, b = 10, l = 110)
  ) +
  labs(title = "P102 & P110 HCV HLA-related mutations")

# ---- Build coloured y-axis tick labels (Host - Date, coloured by host) ----

# We want labels for every YComb level shown on the axis
yl_shown  <- levels(droplevels(both$YComb))

# Map each YComb to Host + Date (extract from SampleLabel)
yl_meta <- both %>%
  distinct(YComb, SampleLabel) %>%
  mutate(
    HostOnly = sub("\\s+-.*$", "", SampleLabel),
    DateOnly = sub("^.*-\\s*", "", SampleLabel)
  )

# ---- Dotted separators & single HLA labels in left margin ----
hla_shown <- sub(" - .*", "", yl_shown)
sep_positions <- which(hla_shown[-1] != hla_shown[-length(hla_shown)]) + 0.5
if (length(sep_positions)) {
  p_main <- p_main + geom_hline(yintercept = sep_positions, linetype = "dotted", color = "grey60")
}

# centers per HLA block
r <- rle(hla_shown)
block_sizes <- r$lengths
block_names <- r$values

# Centers as a fraction of the panel height in NPC coordinates (0..1 from bottom to top)
centers_npc <- (cumsum(block_sizes) - block_sizes/2) / length(yl_shown)

# Ensure only one coord_cartesian and enough left margin
p_main <- p_main +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 80))

# Draw HLA labels in the LEFT OUTER MARGIN using NPC coords
for (i in seq_along(block_names)) {
  p_main <- p_main +
    annotation_custom(
      grob = grid::textGrob(
        label = block_names[i],
        x = grid::unit(-0.06, "npc"),   # fixed left margin position (move more left with -0.08, etc.)
        y = grid::unit(centers_npc[i], "npc"),
        just = "right",
        gp = grid::gpar(fontface = "bold", cex = 0.9)
      )
    )
}

# ---- Genome bar beneath x-axis ----
p_genome <- ggplot(color_map, aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = region)) +
  geom_rect(color = "black", linewidth = 0.2) +
  geom_text(aes(x = (start + end)/2, y = -0.2, label = region),
            inherit.aes = FALSE, size = 3, vjust = 1, fontface = "bold") +
  scale_x_continuous(limits = c(1, 3011), breaks = seq(0,3000,200), expand = c(0,0)) +
  scale_y_continuous(limits = c(-0.5, 1), expand = c(0,0)) +
  scale_fill_manual(values = c(
    core = "#da7b45", E1 = "#00cbc2", E2 = "#18c38c", p7 = "#f87dbd",
    NS2 = "#99992c", NS3 = "#f2e87e", NS4A = "#67ac9c",
    NS4B = "#9ee9dd", NS5A = "#a2f6d1", NS5B = "#60a172"
  )) +
  theme_void() +
  theme(legend.position = "none")

# ---- Coloured header (replaces legend): host ??? alleles in matching colours ----
hdr_df <- data.frame(
  Host = c("P102","P110"),
  Text = c(paste0("P102 - ", paste(host_alleles$P102, collapse = ", ")),
           paste0("P110 - ", paste(host_alleles$P110, collapse = ", ")))
)
p_header <- ggplot(hdr_df, aes(x = 0, y = Host, label = Text, color = Host)) +
  geom_text(hjust = 0, size = 6, fontface = "bold") +
  scale_color_manual(values = c("P102"="#1f77b4","P110"="#d62728")) +
  xlim(c(0,1)) +
  theme_void() +
  theme(plot.margin = margin(1, 1, 1, 1)) +
  guides(color = "none")

# ---- Assemble: header on top, then main + genome bar ----
final_plot <- p_header / (p_main / p_genome + plot_layout(heights = c(5, 0.6))) +
  plot_layout(heights = c(0.6, 6))

ggsave(filename = "results/P102_P110_HLA_epitope_map.png",
       plot = final_plot,
       width = 24,
       height = 16,
       units = "in",
       dpi = 300,
       bg = "white")

# ---- CSV file output ----

outfile <- "results/HLA_variable_epitopes_P102_P110.csv"

# Identify variable epitope regions
variable_regions <- both %>%
                    filter(ColorKey != "Constant") %>%
                    distinct(HLA_Serotype, Coordinate1, Coordinate2)

# Retrieve all observed sequences for those regions
sequence_data <- both %>%
                 semi_join(variable_regions, by = c("HLA_Serotype", "Coordinate1", "Coordinate2")) %>%
                 mutate(Sample = paste0(Host, "_", DateCode),
                        Sequence = as.character(Epitope_Hit)) %>%
                 select(HLA_Serotype,
                        Coordinate1,
                        Coordinate2,
                        Sample,
                        Sequence) %>%
                 distinct()

# Check for multiple different sequences in the same sample at the same coordinates
conflicts <- sequence_data %>%
             count(
               HLA_Serotype, Coordinate1, Coordinate2, Sample,
               name = "n_sequences"
             ) %>%
             filter(n_sequences > 1)

if (nrow(conflicts) > 0) {
  print(conflicts)
  stop("Multiple epitope sequences found for the same sample and coordinates.")
}

# Convert sample rows into individual columns
wide_sequences <- sequence_data %>%
                  pivot_wider(id_cols = c(HLA_Serotype, Coordinate1, Coordinate2),
                              names_from = Sample,
                              values_from = Sequence)

# Add the reference epitope sequence (using the first observed sequence at each region)
reference_sequences <- sequence_data %>%
                       arrange(HLA_Serotype, Coordinate1, Coordinate2, Sample) %>%
                       group_by(HLA_Serotype, Coordinate1, Coordinate2) %>%
                       summarise(Epitope_Hit = first(Sequence),
                         .groups = "drop")

varying_regions <- reference_sequences %>%
                   left_join(wide_sequences,
                             by = c("HLA_Serotype", "Coordinate1", "Coordinate2")
                             ) %>%
                    arrange(HLA_Serotype, Coordinate1, Coordinate2)

# Export CSV
write.csv(varying_regions, outfile, row.names = FALSE, na = "")

cat("Wrote", nrow(varying_regions),
    "epitope regions with",
    ncol(varying_regions) - 4,
    "sample sequence columns to",
    normalizePath(outfile),
    "\n")