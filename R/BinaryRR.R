# Install required packages if not already installed
# install.packages(c("hexSticker", "ggplot2", "dplyr", "showtext"))

library(hexSticker)
library(ggplot2)
library(dplyr)
library(showtext)

# Add Google Fonts
font_add_google("Roboto", "roboto")
font_add_google("Roboto Condensed", "roboto_condensed")
font_add_google("Fira Code", "fira_code")
showtext_auto()

# Create output directory
if (!dir.exists("hex_stickers")) {
  dir.create("hex_stickers")
}

# Load bbssr package to use the actual BinaryRR function
# Make sure bbssr package is installed and loaded
if (!require(bbssr, quietly = TRUE)) {
  stop("bbssr package is required. Please install it first.")
}

# Design: Exact Rejection Region Heatmap
# Shows the precise rejection region for exact tests
create_bbssr_sticker <- function() {
  # Use actual BinaryRR function from bbssr package
  N1 <- 20
  N2 <- 20  # Make it square for better visualization
  alpha <- 0.025
  Test <- 'Boschloo'  # Use Boschloo test like in your example

  RR <- BinaryRR(N1, N2, alpha, Test)

  # Convert to data frame for ggplot
  # Note: RR matrix has dimensions (N1+1) x (N2+1)
  # Row index corresponds to X1 (0 to N1), Column index corresponds to X2 (0 to N2)
  # Need to flip Y-axis to match matrix indexing
  rr_data <- expand.grid(X1 = 0:N1, X2 = 0:N2)

  # Convert matrix to vector with correct indexing
  # Matrix RR[i,j] corresponds to X1=i-1, X2=j-1
  rr_data$reject <- as.vector(RR)

  # Flip Y-axis to match matrix visualization (high X1 values at bottom)
  rr_data$Y_flipped <- N1 - rr_data$X1

  # Set colors: bright green for TRUE (rejection region), very light gray for FALSE (acceptance region)
  rr_data$color <- ifelse(rr_data$reject, "#00CED1", "#DDD0FF")  # Bright green vs very light gray
  rr_data$alpha_val <- ifelse(rr_data$reject, 0.9, 0.8)  # Slightly more opaque for light gray

  rr_plot <- ggplot(rr_data, aes(x = X2, y = Y_flipped)) +
    geom_point(aes(color = color, alpha = alpha_val), size = 1) +  # geom_tile → geom_point
    scale_color_identity() +
    scale_alpha_identity() +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.position = "none"
    ) +
    coord_fixed() +
    xlim(-0.5, N2 + 0.5) +
    ylim(-0.5, N1 + 0.5)

  sticker(
    subplot = rr_plot,
    package = "bbssr",
    p_size = 40,  # Much larger font size
    p_color = "#FFFFFF",
    p_family = "fira_code",
    p_fontface = "bold",
    p_y = 1.5,  # Adjust position slightly for larger font
    s_x = 1,
    s_y = 0.85,  # Move plot down slightly to make room for larger font
    s_width = 1.6,  # Make plot wider
    s_height = 1.1,  # Make plot taller
    h_fill = "#2A0845",
    h_color = "#7209B7",
    h_size = 1.2,
    dpi = 600,
    filename = "hex_stickers/bbssr_sticker.png"
  )

  cat("bbssr sticker created: Exact rejection region heatmap\n")
  cat("- Green tiles: TRUE (rejection region) - Reject H0\n")
  cat("- Light red tiles: FALSE (acceptance region) - Accept H0\n")
  cat("- Light gray background for clean appearance\n")
  cat("- Visualizes exact discrete nature of binary endpoint tests\n")
}

# Create the bbssr sticker
cat("Creating bbssr package hex sticker using actual BinaryRR() function...\n\n")

create_bbssr_sticker()

cat("\nSticker completed! File saved as 'hex_stickers/bbssr_sticker.png'\n")

# Display design concept
cat("\n=== DESIGN CONCEPT ===\n")
cat("Visualizes the exact rejection region as computed by BinaryRR() function\n")
cat("- Each tile represents a possible (X1, X2) outcome combination\n")
cat("- Uses actual BinaryRR(N1=20, N2=20, alpha=0.025, Test='Boschloo')\n")
cat("- Green tiles (bright): TRUE → Reject H0 (rejection region)\n")
cat("- Very light gray tiles: FALSE → Accept H0 (acceptance region)\n")
cat("- High resolution: dpi = 600 for crisp, professional quality\n")
cat("- Should show rejection region concentrated in lower-left (high X1, low X2)\n")
cat("- Light gray background for clean, professional appearance\n")
cat("- Fira Code font for modern, technical look\n")
cat("- Emphasizes the discrete, exact nature of binary endpoint statistical tests\n")
cat("- Directly represents the core functionality of the bbssr package\n")

# Display file information
if (file.exists("hex_stickers/bbssr_sticker.png")) {
  cat("\nFile size:", file.size("hex_stickers/bbssr_sticker.png"), "bytes\n")
  cat("File path: hex_stickers/bbssr_sticker.png\n")
}
