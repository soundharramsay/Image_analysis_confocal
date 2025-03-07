# bioconvert_czi_2_seperate_tiff
bioconvert_czi_2_seperate_tiff

######################################  tiff conversions are good 

#!/bin/bash

# Input and output directories
input_dir="/Users/soundhar/Desktop/18_u20s_stable_halo_z8_stressor_feb18/na2sa_working"
output_dir="/Users/soundhar/Desktop/18_u20s_stable_halo_z8_stressor_feb18/na2sa_working/processed_tiff"

# Channel names (c1=DAIPI, c2=mCherry, c3=647)
channel_names=("DAIPI" "mCherry" "647")

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through all .czi files in the input directory
for file in "$input_dir"/*.czi; do
    # Get the base name of the file (without extension)
    base_name=$(basename "$file" .czi)

    # Run bfconvert to save as a multi-channel TIFF
    ./bfconvert "$file" "$output_dir/${base_name}_multichannel.tiff"

    # Use bfconvert to split the multi-channel TIFF into separate files
    for i in {1..3}; do
        output_file="$output_dir/${base_name}_${channel_names[$((i-1))]}.tiff"
        ./bfconvert -channel $((i-1)) "$output_dir/${base_name}_multichannel.tiff" "$output_file"
    done

    # Remove the temporary multi-channel TIFF file
    rm "$output_dir/${base_name}_multichannel.tiff"

    echo "Processed: $file"
done

echo "All files have been converted."



###################################################### macro script for color conversion
// Define the input and output directories
input_dir = "/Users/sr/Desktop/18_u20s_stable_halo_z8_stressor_feb18/sorbitol_working/output_tiff/";
output_dir = "/Users/sr/Desktop/18_u20s_stable_halo_z8_stressor_feb18/sorbitol_working/output_tiff/";

// Get a list of all .tiff files in the input directory
list = getFileList(input_dir);

// Loop through each file
for (i = 0; i < list.length; i++) {
    // Check if the file is a .tiff file and contains "_647" in its name
    if (endsWith(list[i], ".tif") && indexOf(list[i], "_647_magenta") != -1) {
        // Open the image
        open(input_dir + list[i]);

        // Check if the image is RGB
        if (bitDepth == 24) {
            // Convert RGB image to grayscale
            run("8-bit");
        }

        // Apply the "Cyan Hot" lookup table
        run("Magenta Hot");

        // Save the processed image with a new name
        new_name = replace(list[i], ".tiff", "_magenta.tif");
        saveAs("Tiff", output_dir + new_name);

        // Close the image
        close();
    }
}

// Print a message when done
print("Processing complete!");

##################################################### merging files DAPI and magenta 

// Define the input and output directories
input_dir = "/Users/sr/Desktop/18_u20s_stable_halo_z8_stressor_feb18/sorbitol_working/output_tiff/";
output_dir = "/Users/sr/Desktop/18_u20s_stable_halo_z8_stressor_feb18/sorbitol_working/output_tiff/";

// Get a list of all files in the input directory
list = getFileList(input_dir);

// Loop through each file
for (i = 0; i < list.length; i++) {
    // Check if the file is a "_647_magenta.tif" file
    if (endsWith(list[i], "_647_magenta.tif")) {
        // Extract the base name (e.g., "10_halo_cnt" from "10_halo_cnt_647_magenta.tif")
        base_name = replace(list[i], "_647_magenta.tif", "");

        // Construct the corresponding DAIPI file name
        daipi_file = base_name + "_DAIPI.tiff";

        // Check if the DAIPI file exists
        if (File.exists(input_dir + daipi_file)) {
            // Open the 647_magenta image
            open(input_dir + list[i]);

            // Open the DAIPI image
            open(input_dir + daipi_file);

            // Merge the two images
            run("Merge Channels...", "c3=" + daipi_file + " c6=" + list[i] + " create keep");

            // Flatten the merged image to create a single-page TIFF
            run("Flatten");

            // Save the flattened image
            merge_name = base_name + "_merge.tif";
            saveAs("Tiff", output_dir + merge_name);

            // Close all images
            close();
            close();
        } else {
            print("DAIPI file not found for: " + list[i]);
        }
    }
}

// Print a message when done
print("Processing complete!");
/
######################################## panel making in R
library(magick)
library(ggplot2)
library(cowplot)
library(grid)  # Load grid for rasterGrob()
library(gridExtra)

# Define file names (ensure they are strings)
file_names <- c(
  "1_halo_cnt_merge.tif",
  "1_halo_z8_cnt_merge.tif",
  "1_sor_1hr_halo_cnt_merge.tif",
  "1_sor_halo_z8_cnt_merge.tif",
  "1_sor_halo_z8_cn_removed_30minst_merge.tif",
  "2_halo_cnt_merge.tif",
  "2_halo_z8_cnt_merge.tif",
  "2_sor_1hr_halo_cnt_merge.tif",
  "2_sor_halo_z8_cnt_merge.tif",
  "2_sor_halo_z8_cn_removed_30minst_merge.tif",
  "3_halo_cnt_merge.tif",
  "3_halo_z8_cnt_merge.tif",
  "3_sor_1hr_halo_cnt_merge.tif",
  "3_sor_halo_z8_cnt_merge.tif",
  "3_sor_halo_z8_cn_removed_30minst_merge.tif",
  "4_halo_cnt_merge.tif",
  "4_halo_z8_cnt_merge.tif",
  "4_sor_1hr_halo_cnt_merge.tif",
  "4_sor_halo_z8_cnt_merge.tif",
  "4_sor_halo_z8_cn_removed_30minst_merge.tif",
  "5_halo_cnt_merge.tif",
  "5_halo_z8_cnt_merge.tif",
  "5_sor_1hr_halo_cnt_merge.tif",
  "5_sor_halo_z8_cnt_merge.tif",
  "5_sor_halo_z8_cn_removed_30minst_merge.tif"
)

# Directory containing input files
input_dir <- "./"  # Update this path if necessary

# Ensure the input directory exists
if (!dir.exists(input_dir)) {
  stop("Input directory does not exist:", input_dir)
}

# Create a list to store ggplot images
plot_list <- list()

# Read images and store them as ggplot objects
for (file_name in file_names) {
  input_path <- file.path(input_dir, file_name)
  
  # Check if file exists
  if (!file.exists(input_path)) {
    warning("File does not exist:", input_path)
    next
  }
  
  # Read image
  img <- image_read(input_path)
  
  # Convert image to raster for ggplot
  img_gg <- rasterGrob(as.raster(img), interpolate = TRUE)
  
  # Create a ggplot object with filename as title
  p <- ggplot() +
    annotation_custom(img_gg, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    ggtitle(file_name) +
    theme_void() +  # Remove background/grid
    theme(
      plot.title = element_text(
        hjust = 0.5,  # Center the title
        vjust = 0,    # Adjust vertical position (0 = bottom, 1 = top)
        size = 6,     # Set font size
        face = "bold", # Set font style (bold)
        color = "black" # Set font color
      ),
      plot.margin = margin(b = 20, unit = "pt")  # Add margin at the bottom for the title
    )
  
  # Store the plot
  plot_list[[file_name]] <- p
}

# Arrange images in a panel (5 columns)
panel <- plot_grid(plotlist = plot_list, ncol = 5, align = "hv")

# Show the final panel
#print(panel)

# Save the final panel as an SVG file
output_svg <- file.path(input_dir, "1_final_panel_1_5.svg")
ggsave(output_svg, panel, width = 8, height = 9, dpi = 300)

# Print a message when done
print(paste("Final panel saved as:", output_svg))
