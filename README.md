# bioconvert_czi_2_seperate_tiff
1###### bioconvert_czi_2_seperate_tiff
2######imagej macro recoloring 647 >>> magenta
3##### imagej merging DAPI+magenta
4##### panel making in R
5#### mergeing 647 and g3bp1 
6#### merging 647_dapi_g3bp1
7######color conversion with intensity change 


######################################  tiff conversions are good  ## Users/soundhar/Desktop/softwares/bftools
#/Users/soundhar/Desktop/softwares/bftools

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


################################ mergeing 647 and g3bp1 
// Define the input and output directories
input_dir = "/Users/soundhar/Desktop/20_u20s_stable_g3bp1/poly_I_C_5OONGML_4HRS/processed_tiff/";
output_dir = "/Users/soundhar/Desktop/20_u20s_stable_g3bp1/poly_I_C_5OONGML_4HRS/processed_tiff/";

// Get a list of all files in the input directory
list = getFileList(input_dir);

// Loop through each file
for (i = 0; i < list.length; i++) {
    // Check if the file is a "_647_magenta.tif" file
    if (endsWith(list[i], "_647_magenta.tif")) {
        // Extract the base name (e.g., "10_halo_cnt" from "10_halo_cnt_647_magenta.tif")
        base_name = replace(list[i], "_647_magenta.tif", "");

        // Construct the corresponding G3BP1 file name
        g3bp1_file = base_name + "_G3BP1.tiff";

        // Check if the G3BP1 file exists
        if (File.exists(input_dir + g3bp1_file)) {
            // Open the 647_magenta image
            open(input_dir + list[i]);

            // Open the G3BP1 image
            open(input_dir + g3bp1_file);

            // Merge the two images
            run("Merge Channels...", "c1=" + g3bp1_file + " c2=" + list[i] + " create keep");

            // Flatten the merged image to create a single-page TIFF
            run("Flatten");

            // Save the flattened image with the desired output name
            merge_name = base_name + "_merge_647_g3bp1.tif";
            saveAs("Tiff", output_dir + merge_name);

            // Close all images
            close();
            close();
        } else {
            print("G3BP1 file not found for: " + list[i]);
        }
    }
}

// Print a message when done
print("Processing complete!");

############################## merging 647_dapi_g3bp1
// Define the input and output directories
input_dir = "/Users/soundhar/Desktop/20_u20s_stable_g3bp1/poly_I_C_5OONGML_4HRS/processed_tiff/";
output_dir = "/Users/soundhar/Desktop/20_u20s_stable_g3bp1/poly_I_C_5OONGML_4HRS/processed_tiff/";

// Get a list of all files in the input directory
list = getFileList(input_dir);

// Loop through each file
for (i = 0; i < list.length; i++) {
    // Check if the file is a "_647_magenta.tif" file
    if (endsWith(list[i], "_647_magenta.tif")) {
        // Extract the base name (e.g., "10_halo_cnt" from "10_halo_cnt_647_magenta.tif")
        base_name = replace(list[i], "_647_magenta.tif", "");

        // Construct the corresponding DAIPI and G3BP1 file names
        daipi_file = base_name + "_DAIPI.tiff";
        g3bp1_file = base_name + "_G3BP1.tiff";

        // Check if the DAIPI and G3BP1 files exist
        if (File.exists(input_dir + daipi_file) && File.exists(input_dir + g3bp1_file)) {
            // Open the 647_magenta image
            open(input_dir + list[i]);

            // Open the DAIPI image
            open(input_dir + daipi_file);

            // Open the G3BP1 image
            open(input_dir + g3bp1_file);

            // Merge the three images
            run("Merge Channels...", "c1=" + daipi_file + " c2=" + list[i] + " c3=" + g3bp1_file + " create keep");

            // Flatten the merged image to create a single-page TIFF
            run("Flatten");

            // Save the flattened image with the desired output name
            merge_name = base_name + "_DAPI_647_g3bp1.tif";
            saveAs("Tiff", output_dir + merge_name);

            // Close all images
            close();
            close();
            close();
        } else {
            print("DAIPI or G3BP1 file not found for: " + list[i]);
        }
    }
}

// Print a message when done
print("Processing complete!");


################################################### color conversion with intensity change 


// Define the input directory containing the TIFF files
input_dir = "/Users/soundhar/Desktop/21_u20s_stable_high_low_sorbitol_hexandiol/z8_low_matched/processed_tiff/";

// List all TIFF files in the input directory
list = getFileList(input_dir);

// Loop through each file in the directory
for (i = 0; i < list.length; i++) {
    // Check if the file ends with "_647.tiff"
    if (endsWith(list[i], "_647.tiff")) {
        // Open the image
        open(input_dir + list[i]);

        // Get the image title (filename without extension)
        title = getTitle();
        new_title = replace(title, ".tiff", "_magenta.tiff");

        // Apply the "Magenta Hot" LUT
        run("Magenta Hot");

        // Set the brightness/contrast (optional)
        setMinAndMax(0, 160);
        run("Apply LUT");

        // Save the processed image in the same directory
        saveAs("Tiff", input_dir + new_title);

        // Close the image
        close();
    }
}

// Print a message when done
print("Processing complete. Processed images saved in: " + input_dir);

##################################################################################
##############################################################################################################  z-stack >>>> max projection >>>>tiff conversion 


##### z-stack to tiff 

// Define input and output directories
input_dir = "/Users/soundhar/Desktop/24_u20s_idr_muts/idr_mutants_sorbitol_stress/";
output_dir = "/Users/soundhar/Desktop/24_u20s_idr_muts/idr_mutants_sorbitol_stress/processed_tiff/";

// Create output directory if it doesn't exist
File.makeDirectory(output_dir);

// Get a list of all .czi files in the input directory
list = getFileList(input_dir);

// Loop through each file
for (i = 0; i < list.length; i++) {
    if (endsWith(list[i], ".czi")) {
        // Open the .czi file
        open(input_dir + list[i]);

        // Get the base name of the file (without extension)
        base_name = File.nameWithoutExtension;

        // Perform maximum intensity projection
        run("Z Project...", "projection=[Max Intensity]");

        // Apply the Magenta Hot LUT
        run("Magenta Hot");

        // Save the processed image
        saveAs("Tiff", output_dir + "MAX_" + base_name + ".tif");

        // Close the image
        close();

        print("Processed: " + list[i]);
    }
}

print("All files processed.");

############ 
#################
####################### multi-page tiff to single page seperate tiff

#!/bin/bash

# Input directory containing .tiff files
input_dir="./"

# Output directory for separated channels
output_dir="./"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through all .tiff files in the input directory
for file in "$input_dir"/*.tif; do
    # Get the base name of the file (without extension)
    base_name=$(basename "$file" .tif)

    # Split the channels using bfconvert
    /Users/soundhar/Desktop/softwares/bftools/bfconvert -channel 0 "$file" "$output_dir/${base_name}_dapi.tiff"
    /Users/soundhar/Desktop/softwares/bftools/bfconvert -channel 1 "$file" "$output_dir/${base_name}_mcherry.tiff"
    /Users/soundhar/Desktop/softwares/bftools/bfconvert -channel 2 "$file" "$output_dir/${base_name}_magenta.tiff"
    /Users/soundhar/Desktop/softwares/bftools/bfconvert -channel 3 "$file" "$output_dir/${base_name}_g3bp1.tiff"

    echo "Processed: $file"
done

echo "All files processed."



################################################################### if u have three tiff of DAPI and Magenta, but DAPI is in wrong color. This script change the color of DAPI to blue and create a merged DAPI and magenta 
// Set the input directory
inputDir = "/Users/soundhar/Desktop/24_u20s_idr_muts/z8_low_sorbitol_expo_recovery_time/processed/";

// Get all files in the directory
list = getFileList(inputDir);

// Process files that match the pattern
for (i = 0; i < list.length; i++) {
    if (endsWith(list[i], "dapi.tiff")) {
        // Get base filename
        baseName = substring(list[i], 0, lastIndexOf(list[i], "_dapi.tiff"));
        
        // Construct magenta filename
        magentaFile = inputDir + baseName + "_magenta.tiff";
        
        // Check if matching magenta file exists
        if (File.exists(magentaFile)) {
            // Process the pair
            open(inputDir + list[i]);
            run("Blue");
            open(magentaFile);
            run("Merge Channels...", "c3=" + list[i] + " c6=" + baseName + "_magenta.tiff keep");
            
            // Save merged image
            outputName = inputDir + baseName + "_dapi_magenta_merge.tiff";
            saveAs("Tiff", outputName);
            
            // Close windows
            close();
            close();
            selectImage(list[i]);
            close();
            
            print("Processed: " + baseName);
        }
    }
}

print("Processing complete!");




###################################################################
###################################################################
################################################################### for processing g3bp1 IF images 

################################### z-stack to tiff
// Define input and output directories
input_dir = "/Users/soundhar/Desktop/26_g3bp_if_u20s_z8_low/btach1_stressor_1hrs_exposure/";
output_dir = "/Users/soundhar/Desktop/26_g3bp_if_u20s_z8_low/btach1_stressor_1hrs_exposure/processed_tiff/";

// Create output directory if it doesn't exist
File.makeDirectory(output_dir);

// Get a list of all files in the input directory
list = getFileList(input_dir);

// Loop through each file
for (i = 0; i < list.length; i++) {
    if (endsWith(list[i], ".czi")) {
        // Open the .czi file
        path = input_dir + list[i];
        open(path);

        // Get the base name of the file (without extension)
        base_name = File.nameWithoutExtension;

        // Perform maximum intensity projection
        run("Z Project...", "projection=[Max Intensity]");

        // Apply the Magenta Hot LUT
        run("Magenta Hot");

        // Save the processed image
        saveAs("Tiff", output_dir + "MAX_" + base_name + ".tif");

        // Close all windows
        close("*");

        print("Processed: " + list[i]);
    }
}

print("All files processed.");

############################################ multipage to single page 
#!/bin/bash

# Input directory containing .tiff files
input_dir="./processed_tiff"

# Output directory for separated channels
output_dir="./processed_tiff"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through all .tiff files in the input directory
for file in "$input_dir"/*.tif; do
    # Get the base name of the file (without extension)
    base_name=$(basename "$file" .tif)

    # Split the channels using bfconvert
    /Users/soundhar/Desktop/softwares/bftools/bfconvert -channel 0 "$file" "$output_dir/${base_name}_dapi.tiff"
    /Users/soundhar/Desktop/softwares/bftools/bfconvert -channel 1 "$file" "$output_dir/${base_name}_mcherry.tiff"
    /Users/soundhar/Desktop/softwares/bftools/bfconvert -channel 2 "$file" "$output_dir/${base_name}_magenta.tiff"
    /Users/soundhar/Desktop/softwares/bftools/bfconvert -channel 3 "$file" "$output_dir/${base_name}_g3bp1.tiff"

    echo "Processed: $file"
done

echo "All files processed."



######################### use this on g3bp1 it will take DAPI and gebp1 tiff file do color chnage and merge them and save 
// Set input/output directories
inputDir = "/Users/soundhar/Desktop/26_g3bp_if_u20s_z8_low/btach1_stressor_1hrs_exposure/processed_tiff/";
outputDir = inputDir + "green_merged/";
File.makeDirectory(outputDir);

// Get list of files in the input directory
list = getFileList(inputDir);

// Loop through files
for (i = 0; i < list.length; i++) {
    if (endsWith(list[i], "_dapi.tiff")) {
        // Extract base name for file pairing
        base = substring(list[i], 0, lastIndexOf(list[i], "_dapi.tiff"));
        g3bpFile = inputDir + base + "_g3bp1.tiff";
        
        if (File.exists(g3bpFile)) {
            // Open DAPI image
            open(inputDir + list[i]);
            dapiImage = list[i]; // Store filename
            selectWindow(dapiImage);
            run("Blue"); // Set DAPI channel to Blue
            
            // Open G3BP1 image
            open(g3bpFile);
            g3bpImage = base + "_g3bp1.tiff";
            selectWindow(g3bpImage);
            run("Green"); // Set G3BP1 channel to Green
            
            // Merge channels: DAPI (blue) = c1, G3BP1 (green) = c2
            run("Merge Channels...", "c1=" + dapiImage + " c2=" + g3bpImage + " create");

            // Flatten composite image to RGB (single-page)
            run("Flatten");

            // Save as single-page TIFF
            saveAs("Tiff", outputDir + base + "_green_merged.tiff");

            // Close all images
            close("*");

            // Print confirmation
            print("Processed: " + base);
        }
    }
}
print("Done! Green-merged files in: " + outputDir);

############################################ g3bp1 to green color 
// Set input directory
inputDir = "/Users/soundhar/Desktop/26_g3bp_if_u20s_z8_low/btach1_stressor_1hrs_exposure/processed_tiff/";
list = getFileList(inputDir);

// Loop through all files
for (i = 0; i < list.length; i++) {
    if (endsWith(list[i], "_g3bp1.tiff")) {
        fullPath = inputDir + list[i];

        // Open the image
        open(fullPath);

        // Apply green LUT
        run("Green");

        // Save as TIFF (overwrite original file)
        saveAs("Tiff", fullPath);

        // Close the image
        close();
        
        print("Processed: " + list[i]);
    }
}
print("Done! All _g3bp1.tiff files colorized green and saved.");
