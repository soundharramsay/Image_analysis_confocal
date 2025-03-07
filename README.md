# bioconvert_czi_2_seperate_tiff
bioconvert_czi_2_seperate_tiff

###################################### script is not doing recoloring properly ---but tiff conversions are good 

#!/bin/bash

# Input and output directories
input_dir="/Users/soundhar/Desktop/18_u20s_stable_halo_z8_stressor_feb18/sorbitol_working"
output_dir="/Users/soundhar/Desktop/18_u20s_stable_halo_z8_stressor_feb18/sorbitol_working/output_tiff"

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

    # Recolor the 647 channel to hot magenta using ImageMagick (magick command)
    if [[ -f "$output_dir/${base_name}_647.tiff" ]]; then
        echo "Recoloring 647 channel for: ${base_name}_647.tiff"
        magick "$output_dir/${base_name}_647.tiff" -colorspace RGB -channel R -evaluate Multiply 1.0 -channel G -evaluate Multiply 0.0 -channel B -evaluate Multiply 1.0 -colorspace sRGB "$output_dir/${base_name}_647_magenta.tiff"

        # Verify if the recolored file was created
        if [[ -f "$output_dir/${base_name}_647_magenta.tiff" ]]; then
            echo "Recolored 647 channel to hot magenta: ${base_name}_647_magenta.tiff"
        else
            echo "Error: Recoloring failed for ${base_name}_647.tiff"
        fi
    else
        echo "Warning: 647 channel not found for ${base_name}"
    fi

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

            // Save the merged image
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
