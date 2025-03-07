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
