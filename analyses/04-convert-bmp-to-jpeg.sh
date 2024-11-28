#!/bin/bash

# Start the conversion process in the background
echo "Starting conversion of BMP files to JPEG..."
mogrify -format jpeg ../article-FSIGEN/*.bmp &

# Prompt the user for deleting the original BMP files (while conversion is ongoing)
read -t 10 -p "Do you want to delete the original BMP files after conversion? (Y/N): " delete_bmp

# Convert the response to uppercase and default to 'N' if no response
delete_bmp=$(echo "${delete_bmp:-N}" | tr '[:lower:]' '[:upper:]')

# Prompt the user for removing JPEG metadata
read -t 10 -p "Do you want to remove metadata from the JPEG files? (Y/N): " remove_metadata

# Convert the response to uppercase and default to 'N' if no response
remove_metadata=$(echo "${remove_metadata:-N}" | tr '[:lower:]' '[:upper:]')

# Wait for the conversion process to complete
wait
echo "Conversion completed."

# Handle BMP file deletion
if [[ "$delete_bmp" == "Y" ]]; then
    rm ../article-FSIGEN/*.bmp
    echo "Original BMP files deleted."
else
    echo "Original BMP files retained."
fi

# Handle metadata removal from JPEG files
if [[ "$remove_metadata" == "Y" ]]; then
    echo "Removing metadata from JPEG files..."
    for file in ../article-FSIGEN/*.jpeg; do
        mat2 "$file"  # Create cleaned file with ".cleaned" inserted before extension
        rm "$file"  # Remove the original non-cleaned file
        mv "${file%.jpeg}.cleaned.jpeg" "$file"  # Rename cleaned file to original name
    done
    echo "Metadata removed from all JPEG files."
else
    echo "JPEG metadata retained."
fi

