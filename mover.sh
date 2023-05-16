#!/bin/bash

source_folder="$1"
destination_folder="$2"

# Check if source folder exists
if [ ! -d "$source_folder" ]; then
  echo "Source folder does not exist."
  exit 1
fi

# Check if destination folder exists
if [ ! -d "$destination_folder" ]; then
  echo "Destination folder does not exist."
  exit 1
fi

# Move files from source folder to destination folder
for file in "$source_folder"/*.txt; do
    mv "$file" "$destination_folder"
    echo "Moved file: $file"
done

for file in "$source_folder"/*.pep.xml; do
    mv "$file" "$destination_folder"
    echo "Moved file: $file"
done
