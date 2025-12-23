#!/bin/bash

# Script to find missing files with a given naming pattern

# ==========================================
# CONFIGURATION (modify these values as needed)
# ==========================================

# 1. Directory to search for files
DIR="/eos/user/a/asmirnov/reco_data_lambda_c/analysed_data/signal/root"  # Replace with your directory
# DIR="/eos/user/a/asmirnov/reco_data_lambda_c/sig_data"

# 2. Index range to check
MIN_INDEX=0   # Minimum index
MAX_INDEX=480 # Maximum index

# 3. File name pattern (without the index)
#    Use %d as placeholder for the index
# FILE_PREFIX="reco_full_"
FILE_PREFIX="analysed_signal_"
FILE_SUFFIX=".root"
# Full pattern will be: analysed_signal_%d.root
# A regular expression is used to extract the index

# ==========================================
# MAIN CODE (no changes needed)
# ==========================================

echo "=== Searching for missing files ==="
echo "Directory: $DIR"
echo "File pattern: ${FILE_PREFIX}<index>${FILE_SUFFIX}"
echo "Index range: $MIN_INDEX - $MAX_INDEX"
echo "=========================================="

# Check if directory exists
if [ ! -d "$DIR" ]; then
    echo "Error: Directory '$DIR' does not exist!"
    exit 1
fi

# Verify indices are numbers
if ! [[ "$MIN_INDEX" =~ ^[0-9]+$ ]] || ! [[ "$MAX_INDEX" =~ ^[0-9]+$ ]]; then
    echo "Error: Indices must be integers!"
    exit 1
fi

# Verify minimum index is not greater than maximum index
if [ "$MIN_INDEX" -gt "$MAX_INDEX" ]; then
    echo "Error: Minimum index ($MIN_INDEX) is greater than maximum index ($MAX_INDEX)!"
    exit 1
fi

# Create regex pattern for matching files
# Escape special characters in prefix and suffix
PREFIX_ESC=$(echo "$FILE_PREFIX" | sed 's/[.[\*^$]/\\&/g')
SUFFIX_ESC=$(echo "$FILE_SUFFIX" | sed 's/[.[\*^$]/\\&/g')
REGEX_PATTERN="${PREFIX_ESC}([0-9]+)${SUFFIX_ESC}"

# Array to store found indices
found_indices=()

# Change to target directory
cd "$DIR" || exit 1

echo "Searching for files matching: ${FILE_PREFIX}*${FILE_SUFFIX}"

# Find all files matching the pattern and extract indices
for file in "${FILE_PREFIX}"*"${FILE_SUFFIX}"; do
    # Check if it's a file (not directory) and not an unmatched pattern
    if [ -f "$file" ] && [[ "$file" != "${FILE_PREFIX}*${FILE_SUFFIX}" ]]; then
        # Extract index from filename
        if [[ "$file" =~ $REGEX_PATTERN ]]; then
            index="${BASH_REMATCH[1]}"
            found_indices+=("$index")
        else
            echo "Warning: file '$file' does not match expected format"
        fi
    fi
done

# If no files were found
if [ ${#found_indices[@]} -eq 0 ]; then
    echo "No files found with pattern ${FILE_PREFIX}*${FILE_SUFFIX}"
    echo ""
    echo "All files in range $MIN_INDEX-$MAX_INDEX are missing!"

    # Show all expected indices as missing
    echo "Missing indices:"
    if [ $((MAX_INDEX - MIN_INDEX)) -lt 50 ]; then
        # For small ranges, show all indices
        for (( i=MIN_INDEX; i<=MAX_INDEX; i++ )); do
            echo -n "$i "
        done
        echo ""
    else
        # For large ranges, show only beginning and end
        echo "$MIN_INDEX ... $MAX_INDEX (all $((MAX_INDEX - MIN_INDEX + 1)) files)"
    fi
    exit 0
fi

# Sort indices numerically
IFS=$'\n' sorted_indices=($(sort -n <<< "${found_indices[*]}"))
unset IFS

echo "Files found: ${#sorted_indices[@]}"

# Check for missing files in the specified range
missing_indices=()
found_in_range=0

for (( i=MIN_INDEX; i<=MAX_INDEX; i++ )); do
    found=false
    for idx in "${sorted_indices[@]}"; do
        if [ "$idx" -eq "$i" ]; then
            found=true
            ((found_in_range++))
            break
        fi
    done

    if [ "$found" = false ]; then
        missing_indices+=("$i")
    fi
done

# Output results
echo "Files found in specified range: $found_in_range"
echo ""

if [ ${#missing_indices[@]} -eq 0 ]; then
    echo "✓ All files in range $MIN_INDEX-$MAX_INDEX are present!"
else
    echo "✗ Missing files (indices):"

    # Group consecutive missing indices
    prev_index=-2
    range_start=-1
    range_count=0

    for idx in "${missing_indices[@]}"; do
        if [ $((idx - prev_index)) -eq 1 ]; then
            # Continue the range
            prev_index=$idx
            ((range_count++))
        else
            # Finish previous range (if existed)
            if [ "$range_start" -ne -1 ]; then
                if [ "$range_start" -eq "$prev_index" ]; then
                    echo "  $range_start"
                else
                    echo "  $range_start-$prev_index ($((prev_index - range_start + 1)) files)"
                fi
            fi

            # Start new range
            range_start=$idx
            prev_index=$idx
            range_count=1
        fi
    done

    # Output last range
    if [ "$range_start" -ne -1 ]; then
        if [ "$range_start" -eq "$prev_index" ]; then
            echo "  $range_start"
        else
            echo "  $range_start-$prev_index ($((prev_index - range_start + 1)) files)"
        fi
    fi

    echo ""
    echo "Total missing: ${#missing_indices[@]} files"
fi

# Additional statistics
echo ""
echo "=========================================="
echo "Statistics:"
echo "Checked range: $MIN_INDEX - $MAX_INDEX"
echo "Expected number of files: $((MAX_INDEX - MIN_INDEX + 1))"
echo "Files found in range: $found_in_range"
echo "Missing files: ${#missing_indices[@]}"
echo "Completion percentage: $((found_in_range * 100 / (MAX_INDEX - MIN_INDEX + 1)))%"

# Show minimum and maximum found indices
if [ ${#sorted_indices[@]} -gt 0 ]; then
    echo ""
    echo "Additional information:"
    echo "Minimum found index: ${sorted_indices[0]}"
    echo "Maximum found index: ${sorted_indices[-1]}"

    # Check for files outside specified range
    outside_min=()
    outside_max=()

    for idx in "${sorted_indices[@]}"; do
        if [ "$idx" -lt "$MIN_INDEX" ]; then
            outside_min+=("$idx")
        elif [ "$idx" -gt "$MAX_INDEX" ]; then
            outside_max+=("$idx")
        fi
    done

    if [ ${#outside_min[@]} -gt 0 ]; then
        echo "Files below minimum index: ${outside_min[*]}"
    fi

    if [ ${#outside_max[@]} -gt 0 ]; then
        echo "Files above maximum index: ${outside_max[*]}"
    fi
fi
