#!/bin/bash

TARGET_DIR="/eos/user/a/asmirnov/reco_data_lambda_c/analysed_data/background/root"
FILE_PATTERN="analysed_background_*.root"
# TARGET_DIR="/eos/user/a/asmirnov/reco_data_lambda_c/sig_data"
# FILE_PATTERN="reco_full_*.root"

echo "Checking for empty $FILE_PATTERN files in: $TARGET_DIR"
echo "------------------------------------------------------"

# Counters
total_count=0
empty_count=0
non_empty_count=0
total_events=0
files_with_events=0

# Arrays to store files
declare -a empty_files=()
declare -a files_event_counts=()

# Enable nullglob to handle case with no matching files
shopt -s nullglob

# Find and check files
for file in "$TARGET_DIR"/$FILE_PATTERN; do
    total_count=$((total_count + 1))

    # Check if file is empty (size is 0)
    if [ ! -s "$file" ]; then
        empty_count=$((empty_count + 1))
        empty_files+=("$file")
        echo "[EMPTY] $(basename "$file") (0 bytes)"
    else
        non_empty_count=$((non_empty_count + 1))
        
        # Try to get event count from ROOT file
        event_count=0
        # Method 1: Try using root-ls to count entries in a tree (e.g., "events" or "tree")
        # Adjust the tree name based on your actual ROOT file structure
        TREE_NAME="events"  # Change this to your actual tree name
        
        # Try different methods to extract event count
        if command -v root-ls &>/dev/null; then
            # Using root-ls (usually available in ROOT installations)
            event_count=$(root-ls -l "$file":"$TREE_NAME" 2>/dev/null | grep -E "^Tree.*entries.*$" | head -1 | sed -E 's/.*entries[[:space:]]+([0-9]+).*/\1/')
        fi
        
        # If root-ls failed or not available, try rootdump or python
        if [[ -z "$event_count" || "$event_count" == "0" || ! "$event_count" =~ ^[0-9]+$ ]]; then
            # Try using Python with uproot (if available)
            if command -v python3 &>/dev/null; then
                python_script=$(cat << 'EOF'
import uproot
import sys
try:
    file = uproot.open(sys.argv[1])
    # Try common tree names
    tree_names = ["events", "tree", "T", "Data", "AnalysisTree"]
    for tree_name in tree_names:
        if tree_name in file:
            print(file[tree_name].num_entries)
            sys.exit(0)
    # If no standard tree name found, try the first tree
    for key in file.keys():
        if ";" in key:  # ROOT object with class type
            tree_name = key.split(";")[0]
            if tree_name in file:
                print(file[tree_name].num_entries)
                sys.exit(0)
    print("0")
except Exception as e:
    print("0")
EOF
                )
                event_count=$(python3 -c "$python_script" "$file" 2>/dev/null)
            fi
        fi
        
        # If still no event count, try using root -b -q (ROOT macro)
        if [[ -z "$event_count" || ! "$event_count" =~ ^[0-9]+$ ]]; then
            if command -v root &>/dev/null; then
                root_macro=$(cat << 'EOF'
{
    TFile* file = TFile::Open("%s");
    if (!file || file->IsZombie()) {
        cout << 0 << endl;
        return;
    }
    
    // Try to find a TTree
    TIter next(file->GetListOfKeys());
    TKey* key;
    Long64_t max_entries = 0;
    while ((key = (TKey*)next())) {
        TClass* cl = gROOT->GetClass(key->GetClassName());
        if (cl->InheritsFrom("TTree")) {
            TTree* tree = (TTree*)file->Get(key->GetName());
            if (tree->GetEntries() > max_entries) {
                max_entries = tree->GetEntries();
            }
        }
    }
    cout << max_entries << endl;
    file->Close();
}
EOF
                )
                # Create temporary macro
                tmp_macro=$(mktemp)
                printf "$root_macro" "$file" > "$tmp_macro"
                event_count=$(root -b -q -l "$tmp_macro" 2>/dev/null | tail -1)
                rm -f "$tmp_macro"
            fi
        fi
        
        # Default to 0 if still not found
        if [[ -z "$event_count" || ! "$event_count" =~ ^[0-9]+$ ]]; then
            event_count=0
        fi
        
        # Update counters
        if [ "$event_count" -gt 0 ]; then
            files_with_events=$((files_with_events + 1))
            total_events=$((total_events + event_count))
            files_event_counts+=("$(basename "$file") : $event_count events")
        fi
        
        # Get file size
        file_size=$(stat -c%s "$file" 2>/dev/null || stat -f%z "$file" 2>/dev/null)
        echo "[OK]     $(basename "$file") (${file_size} bytes, ${event_count} events)"
    fi
done

# Disable nullglob
shopt -u nullglob

echo "------------------------------------------------------"
echo "Summary:"
echo "  Total files found: $total_count"
echo "  Empty files: $empty_count"
echo "  Non-empty files: $non_empty_count"
echo "  Files with events: $files_with_events"
echo "  Total events in all files: $total_events"
if [ $files_with_events -gt 0 ]; then
    avg_events=$(echo "scale=2; $total_events / $files_with_events" | bc)
    echo "  Average events per file: $avg_events"
fi

# Print empty files
if [ $empty_count -gt 0 ]; then
    echo ""
    echo "List of all empty files:"
    echo "------------------------"
    for empty_file in "${empty_files[@]}"; do
        echo "$(basename "$empty_file")"
    done
    echo "------------------------"
    echo "Total empty files: $empty_count"
else
    echo ""
    echo "No empty files found."
fi

# Print event counts for files with events
if [ ${#files_event_counts[@]} -gt 0 ]; then
    echo ""
    echo "Files with event counts:"
    echo "------------------------"
    for event_info in "${files_event_counts[@]}"; do
        echo "$event_info"
    done
    echo "------------------------"
    echo "Total events in all files: $total_events"
fi

# If no files were found
if [ $total_count -eq 0 ]; then
    echo "Warning: No files matching pattern '$FILE_PATTERN' found in $TARGET_DIR"
fi
