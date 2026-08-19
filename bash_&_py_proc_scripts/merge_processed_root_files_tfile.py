import ROOT
import glob
import os

# Define paths
input_pattern = "/home/ome123/Рабочий стол/SPD/data/analysed_background_*.root"
output_file = "/home/ome123/Рабочий стол/SPD/data/analysed_background.root"

# Get input files EXCLUDING the output file
all_files = glob.glob(input_pattern)

# Restrict files amount (currently using all files, can be modified if needed)
input_files = all_files[:]

print(f"Found {len(input_files)} input files")

# Remove existing output file if it exists
if os.path.exists(output_file):
    print(f"Removing existing output file: {output_file}")
    os.remove(output_file)

# Use ROOT's TFileMerger
print(f"\nMerging using ROOT's TFileMerger...")
merger = ROOT.TFileMerger(False)  # False = not verbose
merger.OutputFile(output_file)

for input_file in input_files:
    merger.AddFile(input_file)

# Merge the files
success = merger.Merge()

if success:
    print(f"Successfully merged files to: {output_file}")
else:
    print(f"Failed to merge files")
    exit(1)

# ===================== VERIFICATION FUNCTIONS =====================

def check_merge_with_root(merged_file, input_files):
    """Check merge using ROOT's built-in classes"""

    print("\n" + "=" * 60)
    print("Checking merge with ROOT...")
    print("=" * 60)

    # Open merged file
    f_merged = ROOT.TFile.Open(merged_file)
    if not f_merged or f_merged.IsZombie():
        print("❌ Error: Cannot open merged file")
        return False

    print("✓ Merged file opened successfully")

    # Get list of trees in merged file
    merged_trees = []
    merged_objects = {}
    for key in f_merged.GetListOfKeys():
        obj = key.ReadObj()
        obj_name = key.GetName()
        obj_class = obj.ClassName()
        merged_objects[obj_name] = obj_class

        if "TTree" in obj_class:
            merged_trees.append(obj_name)

    print(f"\nObjects in merged file:")
    for obj_name, obj_class in merged_objects.items():
        print(f"  {obj_name} ({obj_class})")

    # Count events in merged file
    total_merged_events = 0
    for tree_name in merged_trees:
        tree = f_merged.Get(tree_name)
        events = tree.GetEntries()
        total_merged_events += events
        print(f"\n  Tree '{tree_name}':")
        print(f"    Events: {events}")
        print(f"    Branches: {tree.GetListOfBranches().GetEntries()}")

    # Count events in input files
    total_input_events = 0
    input_trees_info = {}

    print(f"\n{'='*60}")
    print("Checking input files...")
    print(f"{'='*60}")

    for i, input_file in enumerate(input_files):
        f_input = ROOT.TFile.Open(input_file)
        if not f_input or f_input.IsZombie():
            print(f"❌ Warning: Cannot open input file {input_file}")
            continue

        file_events = 0
        for key in f_input.GetListOfKeys():
            obj = key.ReadObj()
            if "TTree" in obj.ClassName():
                events = obj.GetEntries()
                file_events += events
                tree_name = key.GetName()
                if tree_name not in input_trees_info:
                    input_trees_info[tree_name] = 0
                input_trees_info[tree_name] += events

        total_input_events += file_events
        print(f"  File {i+1:3d}/{len(input_files)}: {os.path.basename(input_file):30s} - {file_events:8d} events")
        f_input.Close()

    print(f"\n{'='*60}")
    print("SUMMARY:")
    print(f"{'='*60}")
    print(f"Total input events: {total_input_events}")
    print(f"Total merged events: {total_merged_events}")

    # Check for each tree individually
    all_trees_match = True
    for tree_name in set(list(input_trees_info.keys()) + merged_trees):
        input_count = input_trees_info.get(tree_name, 0)
        if tree_name in merged_trees:
            tree = f_merged.Get(tree_name)
            merged_count = tree.GetEntries() if tree else 0
        else:
            merged_count = 0

        status = "✓" if input_count == merged_count else "❌"
        match = "MATCH" if input_count == merged_count else f"MISMATCH (input: {input_count}, merged: {merged_count})"
        print(f"{status} Tree '{tree_name}': {match}")

        if input_count != merged_count:
            all_trees_match = False

    print(f"\nOverall events match: {total_input_events == total_merged_events}")

    # Check structure consistency
    print(f"\n{'='*60}")
    print("STRUCTURE CHECK:")
    print(f"{'='*60}")

    # Check first input file structure vs merged file
    if input_files:
        first_input = ROOT.TFile.Open(input_files[0])
        if first_input:
            print("Comparing structure with first input file...")

            # Compare number of objects
            input_keys = [key.GetName() for key in first_input.GetListOfKeys()]
            merged_keys = [key.GetName() for key in f_merged.GetListOfKeys()]

            if set(input_keys) == set(merged_keys):
                print("✓ Object names match")
            else:
                print("⚠️  Object names differ:")
                print(f"   Input objects: {sorted(input_keys)}")
                print(f"   Merged objects: {sorted(merged_keys)}")

            # Compare tree branches for each tree
            for tree_name in merged_trees:
                if tree_name in input_keys:
                    input_tree = first_input.Get(tree_name)
                    merged_tree = f_merged.Get(tree_name)

                    if input_tree and merged_tree:
                        input_branches = [br.GetName() for br in input_tree.GetListOfBranches()]
                        merged_branches = [br.GetName() for br in merged_tree.GetListOfBranches()]

                        if set(input_branches) == set(merged_branches):
                            print(f"✓ Tree '{tree_name}' branches match ({len(input_branches)} branches)")
                        else:
                            print(f"⚠️  Tree '{tree_name}' branches differ:")
                            missing_in_merged = set(input_branches) - set(merged_branches)
                            extra_in_merged = set(merged_branches) - set(input_branches)
                            if missing_in_merged:
                                print(f"   Missing in merged: {missing_in_merged}")
                            if extra_in_merged:
                                print(f"   Extra in merged: {extra_in_merged}")

            first_input.Close()

    f_merged.Close()

    return all_trees_match and (total_input_events == total_merged_events)

def basic_file_checks(merged_file, input_files):
    """Perform basic file existence and size checks"""

    print("\n" + "=" * 60)
    print("Performing basic file checks...")
    print("=" * 60)

    # Check if merged file exists
    if not os.path.exists(merged_file):
        print("❌ Merged file does not exist")
        return False
    print("✓ Merged file exists")

    # Check file size
    merged_size = os.path.getsize(merged_file)
    print(f"\nMerged file:")
    print(f"  Path: {merged_file}")
    print(f"  Size: {merged_size:,} bytes ({merged_size / (1024**3):.3f} GB)")

    # Check input files
    total_input_size = 0
    existing_files = 0
    missing_files = []

    print(f"\nInput files:")
    for i, input_file in enumerate(input_files):
        if os.path.exists(input_file):
            size = os.path.getsize(input_file)
            total_input_size += size
            existing_files += 1
            print(f"  {i+1:3d}/{len(input_files)}: {os.path.basename(input_file):30s} - {size/1024**2:8.2f} MB")
        else:
            missing_files.append(input_file)
            print(f"  {i+1:3d}/{len(input_files)}: {os.path.basename(input_file):30s} - ❌ NOT FOUND")

    print(f"\nSummary:")
    print(f"  Existing input files: {existing_files}/{len(input_files)}")
    if missing_files:
        print(f"  Missing files: {len(missing_files)}")
        for i, missing in enumerate(missing_files[:3]):  # Show first 3 missing
            print(f"    {i+1}. {missing}")
        if len(missing_files) > 3:
            print(f"    ... and {len(missing_files)-3} more")

    print(f"  Total input size: {total_input_size:,} bytes ({total_input_size / (1024**3):.3f} GB)")
    print(f"  Merged file size: {merged_size:,} bytes ({merged_size / (1024**3):.3f} GB)")

    # Size ratio check
    if total_input_size > 0:
        ratio = merged_size / total_input_size
        print(f"\nSize analysis:")
        print(f"  Size ratio (merged/input): {ratio:.3f}")

        # ROOT files typically have some overhead
        if 0.85 < ratio < 1.15:  # Slightly wider range for safety
            print(f"  ✓ File size ratio is reasonable")
            size_ok = True
        else:
            print(f"  ⚠️  File size ratio seems unusual")
            print(f"     Expected: ~0.95-1.05 (due to ROOT compression/overhead)")
            size_ok = False
    else:
        print("❌ Total input size is 0")
        size_ok = False

    return existing_files == len(input_files) and size_ok

def check_file_integrity(file_path):
    """Check if a ROOT file is not corrupted"""
    try:
        f = ROOT.TFile.Open(file_path)
        if not f or f.IsZombie():
            return False
        # Try to read something
        keys = f.GetListOfKeys()
        if keys:
            # Try to read the first object
            key = keys.At(0)
            if key:
                obj = key.ReadObj()
                if obj:
                    return True
        f.Close()
        return True
    except:
        return False

def complete_verification(merged_file, input_files):
    """Complete verification of merged file"""

    print("\n" + "=" * 80)
    print("COMPLETE MERGE VERIFICATION")
    print("=" * 80)

    overall_success = True
    results = {}

    # 1. Basic file checks
    print("\n" + "=" * 60)
    print("STEP 1: BASIC FILE CHECKS")
    print("=" * 60)
    basic_ok = basic_file_checks(merged_file, input_files)
    results["Basic checks"] = basic_ok
    if not basic_ok:
        print("❌ Basic file checks failed")
        overall_success = False
    else:
        print("✓ Basic file checks passed")

    # 2. File integrity check
    print("\n" + "=" * 60)
    print("STEP 2: FILE INTEGRITY CHECK")
    print("=" * 60)
    print("Checking merged file integrity...")
    merged_integrity = check_file_integrity(merged_file)
    if merged_integrity:
        print("✓ Merged file integrity check passed")
    else:
        print("❌ Merged file appears to be corrupted")
        overall_success = False
    results["Merged file integrity"] = merged_integrity

    # 3. ROOT-based verification
    print("\n" + "=" * 60)
    print("STEP 3: ROOT-BASED VERIFICATION")
    print("=" * 60)
    try:
        root_ok = check_merge_with_root(merged_file, input_files)
        results["ROOT verification"] = root_ok
        if not root_ok:
            overall_success = False
    except Exception as e:
        print(f"❌ ROOT verification failed with error: {e}")
        import traceback
        traceback.print_exc()
        results["ROOT verification"] = False
        overall_success = False

    # 4. Final summary
    print("\n" + "=" * 80)
    print("FINAL VERIFICATION SUMMARY")
    print("=" * 80)

    print("\nCheck Results:")
    for check_name, check_result in results.items():
        status = "✓ PASS" if check_result else "❌ FAIL"
        print(f"  {status} - {check_name}")

    print("\n" + "=" * 80)
    if overall_success:
        print("🎉 MERGE VERIFICATION COMPLETED SUCCESSFULLY!")
        print(f"   Merged file: {merged_file}")
        print(f"   Total input files: {len(input_files)}")
    else:
        print("❌ MERGE VERIFICATION FAILED!")
        print("   Please check the errors above.")

    print("=" * 80 + "\n")

    return overall_success

# ===================== RUN VERIFICATION =====================

print("\n" + "=" * 80)
print("STARTING VERIFICATION PROCESS")
print("=" * 80)

# Run complete verification
success = complete_verification(output_file, input_files)

if success:
    print("✅ Merge and verification completed successfully!")
else:
    print("❌ Merge verification failed. Check the errors above.")
    exit(1)
