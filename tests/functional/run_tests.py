import hashlib
import os
import shutil
import subprocess

from pathlib import Path

def compute_hash(file_path: Path) -> str:
    """
    Function to compute the SHA256 hash of a file.

    Parameters
    ----------
    file_path : Path
        The ``Path`` object of the file to process.

    Returns
    -------
    str
        The SHA256 hash of a given file.
    """
    sha256 = hashlib.sha256()
    with open(file_path, "rb") as f:
        for block in iter(lambda: f.read(4096), b""):
            sha256.update(block)
    return sha256.hexdigest()

# Folder containing the test scripts
script_folder = Path(os.getcwd() + "/tests/functional/")
script_files = sorted(script_folder.glob("test_*.py"))
ref_folder = script_folder / "reference_dat"
output_folder = script_folder / "generated_dat"
# Create the folder where the generated .dat files will be saved, if it does
# not exist
if not output_folder.exists():
    os.mkdir(output_folder)
else:
    # Check if the path is a folder; raise an exception if it is not
    if not output_folder.is_dir():
        raise RuntimeError(
            f"The indicated '{output_folder}' path is not a folder.")

# Counters
total = 0
passed = 0
failed = 0
skipped = 0

print("\nStarting tests execution...\n")
# Loop over all the found script files
for script in script_files:
    total += 1
    base_name = script.stem
    # Build the paths for the reference .dat and the one to be generated
    output_file = script_folder / f"{base_name}.dat"
    reference_file = ref_folder / f"{base_name}_ref.dat"

    print(f"[{total}] Running test: {script.name}")
    # Skip the test if no reference .dat file is present
    if not reference_file.exists():
        print(f"Reference file missing: {reference_file.name} — "
              + "skipping test.")
        skipped += 1
        continue

    # Build the command to run in the SALOME shell environment
    command = f"salome shell {script}"
    try:
        result = subprocess.run(
            command,
            shell=True,
            check=True,
            executable="/bin/bash",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
    except subprocess.CalledProcessError as e:
        print(f"Script failed: {script.name}")
        print(f"stderr: {e.stderr.strip()}")
        failed += 1
        # Save log for the failed script
        log_file = script.with_suffix(".log")
        log_file.write_text(
            result.stdout + "\n--- STDERR ---\n" + result.stderr)
        continue

    # Verify the output .dat file has been generated in the target folder
    if not output_file.exists():
        print(f"Output file missing: {output_file.name}")
        failed += 1
        continue
    # Move the generated .dat in the output folder
    try:
        shutil.move(output_file, output_folder / output_file.name)
        print(f"Geometry generated successfully: {output_file.name}")
    except:
        raise RuntimeError(
            f"Error in moving generated '{output_file.name}' file.")

    # Compute and compare hashes of reference and generated .dat files
    hash_generated = compute_hash(output_folder / output_file.name)
    hash_reference = compute_hash(reference_file)
    # The test is passed if the hashes are the same, meaning the two files
    # are identical
    if hash_generated == hash_reference:
        print(f"Test passed.\n")
        passed += 1
    else:
        print(f"Test failed.")
        print(f"  Generated geometry hash:   {hash_generated}")
        print(f"  Reference geometry hash:   {hash_reference}\n")
        failed += 1

# Log the summary of the functional tests
print("\nTests Summary")
print("────────────────────────────")
print(f"Total tests run:     {total}")
print(f"Passed:              {passed}")
print(f"Failed:              {failed}")
print(f"Skipped (no ref):    {skipped}")
print("────────────────────────────\n")
