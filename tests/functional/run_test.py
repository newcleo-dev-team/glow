import subprocess
from pathlib import Path
import hashlib
import os

# Function to compute the SHA256 hash of a file
def compute_hash(file_path):
    sha256 = hashlib.sha256()
    with open(file_path, "rb") as f:
        for block in iter(lambda: f.read(4096), b""):
            sha256.update(block)
    return sha256.hexdigest()

# Folder containing the test scripts
script_folder = Path(os.getcwd() + "/tests/functional/")
script_files = sorted(script_folder.glob("test_*.py"))

# Counters
total = 0
passed = 0
failed = 0
skipped = 0

print("\nStarting tests execution...\n")

for script in script_files:
    total += 1
    base_name = script.stem

    output_file = script_folder / f"{base_name}.dat"
    reference_file = script_folder / f"{base_name}_ref.dat"

    print(f"[{total}] Running test: {script.name}")

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
        print(f"Geometry generated successfully: {output_file.name}")

        # Save log
        log_file = script.with_suffix(".log")
        log_file.write_text(result.stdout + "\n--- STDERR ---\n" + result.stderr)

    except subprocess.CalledProcessError as e:
        print(f"Script failed: {script.name}")
        print(f"stderr: {e.stderr.strip()}")
        failed += 1
        continue

    if not output_file.exists():
        print(f"Output file missing: {output_file.name}")
        failed += 1
        continue

    if not reference_file.exists():
        print(f"Reference file missing: {reference_file.name} — skipping hash comparison.")
        skipped += 1
        continue

    # Compare hashes
    hash_generated = compute_hash(output_file)
    hash_reference = compute_hash(reference_file)

    if hash_generated == hash_reference:
        print(f"Test passed.\n")
        passed += 1
    else:
        print(f"Test failed.")
        print(f"  Generated geometry hash:   {hash_generated}")
        print(f"  Reference geometry hash:   {hash_reference}")
        failed += 1

print("\nTests Summary")
print("────────────────────────────")
print(f"Total tests run:     {total}")
print(f"Passed:              {passed}")
print(f"Failed:              {failed}")
print(f"Skipped (no ref):    {skipped}")
print("────────────────────────────\n")
    

