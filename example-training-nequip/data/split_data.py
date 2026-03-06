import os
import random

def read_xyz_blocks(filename):
    blocks = []
    with open(filename, 'r') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        n_atoms = int(lines[i].strip())
        block = lines[i:i + n_atoms + 2]
        blocks.append(block)
        i += n_atoms + 2
    return blocks

def write_xyz_blocks(filename, blocks):
    with open(filename, 'w') as f:
        for block in blocks:
            f.writelines(block)

def split_blocks(blocks, train_ratio=0.70, val_ratio=0.15, test_ratio=0.15):
    assert abs(train_ratio + val_ratio + test_ratio - 1.0) < 1e-6
    random.shuffle(blocks)
    total = len(blocks)
    train_end = int(train_ratio * total)
    val_end = train_end + int(val_ratio * total)
    return blocks[:train_end], blocks[train_end:val_end], blocks[val_end:]

# === Main ===
random.seed(13)  # For reproducibility

input_xyz = "nequip_dbh_404.xyz"  # Replace with your input file
output_dir = "split_two_state_data"

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Read and split
blocks = read_xyz_blocks(input_xyz)
train, val, test = split_blocks(blocks)

# Write to output directory
write_xyz_blocks(os.path.join(output_dir, "train.xyz"), train)
write_xyz_blocks(os.path.join(output_dir, "val.xyz"), val)
write_xyz_blocks(os.path.join(output_dir, "test.xyz"), test)

print("✅ Splitting complete. Files saved to:", output_dir)

