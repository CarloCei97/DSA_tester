import os

# Directory containing screenshots
base_dir = "screenshots"

# Supported image file extensions
image_extensions = (".png", ".jpg", ".jpeg", ".bmp", ".tif", ".tiff")

# Function to create LaTeX figure block
def create_figure_block(image_path, folder_name, label, caption):
    return f"""\\begin{{figure}}
    \\centering
    \\includegraphics[width=1\\linewidth]{{{image_path}}}
    \\caption{{Dataset {folder_name} - {caption}}}
    \\label{{fig:{label}}}
\\end{{figure}}\n"""

# Function to edit the file name for the caption
def edit_file_name(file_name):
    replacements = {
        "_": " ",
        "curr": "current",
        "volt": "voltage",
        "soc": "state of charge",
        "screenshot": ""
    }
    for old, new in replacements.items():
        file_name = file_name.replace(old, new)
    return file_name

# Function to edit the folder name for LaTeX compatibility
def edit_folder_name(folder_name):
    return folder_name.replace("_", "\\_")

# Walk through the directory and collect LaTeX blocks
latex_blocks = []
for root, _, files in os.walk(base_dir):
    for file in files:
        if file.endswith(image_extensions):
            relative_path = os.path.join(root, file).replace("\\", "/")
            folder_name = os.path.basename(root)
            
            # Edit folder name for LaTeX
            edited_folder_name = edit_folder_name(folder_name)
            
            # Create label
            image_name = os.path.splitext(file)[0]
            label = f"{folder_name}_{image_name}"
            
            # Create caption
            edited_file_name = edit_file_name(image_name)
            caption = f"{edited_file_name}"
            
            # Generate LaTeX block
            latex_blocks.append(create_figure_block(relative_path, edited_folder_name, label, caption))

# Save the LaTeX code to a file
output_file = "figures.tex"
with open(output_file, "w", encoding="utf-8") as f:
    f.writelines(latex_blocks)

print(f"LaTeX code for figures saved to {output_file}")
