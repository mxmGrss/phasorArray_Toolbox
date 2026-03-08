import os
import subprocess
import pathlib
from datetime import datetime

def get_loc(file_path):
    try:
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            return sum(1 for _ in f)
    except:
        return 0

def get_size_str(file_path):
    size_bytes = os.path.getsize(file_path)
    if size_bytes >= 1024 * 1024:
        return f"{size_bytes / (1024 * 1024):.1f}MB"
    elif size_bytes >= 1024:
        return f"{size_bytes / 1024:.1f}KB"
    else:
        return f"{size_bytes}B"

def get_dependencies(file_name, search_dir):
    """Finds files that call this function (filename without extension)."""
    func_name = pathlib.Path(file_name).stem
    ext = pathlib.Path(file_name).suffix
    if ext not in ['.m', '.py']:
        return "—"
    
    try:
        # Using Get-ChildItem -Recurse | Select-String for compatibility
        ps_cmd = f'Get-ChildItem -Path "{search_dir}\\*" -Include *.m,*.py -Recurse | Select-String -Pattern "\\b{func_name}\\b" -List | % {{ $_.Path }}'
        cmd = ["powershell", "-Command", ps_cmd]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            lines = result.stdout.strip().split('\n')
            files = []
            for line in lines:
                p = line.strip()
                if not p: continue
                b = os.path.basename(p)
                if b != file_name:
                    files.append(b)
            
            if not files:
                return "None"
            # Remove duplicates and format
            files = list(dict.fromkeys(files)) 
            return ", ".join(files[:5]) + ("..." if len(files) > 5 else "")
        return "None"
    except Exception as e:
        return f"Error: {str(e)}"

def generate_audit_report(target_folder, output_file, search_root=None):
    target_path = pathlib.Path(target_folder).resolve()
    if not search_root:
        search_root = target_path
    search_root = pathlib.Path(search_root).resolve()
    
    report_lines = [f"# Systematic Audit: {target_path.name}\n", 
                   f"Target Folder: `{target_path}`\n",
                   f"Search Root (Dependencies): `{search_root}`\n",
                   f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"]
    
    # Track files grouped by subfolder
    grouped_files = {}

    for root, dirs, files in os.walk(target_path):
        # Skip hidden directories
        dirs[:] = [d for d in dirs if not d.startswith('.')]
        
        rel_root = os.path.relpath(root, target_path)
        if rel_root == ".":
            rel_root = "/"
        
        grouped_files[rel_root] = []
        for file in files:
            file_path = os.path.join(root, file)
            size = get_size_str(file_path)
            loc = get_loc(file_path) if not file.endswith(('.mat', '.pdf', '.png', '.jpg', '.eps')) else "—"
            # Search dependencies in search_root
            deps = get_dependencies(file, search_root)
            
            # Calculate relative path from the output report's directory to the file
            # Assuming output_file is in current directory or specified path
            out_dir = os.path.dirname(os.path.abspath(output_file))
            rel_link_path = os.path.relpath(file_path, out_dir).replace('\\', '/')
            
            grouped_files[rel_root].append({
                "File": file,
                "Link": rel_link_path,
                "Size": size,
                "LOC": loc,
                "Deps": deps
            })

    for folder, files_list in grouped_files.items():
        if not files_list:
            continue
            
        report_lines.append(f"## {folder} ({len(files_list)} files)\n")
        report_lines.append("| File | Size | LOC | Status | Dependencies | Notes |")
        report_lines.append("|------|------|-----|--------|--------------|-------|")
        
        for f in files_list:
            report_lines.append(f"| [`{f['File']}`]({f['Link']}) | {f['Size']} | {f['LOC']} | `TO READ` | {f['Deps']} | — |")
        report_lines.append("\n")

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("\n".join(report_lines))
    print(f"Report generated: {output_file}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python audit_toolbox.py <folder_path> [search_root]")
    else:
        folder = sys.argv[1]
        s_root = sys.argv[2] if len(sys.argv) > 2 else None
        out = "file_triage.md"
        generate_audit_report(folder, out, s_root)
