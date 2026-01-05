#!/usr/bin/env python3
"""
Script to find remaining hard-coded paths in the codebase.
Searches for the pattern '/cellar/users/zkoch/dream' in all Python files and notebooks.

Usage:
    python scripts/find_hardcoded_paths.py [--fix]

Options:
    --fix    Attempt to automatically fix the paths found (use with caution)
"""

import os
import re
import json
import glob
import argparse
from pathlib import Path


def find_hardcoded_paths_in_python(file_path: str) -> list:
    """Find hard-coded paths in a Python file."""
    findings = []
    with open(file_path, 'r') as f:
        for i, line in enumerate(f, 1):
            if '/cellar/users/zkoch/dream' in line:
                findings.append({
                    'file': file_path,
                    'line': i,
                    'content': line.strip(),
                    'type': 'code'
                })
    return findings


def find_hardcoded_paths_in_notebook(file_path: str) -> list:
    """Find hard-coded paths in a Jupyter notebook."""
    findings = []
    with open(file_path, 'r') as f:
        nb = json.load(f)

    for cell_idx, cell in enumerate(nb.get('cells', [])):
        cell_type = cell.get('cell_type', '')
        source = cell.get('source', [])
        if isinstance(source, list):
            source = ''.join(source)

        for line_idx, line in enumerate(source.split('\n'), 1):
            if '/cellar/users/zkoch/dream' in line:
                findings.append({
                    'file': file_path,
                    'cell': cell_idx,
                    'line': line_idx,
                    'content': line.strip()[:100] + ('...' if len(line.strip()) > 100 else ''),
                    'type': 'code' if cell_type == 'code' else 'markdown'
                })

        # also check outputs
        for output in cell.get('outputs', []):
            text = output.get('text', [])
            if isinstance(text, list):
                text = ''.join(text)
            if '/cellar/users/zkoch/dream' in text:
                findings.append({
                    'file': file_path,
                    'cell': cell_idx,
                    'content': '[output]',
                    'type': 'output'
                })

    return findings


def fix_path_in_python_file(file_path: str) -> int:
    """Fix hard-coded paths in a Python file. Returns number of fixes."""
    with open(file_path, 'r') as f:
        content = f.read()

    original = content

    # check if REPO_ROOT is defined
    has_repo_root = 'REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))' in content

    if not has_repo_root and '/cellar/users/zkoch/dream' in content:
        # add REPO_ROOT after imports
        import_match = re.search(r'^(import\s+os\s*$)', content, re.MULTILINE)
        if import_match:
            content = content[:import_match.end()] + '\n\n# repo root for relative paths\nREPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))' + content[import_match.end():]
        else:
            # add import os and REPO_ROOT at the top
            content = 'import os\n\n# repo root for relative paths\nREPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))\n\n' + content

    # replace paths
    content = content.replace("'/cellar/users/zkoch/dream/", "os.path.join(REPO_ROOT, '")
    content = content.replace('"/cellar/users/zkoch/dream/', 'os.path.join(REPO_ROOT, "')

    # fix closing parens
    lines = content.split('\n')
    fixed_lines = []
    for line in lines:
        if 'os.path.join(REPO_ROOT,' in line:
            open_parens = line.count('(')
            close_parens = line.count(')')
            diff = open_parens - close_parens
            if diff > 0:
                match = re.search(r"os\.path\.join\(REPO_ROOT, '([^']+)'", line)
                if match:
                    line = line[:match.end()] + ')' * diff + line[match.end():]
                else:
                    match = re.search(r'os\.path\.join\(REPO_ROOT, "([^"]+)"', line)
                    if match:
                        line = line[:match.end()] + ')' * diff + line[match.end():]
        fixed_lines.append(line)

    content = '\n'.join(fixed_lines)

    if content != original:
        with open(file_path, 'w') as f:
            f.write(content)
        return 1
    return 0


def main():
    parser = argparse.ArgumentParser(description='Find hard-coded paths in the codebase')
    parser.add_argument('--fix', action='store_true', help='Attempt to fix the paths')
    args = parser.parse_args()

    # find the repo root
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(script_dir)
    os.chdir(repo_root)

    print(f"Searching in: {repo_root}")
    print()

    all_findings = []

    # search Python files
    for pattern in ['source/*.py', 'python_scripts/*.py', 'notebooks/*.py']:
        for file_path in glob.glob(pattern):
            findings = find_hardcoded_paths_in_python(file_path)
            all_findings.extend(findings)

    # search notebooks
    for file_path in glob.glob('notebooks/*.ipynb'):
        findings = find_hardcoded_paths_in_notebook(file_path)
        all_findings.extend(findings)

    # report findings
    code_findings = [f for f in all_findings if f['type'] == 'code']
    output_findings = [f for f in all_findings if f['type'] == 'output']

    if code_findings:
        print(f"Found {len(code_findings)} hard-coded path(s) in code:")
        print("-" * 60)
        for f in code_findings:
            if 'line' in f and 'cell' not in f:
                print(f"  {f['file']}:{f['line']}: {f['content']}")
            else:
                print(f"  {f['file']} (cell {f.get('cell', '?')}, line {f.get('line', '?')}): {f['content']}")
        print()
    else:
        print("No hard-coded paths found in code cells/files!")
        print()

    if output_findings:
        print(f"Found {len(output_findings)} hard-coded path(s) in notebook outputs (these are OK):")
        files_with_outputs = set(f['file'] for f in output_findings)
        for f in files_with_outputs:
            print(f"  {f}")
        print()

    # fix if requested
    if args.fix and code_findings:
        print("Attempting to fix Python files...")
        fixed_count = 0
        for pattern in ['source/*.py', 'python_scripts/*.py']:
            for file_path in glob.glob(pattern):
                fixed_count += fix_path_in_python_file(file_path)
        print(f"Fixed {fixed_count} file(s)")
        print("Note: Notebooks need manual fixing or use the fix_notebooks.py script")

    # summary
    if not code_findings:
        print("All code paths have been converted to relative paths!")
        return 0
    else:
        return 1


if __name__ == '__main__':
    exit(main())
