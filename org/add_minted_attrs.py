#!/usr/bin/env python

import re


def get_default_tangle_file(lines):
    """
    Extracts the default tangle file name from the #+PROPERTY: line.
    """
    for line in lines:
        match = re.match(
            r"^\s*#\+PROPERTY:.*:tangle\s+([^\s]+)", line, re.IGNORECASE
        )
        if match:
            tangle_file = match.group(1)
            return tangle_file
    return None


def remove_attr_option_strings(fname):
    with open(fname, 'r') as file:
        lines = file.readlines()

    with open(fname, 'w') as file:
        for line in lines:
            if "#+attr_latex: :options label=" in line.lower():
                continue
            file.write(line)


def process_org_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    default_tangle_file = get_default_tangle_file(lines)
    if not default_tangle_file:
        print(
            "The default  file to which to tangle to is not set in the propertiess."
        )
        quit()

    output_lines = []
    output_lines.append(lines[0])
    i = 1  # skip first line of org mode file

    while i < len(lines):
        line = lines[i]
        match = re.match(r"^\s*#\+begin_src sage(.*)", line)
        if match:
            if ":tangle" not in match.group(1):
                tangle_fname = f"{default_tangle_file}"
            else:
                match = re.match(r".*:tangle\s+([^\s]+)", line)
                tangle_fname = match.group(1)
                if tangle_fname == "no":
                    tangle_fname = "don't tangle"
            tangle_fname = tangle_fname.replace("_", r"\_")
            attr_line = f"#+attr_latex: :options label={tangle_fname}\n"
            output_lines.append(attr_line)

        output_lines.append(line)
        i += 1

    with open(filename, 'w') as file:
        file.writelines(output_lines)


files = [
    "show_expression.org",
    "tuples.org",
    "functions.org",
    "differentiation.org",
    "section1.4.org",
    "section1.5.org",
    "section1.6.org",
    "section1.7.org",
    "section1.8.org",
    "section1.9.org",
    "section3.1.org",
    "section3.2.org",
    "section3.4.org",
    "section3.5.org",
    "section3.9.org",
    "section5.1.org",
    "section5.2.org",
]


if __name__ == "__main__":
    for fname in files:
        remove_attr_option_strings(fname)
        process_org_file(fname)
