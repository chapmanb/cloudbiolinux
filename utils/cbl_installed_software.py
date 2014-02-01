#!/usr/bin/env python
"""Provide dump of software and libraries installed on CloudBioLinux image.

Run from the top level of the cloudbiolinux source directory:
    python utils/cbl_installed_software.py
"""
import os

from cloudbio import manifest

def main():
    out_dir = os.path.join(os.getcwd(), "manifest")
    manifest.create(out_dir)

if __name__ == "__main__":
    main()
