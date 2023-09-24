# rdkitugm2023_hackathon
code from RDKitUGM Hackathon

## Overview
Python knime node for cheminformatics tasks

## Requirement
Please see ./python_knime_node/knime5_node_dev.yml
Don't use newest version of RDKit! (at 202309)

## Usage
Modify knime.ini in your knime folder and config.yml in this repo.

## Features
Knime node for cheminformatics tasks
 - Enumerate tautomer
 - Enumerate stereoisomer

## Reference
 - https://docs.knime.com/latest/pure_python_node_extensions_guide/#introduction
 - https://github.com/greglandrum/python_knime_nodes/tree/main/basic/tutorial_extension

## Author
 [twitter](https://twitter.com/iwatobipen)
 [bluesky](@iwatobipen.bsky.social)

## Issue
Wedge bond of Enumerated stereoisomers not shown evenif the stereochemistry is enumereted.
