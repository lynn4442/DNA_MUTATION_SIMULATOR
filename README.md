# DNA Mutation Simulator

An educational program that simulates DNA mutations and demonstrates their effects on RNA and protein synthesis, following the Central Dogma of Molecular Biology (DNA → RNA → Protein).

## Table of Contents
- [Overview](#overview)
- [Features](#features)
- [How It Works](#how-it-works)
- [Usage](#usage)
- [Program Components](#program-components)
- [Mutation Types](#mutation-types)
- [Visualization](#visualization)
- [Educational Value](#educational-value)
- [Technical Details](#technical-details)

## Overview

This program is designed for educational purposes to help students and enthusiasts understand basic concepts in genetics and molecular biology. It simulates how DNA mutations occur and demonstrates their cascading effects through the biological process of protein synthesis.

The program follows the real biological pathway:
1. **DNA** → **RNA** (Transcription)
2. **RNA** → **Protein** (Translation)

## Features

### Core Functionality
- ✨ Simulates three basic types of DNA mutations
- 🧬 Demonstrates transcription (DNA to RNA) and translation (RNA to protein)
- 📊 Provides comprehensive data visualization
- 🎯 Interactive menu system for user control
- 📚 Includes sample genes for common proteins (Insulin, Hemoglobin, etc.)
- 🔄 Shows complete molecular flow from DNA changes to protein changes

### Visualization Options
- **Graphical Mode**: Beautiful matplotlib histograms (when available)
- **Text Mode**: ASCII bar charts for universal compatibility
- **Statistics**: Comprehensive mutation pattern analysis

### Educational Benefits
- Demonstrates the Central Dogma of Molecular Biology
- Shows how mutations can be silent, harmful, or beneficial
- Visualizes statistical patterns in mutation data
- Accessible to users regardless of technical setup
- Helps understand concepts like truncated proteins and amino acid substitutions

## How It Works

### What is DNA?
DNA is like a long instruction book made up of 4 letters: **A**, **T**, **G**, and **C**. These letters are arranged in groups of 3 (called **codons**), and each group tells the cell to make a specific amino acid. A series of amino acids forms a **protein**, which performs important functions in your body.

### What is RNA?
RNA is a copy of DNA that the cell uses to make proteins. When DNA is copied to RNA, the letter **T** (thymine) changes to **U** (uracil). So RNA has the letters **A**, **U**, **G**, and **C**.

### The Process
1. Creates DNA sequences (random or from sample genes)
2. Applies random mutations to the DNA
3. Transcribes DNA to RNA (T → U conversion)
4. Translates RNA to protein using the genetic code
5. Compares original and mutated sequences
6. Visualizes mutation statistics and effects

## Usage

Run the program:
```bash
python dna_simulator.py
```

### Menu Options
1. **Test Random Mutations**: Apply mutations to randomly generated DNA
2. **Select Specific Gene**: Choose from predefined genes (Insulin, Hemoglobin, etc.)
3. **Custom Settings**: Adjust number of mutations and visualization preferences

### Sample Output
```
Original DNA: ATGAAGTTTGGCGAA...
Mutated DNA:  ATGAAGTCTGGCGAA...
Original RNA: AUGAAGUUUGGCGAA...
Mutated RNA:  AUGAAGUCUGGCGAA...
Original Protein: M-K-F-G-E...
Mutated Protein:  M-K-S-G-E...
Effect: Amino acid changed at position 3 (F → S)
```

## Program Components

### Key Functions

#### `dna_to_rna(dna_sequence)`
Converts DNA to RNA by replacing all T's with U's (transcription process).

#### `dna_to_protein(dna_sequence)`
Converts DNA to protein through a two-step process:
1. DNA → RNA (transcription)
2. RNA → Protein (translation using genetic code)

#### `create_mutation(dna_sequence)`
Applies random mutations and tracks statistics:
- Selects mutation type and position
- Applies the mutation
- Compares original and mutated proteins
- Records detailed statistics

#### `visualize_mutations(mutation_stats)`
Creates comprehensive visualizations:
- Mutation type distribution
- Position frequency analysis
- Protein effect categorization
- Amino acid change statistics

#### `text_bar_chart(data, title)`
Generates ASCII bar charts for environments without matplotlib.

### Data Structures

#### `rna_to_amino` Dictionary
Translation table converting RNA codons to amino acids:
```python
'AUG': 'M',  # Methionine (Start codon)
'UUU': 'F',  # Phenylalanine
'UAG': '_',  # Stop codon
# ... and more
```

## Mutation Types

### 1. Point Mutations
- **Description**: Single base changes (e.g., A → G)
- **Effect**: May cause amino acid substitution or be silent
- **Example**: `ATG` → `AGG` (Start codon → Arginine)

### 2. Insertions
- **Description**: Addition of new bases
- **Effect**: Often causes frameshift, changing all downstream amino acids
- **Example**: `ATGAAA` → `ATGCAAA` (adds C)

### 3. Deletions
- **Description**: Removal of existing bases
- **Effect**: Often causes frameshift or protein truncation
- **Example**: `ATGAAA` → `ATAAA` (removes G)

## Visualization

### Graphical Mode (with matplotlib)
Four comprehensive histograms showing:
1. **Mutation Types**: Distribution of point, insertion, and deletion mutations
2. **Mutation Positions**: Frequency of mutations at different DNA positions
3. **Protein Effects**: Categorized impact on protein structure
4. **Amino Acid Changes**: Number of amino acid substitutions per mutation

## Educational Value

### Concepts Demonstrated
- **Central Dogma**: DNA → RNA → Protein flow
- **Transcription**: DNA to RNA conversion
- **Translation**: RNA to protein synthesis
- **Genetic Code**: Codon to amino acid mapping
- **Mutation Effects**: Silent, missense, nonsense, and frameshift mutations

### Learning Outcomes
- Understand how DNA mutations propagate through molecular processes
- Visualize statistical patterns in genetic variation
- Comprehend the relationship between genotype and phenotype
- Appreciate the complexity of molecular biology

## Technical Details

### Code Quality
- **Simplicity**: Uses basic Python concepts for accessibility
- **Clarity**: Prioritizes understanding over optimization
- **Modularity**: Well-organized functions with clear purposes
- **Robustness**: Graceful handling of missing dependencies
- **Compatibility**: Works across different Python environments

### Performance
- Suitable for educational use with moderate-sized sequences
- Efficient for typical classroom demonstrations
- Scalable for different numbers of mutations

### Development Setup
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

*This program is designed for educational purposes to help understand the fascinating world of genetics and molecular biology.*
