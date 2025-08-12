import random
import time
import matplotlib.pyplot as plt

rna_to_amino = {
    'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
    'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
    'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
    'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
    'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
    'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
    'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
    'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
    'UAC': 'Y', 'UAU': 'Y', 'UAA': '_', 'UAG': '_',
    'UGC': 'C', 'UGU': 'C', 'UGA': '_', 'UGG': 'W',
}

sample_genes = {
    "Insulin": "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTC",
    "Hemoglobin": "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAG",
    "BRCA1": "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTG",
    "p53": "ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAA",
    "CFTR": "ATGCAGAGGTCGCCTCTGGAAAAGGCCAGCGTTGTCTCCAAACTTTTTTTCAGCTGGACCAGACCAATTTTGAGGAAAGGATACAGACAGCGCCTGGAATTGTCAGACATATACCA"
}

mutation_stats = {
    "types": [],  
    "positions": [],  
    "protein_effects": [],  
    "amino_changes": []  
}

def make_random_dna(length=100):
    bases = ['A', 'T', 'G', 'C']
    dna = "ATG"
    
    for i in range(length - 3):
        dna = dna + random.choice(bases)
        
    return dna

def dna_to_rna(dna):
    rna = ""
    for base in dna:
        if base == 'T':
            rna = rna + 'U'
        else:
            rna = rna + base
    return rna

def dna_to_protein(dna):
    dna = dna.upper()
    
    rna = dna_to_rna(dna)
    
    start = rna.find('AUG')
    
    if start == -1:
        return "No start codon found"
    
    protein = ""
    
    i = start
    while i < len(rna) - 2:
        codon = rna[i:i+3]
        
        if len(codon) < 3:
            break
            
        if codon in rna_to_amino:
            amino = rna_to_amino[codon]
        else:
            amino = 'X'
        
        if amino == '_':
            break
            
        protein = protein + amino
        i = i + 3
    
    return protein

def create_mutation(dna):
    types = ["point", "insertion", "deletion"]
    mutation_type = random.choice(types)
    
    original_dna = dna
    
    if mutation_type == "point":
        position = random.randint(0, len(dna) - 1)
        old_base = dna[position]
        
        bases = ['A', 'T', 'G', 'C']
        bases.remove(old_base)
        new_base = random.choice(bases)
        
        mutated_dna = dna[:position] + new_base + dna[position+1:]
        description = f"Changed base at position {position} from {old_base} to {new_base}"
    
    elif mutation_type == "insertion":
        position = random.randint(0, len(dna))
        
        new_bases = ""
        for i in range(random.randint(1, 3)):
            new_bases = new_bases + random.choice(['A', 'T', 'G', 'C'])
        
        mutated_dna = dna[:position] + new_bases + dna[position:]
        description = f"Added {new_bases} at position {position}"
    
    else:  # deletion
        position = random.randint(0, len(dna) - 1)
        
        num_to_delete = min(3, len(dna) - position)
        num_to_delete = random.randint(1, num_to_delete)
        
        deleted = dna[position:position+num_to_delete]
        mutated_dna = dna[:position] + dna[position+num_to_delete:]
        description = f"Removed {deleted} at position {position}"
    
    old_protein = dna_to_protein(original_dna)
    new_protein = dna_to_protein(mutated_dna)
    
    old_rna = dna_to_rna(original_dna)
    new_rna = dna_to_rna(mutated_dna)
    
    mutation_stats["types"].append(mutation_type)
    mutation_stats["positions"].append(position)
    
    if old_protein == new_protein:
        effect = "unchanged"
    elif len(new_protein) < len(old_protein):
        effect = "shorter"
    elif len(new_protein) > len(old_protein):
        effect = "longer"
    else:
        effect = "changed"
    
    mutation_stats["protein_effects"].append(effect)
    
    changes = 0
    for i in range(min(len(old_protein), len(new_protein))):
        if old_protein[i] != new_protein[i]:
            changes += 1
    
    mutation_stats["amino_changes"].append(changes)
    
    return {
        "original_dna": original_dna,
        "mutated_dna": mutated_dna,
        "original_rna": old_rna,
        "mutated_rna": new_rna,
        "description": description,
        "mutation_type": mutation_type,
        "original_protein": old_protein,
        "mutated_protein": new_protein
    }

def format_dna(dna, width=50):
    result = ""
    
    for i in range(0, len(dna), width):
        result = result + dna[i:i+width] + "\n"
        
    return result

def show_mutation(mutation_data):
    print("\n=== " + mutation_data['description'] + " ===\n")
    
    print("Original DNA:")
    print(format_dna(mutation_data['original_dna']))
    print("Mutated DNA:")
    print(format_dna(mutation_data['mutated_dna']))
    
    print("Original RNA:")
    print(format_dna(mutation_data['original_rna']))
    print("Mutated RNA:")
    print(format_dna(mutation_data['mutated_rna']))
    
    print("Original protein:")
    print(mutation_data['original_protein'])
    print("Mutated protein:")
    print(mutation_data['mutated_protein'])
    
    orig = mutation_data['original_protein']
    mutated = mutation_data['mutated_protein']
    
    if orig == mutated:
        print("RESULT: No change to protein")
    elif len(mutated) < len(orig):
        print("RESULT: Protein got shorter")
    elif len(mutated) > len(orig):
        print("RESULT: Protein got longer")
    else:
        changes = 0
        for i in range(len(orig)):
            if i < len(mutated) and orig[i] != mutated[i]:
                changes = changes + 1
                print(f"Change at position {i+1}: {orig[i]} changed to {mutated[i]}")
        
        print(f"RESULT: {changes} amino acid changes")

def test_random_mutations(count=3):
    for key in mutation_stats:
        mutation_stats[key] = []
    
    try:
        gene_names = list(sample_genes.keys())
        gene_name = random.choice(gene_names)
        start_dna = sample_genes[gene_name]
        
        print(f"Starting with {gene_name} gene ({len(start_dna)} bases)")
        print(format_dna(start_dna))
        print(f"RNA: {dna_to_rna(start_dna)[:50]}...")
        print(f"Original protein: {dna_to_protein(start_dna)}")
        
        current_dna = start_dna
        
        for i in range(count):
            print(f"\n--- Mutation {i+1} ---")
            mutation = create_mutation(current_dna)
            show_mutation(mutation)
            
            current_dna = mutation['mutated_dna']
            
            if i < count - 1:
                print("Applying next mutation...")
                time.sleep(1)
        
        visualize_mutations()
    except (EOFError, KeyboardInterrupt):
        print("\n\nProgram terminated by user.")
        return

def mutate_specific_gene():
    for key in mutation_stats:
        mutation_stats[key] = []
    
    print("\nGenes available:")
    
    gene_names = list(sample_genes.keys())
    for i in range(len(gene_names)):
        gene_name = gene_names[i]
        print(f"{i+1}. {gene_name} ({len(sample_genes[gene_name])} bases)")
    
    valid_choice = False
    while not valid_choice:
        try:
            choice = input("\nChoose a gene (1-5) or type 'custom' for your own DNA: ")
            
            if choice.lower() == "custom":
                valid_choice = True
                dna = input("Enter DNA sequence (or press Enter for random): ").upper()
                
                if dna == "":
                    valid_length = False
                    while not valid_length:
                        try:
                            length_input = input("How long? (default 100, range 10-1000): ")
                            if length_input == "":
                                length = 100
                                valid_length = True
                            else:
                                try:
                                    length = int(length_input)
                                    if 10 <= length <= 1000:
                                        valid_length = True
                                    else:
                                        print("Please enter a length between 10 and 1000.")
                                except ValueError:
                                    print("Please enter a valid number.")
                        except (EOFError, KeyboardInterrupt):
                            print("\n\nProgram terminated by user.")
                            return
                        
                    dna = make_random_dna(length)
                else:
                    clean_dna = ""
                    for c in dna:
                        if c in "ATGC":
                            clean_dna = clean_dna + c
                    
                    dna = clean_dna
                    
                    if not dna.startswith("ATG"):
                        dna = "ATG" + dna
            else:
                try:
                    index = int(choice) - 1
                    if 0 <= index < len(gene_names):
                        dna = sample_genes[gene_names[index]]
                        valid_choice = True
                    else:
                        print(f"Please enter a number between 1 and {len(gene_names)} or 'custom'.")
                except ValueError:
                    print("Please enter a valid number or 'custom'.")
        except (EOFError, KeyboardInterrupt):
            print("\n\nProgram terminated by user.")
            return
    
    valid_num = False
    while not valid_num:
        try:
            num_input = input("\nHow many mutations? (1-10): ")
            try:
                num = int(num_input)
                if 1 <= num <= 10:
                    valid_num = True
                else:
                    print("Please enter a number between 1 and 10.")
            except ValueError:
                print("Please enter a valid number.")
        except (EOFError, KeyboardInterrupt):
            print("\n\nProgram terminated by user.")
            return
    
    current_dna = dna
    
    print(f"\nStarting DNA ({len(current_dna)} bases):")
    print(format_dna(current_dna))
    print(f"Starting RNA: {dna_to_rna(current_dna)[:50]}...")
    print(f"Starting protein: {dna_to_protein(current_dna)}")
    
    for i in range(num):
        print(f"\n--- Mutation {i+1} ---")
        mutation = create_mutation(current_dna)
        show_mutation(mutation)
        
        current_dna = mutation['mutated_dna']
        
        if i < num - 1:
            print("Applying next mutation...")
            time.sleep(1)
    
    visualize_mutations()

def visualize_mutations():
    if len(mutation_stats["types"]) == 0:
        print("No mutation data to visualize.")
        return
    
    print("\n--- Mutation Statistics Visualization ---")
    
    try:
        plt.figure(figsize=(14, 10))
        
        # Plot 1: Mutation Types
        plt.subplot(2, 2, 1)
        types_count = {"point": 0, "insertion": 0, "deletion": 0}
        for t in mutation_stats["types"]:
            types_count[t] += 1
        
        plt.bar(types_count.keys(), types_count.values(), color=['blue', 'green', 'red'])
        plt.title('Types of Mutations', pad=15)
        plt.ylabel('Count')
        
        # Plot 2: Mutation Positions
        plt.subplot(2, 2, 2)
        plt.hist(mutation_stats["positions"], bins=10, color='orange')
        plt.title('Mutation Positions', pad=15)
        plt.xlabel('Position in DNA')
        plt.ylabel('Frequency')
        
        # Plot 3: Protein Effects
        plt.subplot(2, 2, 3)
        effects_count = {"unchanged": 0, "shorter": 0, "longer": 0, "changed": 0}
        for e in mutation_stats["protein_effects"]:
            effects_count[e] += 1
        
        plt.bar(effects_count.keys(), effects_count.values(), color=['green', 'red', 'blue', 'purple'])
        plt.title('Effects on Protein', pad=15)
        plt.ylabel('Count')
        
        # Plot 4: Amino Acid Changes
        plt.subplot(2, 2, 4)
        if len(mutation_stats["amino_changes"]) == 0:
            plt.text(0.5, 0.5, "No amino acid changes to display", 
                    horizontalalignment='center', verticalalignment='center')
            plt.title('Number of Amino Acid Changes', pad=15)
        else:
            max_change = max(mutation_stats["amino_changes"]) + 1
            plt.hist(mutation_stats["amino_changes"], bins=range(max_change + 1), color='purple')
            plt.title('Number of Amino Acid Changes', pad=15)
            plt.xlabel('Number of Changes')
            plt.ylabel('Frequency')
        
        plt.tight_layout(pad=4.0)
        plt.show()
        
    except Exception as e:
        print(f"Could not create visualization: {e}")
        print("You may need to install matplotlib or run this in an environment that supports plotting.")

def main():
    print("\n" + "=" * 50)
    print(" DNA MUTATION SIMULATOR ")
    print("=" * 50)
    print("\nThis program shows what happens when DNA changes")
    print("It converts DNA → RNA → Protein")
    
    try:
        while True:
            print("\n" + "-" * 30)
            print("MENU")
            print("-" * 30)
            print("1. Test random mutations")
            print("2. Pick a gene to mutate")
            print("3. Exit")
            
            try:
                choice = input("\nWhat do you want to do? (1-3): ")
                
                if choice == "1":
                    test_random_mutations()
                elif choice == "2":
                    mutate_specific_gene()
                elif choice == "3":
                    print("\nThanks for using the DNA Mutation Simulator!")
                    break
                else:
                    print("Please choose 1, 2, or 3.")
            except (EOFError, KeyboardInterrupt):
                print("\n\nProgram terminated by user.")
                break
    except (EOFError, KeyboardInterrupt):
        print("\n\nProgram terminated by user.")

random.seed(time.time())
main() 