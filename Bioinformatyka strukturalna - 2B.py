import os
from Bio import PDB
from Bio.PDB import PDBParser
from Bio.PDB import Polypeptide
import numpy as np
import matplotlib.pyplot as plt
 
 
def calculate_phi_psi(structure):
 
    phi_values = []
    psi_values = []
    secondary_structure = []
 
    for model in structure:
        for chain in model:
            polypeptides = Polypeptide.Polypeptide(chain)       # obiekt Polypeptide -> klasa do analizy sekwencji aminokwasów
            phi_psi_angles = polypeptides.get_phi_psi_list()    # pobranie listy kątów phi i psi -> funkcja Polypeptide
            dssp = PDB.DSSP(model, fr)                          # obiekt DSSP -> analiza struktury wtórnej
 
            residues = list(chain.get_residues())                                       # lista reszt aminokwasowych
            for phi, psi in phi_psi_angles:                                             # iteracja po parach kątów phi i psi
                if phi is not None and psi is not None:                                 # czy kąty są zdefiniowane?
                    residue_index = polypeptides.get_phi_psi_list().index((phi, psi))   # indeks reszty aminokwasowej w Polypeptide
                    residue_id = (chain.id, residues[residue_index].get_id()[1])        # identyfikator reszty aminokwasowej -> identyfikator łańcucha + indeks reszty
                    phi_values.append(np.degrees(phi))                                  # dodanie kątów phi w stopniach do list
                    psi_values.append(np.degrees(psi))                                  # dodanie kątów psi w stopniach do list
                    secondary_structure.append(dssp[residue_id][2])                     # dodanie informacji o strukturze wtórnej do listy (trzeci element -> informacja o strukturze wtórnej)
 
    return phi_values, psi_values, secondary_structure
 
 
def plot_ramachandran(phi_values, psi_values, secondary_structure):
 
    plt.figure(figsize=(8, 6))  # rozmiar plot
 
    # zdefiniowanie kolorów dla typów struktury wtórnej
    color_map = {'H': 'red',
                 'B': 'orange',
                 'E': 'yellow',
                 'G': 'green',
                 'I': 'blue',
                 'T': 'purple',
                 'S': 'blue',
                 '-': 'grey'}
 
    # mapowanie kodów struktur wtórnych do ich nazw (na potrzeby legendy)
    full_names = {
        'H': 'Alpha helix',
        'B': 'Beta bridge',
        'E': 'Strand',
        'G': 'Helix-3',
        'I': 'Helix-5',
        'T': 'Turn',
        'S': 'Bend',
        '-': 'Unknown'
    }
 
    handles = []
 
    for label, color in color_map.items():
        handles.append(plt.Line2D([0], [0],
                                  marker='o', color='w', markerfacecolor=color, markersize=10,  # wygląd punktów na legendzie
                                  label=f'{full_names[label]} ({label})'))                      # etykieta dla legendy -> pełna nazwa (skrót)
 
    for phi, psi, ss in zip(phi_values, psi_values, secondary_structure):    # zip -> łączenie w pary elementów
        plt.scatter(phi, psi, s=10, color = color_map[ss])                     # rysowanie pojedyńczych punktów
 
    plt.title('Ramachandran Plot')
    plt.xlabel('Phi Angle (degrees)')
    plt.ylabel('Psi Angle (degrees)')
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.grid(color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
 
    plt.legend(handles=handles, title='Secondary Structure', loc='lower right')  # Adjust the legend position
 
    plt.savefig('ramachandran_plot.pdf')
    plt.show()
 
 
# wpisywanie nazwy pliku
fr = input("PDB FILE NAME: ")
structure = PDBParser(PERMISSIVE=False, QUIET=False).get_structure(fr[-8:-4], os.path.basename(fr))
 
# wywołanie funkcji -> obliczenie kątów i stuktur wtórnych dla podanego pliku PDB
phi_values, psi_values, secondary_structure = calculate_phi_psi(structure)
 
# wywołanie funkcji -> stworzenie wykresu Ramachandrana
plot_ramachandran(phi_values, psi_values, secondary_structure)
 