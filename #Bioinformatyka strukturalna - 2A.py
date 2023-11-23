#Bioinformatyka strukturalna - 2.A
import os
import numpy as np
from Bio.PDB import calc_dihedral
import Bio.PDB.PDBParser as PDBParser
import Bio.Seq as Seq

#RNA

fr = input("Podaj nazwę pliku w formacie PDB zawierającego RNA: ") #Wpisywanie nazwy pliku
structure = PDBParser(PERMISSIVE = False, QUIET = False).get_structure(fr[-8:-4], os.path.basename(fr)) #Podzielenie sobie struktury na łańcuchy, reszty, atomy itd...

def vector(residue, atom): #Funkcja pobierająca vector współrzędnych konkretnego atomu ze struktury
    if atom in residue:
        return residue[atom].get_vector()
    else:
        None

def calculate_dihedral(vector1, vector2, vector3, vector4): #Funkcja licząca kąt pomiędzy 4 atomami
    if vector1 is not None and vector2 is not None and vector3 is not None and vector4 is not None:
        return calc_dihedral(vector1, vector2, vector3, vector4)
    else:
        return "NULL"

headers = ["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "chi"] #Nazwy kątów

with open("RNAOut.txt", "w") as f: #Zapisywanie bezpośrednio do pliku otrzymanych kątów
    f.write("|".join(headers) + "\n")
    all_residues = structure.get_residues()
    i = 0
    residues = list(all_residues)
    while i < len(residues): #iterowanie po wszystkich resztach
        residue_name = residues[i].get_resname()
        if(residue_name == "HOH"): #Jeżeli woda to skip
            i+=1
            continue
        if i > 0: #Jeżeli i == 0 to znaczy, że jest to pierwsz reszta i nie ma żadnej "poprzedniej", więc atom z poprzedniej reszty nie może istnieć i ma ustawioną wartość na None
            p_residue = residues[i-1]
            if "O3'" in p_residue: #Jeżeli w poprzedniej reszcie jest atom O3'
                O3P = p_residue["O3'"].get_vector()
            else:
                O3P = None
        else:
            O3P = None
        if i < len(residues) - 1: #Jeżeli i doszło już do ostatniej reszty to niemożliwe jest pobranie dwóch atomów z następnych reszt, więc będą miały wartość None
            n_residue = residues[i+1]
            if "O5'" in n_residue: #Sprawdzanie czy 05' jest w kolejnej reszcie
                O5N = n_residue["O5'"].get_vector()
            else: 
                O5N = None
            if "P" in n_residue: #Sprawdzanie czy P jest w kolejnej reszcie
                PN = n_residue["P"].get_vector()
            else:
                PN = None
        else:
            O5N = None
            PN = None

        #przyspisywanie konkretnym wektorom współrzędnych danego atomu

        P = vector(residues[i], "P")
        O5 = vector(residues[i], "O5'")
        C5 = vector(residues[i], "C5'")
        C4 = vector(residues[i], "C4'")
        C3 = vector(residues[i], "C3'")
        O3 = vector(residues[i], "O3'")
        C1 = vector(residues[i], "C1'")
        O4 = vector(residues[i], "O4'")
        C2C = vector(residues[i], "C2")
        C4C = vector(residues[i], "C4")
        if(residue_name == "A" or residue_name == "G"):
            N = vector(residues[i], "N9")
        else:
            N = vector(residues[i], "N1")

        
        #Obliczanie kątów na podstawie współrzędnych atomów
        alpha = calculate_dihedral(O3P, P, O5, C5)
        beta = calculate_dihedral(P, O5, C5, C4)
        gamma = calculate_dihedral(O5, C5, C4, C3)
        delta = calculate_dihedral(C5, C4, C3, O3)
        epsilon = calculate_dihedral(C4, C3, O3, PN)
        zeta = calculate_dihedral(C3, O3, PN, O5N)
        if(residue_name == "A" or residue_name == "G"):
            chi = calculate_dihedral(O4, C1, N, C4C)
        else:
            chi = calculate_dihedral(O4, C1, N, C2C)

        angles = [alpha, beta, gamma, delta, epsilon, zeta, chi] #Nazwy kątów
        f.write(" ".join(map(str, angles)) + "\n") #Zapis do pliku obliczonych kątów torsyjnych
        i+=1

