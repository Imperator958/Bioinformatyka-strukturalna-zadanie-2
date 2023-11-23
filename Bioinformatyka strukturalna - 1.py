import Bio.PDB
import matplotlib.pyplot as plt
import numpy
import pylab

#Odleglosc pomiedzy weglami alfa par reszt
def residue_distance(res_1, res_2):
    if "CA" in res_1 and "CA" in res_2:
        distance = res_1["CA"].coord - res_2["CA"].coord
        return numpy.sqrt(numpy.sum(distance * distance))


#wypelnianie macierzy odleglosci wartosciami w Da
def fill_distance_matrix(ch_1, ch_2):
    distance_matrix = numpy.zeros((len(ch_1), len(ch_2)), numpy.single)
    for i, res_1 in enumerate(ch_1):
        for j, res_2 in enumerate(ch_2):
            distance_matrix[i, j] = residue_distance(res_1, res_2)

    return distance_matrix

#Odczyt z pliku PDB
id = "1HHB"
file = "1HHB.pdb"
structure = Bio.PDB.PDBParser().get_structure(id, file)
chain = structure[0]

#Przypisywanie wartosci w macierzy odleglosci z 2 argumentami dla lancuchow
distance_matrix = fill_distance_matrix(chain["A"], chain["B"])

#Wypelnienie mapy kontaktow zerami
protein_contact_map = numpy.zeros((len(distance_matrix), len(distance_matrix)), numpy.single)

#Wypelnienie mapy kontaktow 1, jesli spelniony jest warunek odleglosci < 8 Da
for i, a in enumerate(distance_matrix):
    for j, b in enumerate(distance_matrix):
        if distance_matrix[i][j] < 8.0:
            protein_contact_map[i][j] = 1


#Wykres dla mapy kontaktow
plt.figure(figsize=(6, 6))
plt.scatter(*numpy.where(protein_contact_map), color='red', marker='.')
plt.show()

#Wykres dla macierzy odleglosci
pylab.figure(figsize=(6, 6))
pylab.imshow(numpy.transpose(distance_matrix))
pylab.colorbar()
pylab.show()
