import Bio.PDB
import matplotlib.pyplot as plt
import numpy
import pylab

#Odleglosc pomiedzy weglami alfa par reszt
def distance(res_1, res_2):
    if "CA" in res_1 and "CA" in res_2:
        distance = res_1["CA"].coord - res_2["CA"].coord
        return numpy.sqrt(numpy.sum(distance * distance))


#wypelnianie macierzy odleglosci wartosciami w Da
def matrix(ch_1, ch_2):
    matrix_dist = numpy.zeros((len(ch_1), len(ch_2)), numpy.single)
    for i, res_1 in enumerate(ch_1):
        for j, res_2 in enumerate(ch_2):
            matrix_dist[i, j] = distance(res_1, res_2)

    return matrix_dist

#Odczyt z pliku PDB
id = "1HHB"
file = "1HHB.pdb"
structure = Bio.PDB.PDBParser().get_structure(id, file)
chain = structure[0]

#Przypisywanie wartosci w macierzy odleglosci z 2 argumentami dla lancuchow
matrix_dist = matrix(chain["A"], chain["B"])

#Wypelnienie mapy kontaktow zerami
map = numpy.zeros((len(matrix_dist), len(matrix_dist)), numpy.single)

#Wypelnienie mapy kontaktow 1, jesli spelniony jest warunek odleglosci < 8 Da
for i, a in enumerate(matrix_dist):
    for j, b in enumerate(matrix_dist):
        if matrix_dist[i][j] < 8.0:
            map[i][j] = 1


#Wykres dla mapy kontaktow
plt.figure(figsize=(6, 6))
plt.scatter(*numpy.where(map), color='red', marker='.')
plt.show()

#Wykres dla macierzy odleglosci
pylab.figure(figsize=(6, 6))
pylab.imshow(numpy.transpose(matrix_dist))
pylab.colorbar()
pylab.show()
