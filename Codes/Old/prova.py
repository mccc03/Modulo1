Directory = '/home/marco/Documents/University/Physics/MetodiNumerici/Modulo1/_data/'
out_file = open(Directory+'prova.txt', "w+")

a=[[1,2],[1,2,3,4]]

out_file.write(str(a))
out_file.close()