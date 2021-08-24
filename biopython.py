
'''
Fuente de consulta: http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec118'''

# Librerías
## biopython librería para realizar cálculos bioinformáticos
## Se utiliza para datos biológicos, principalmente  de ADN y de proteínas
!pip install biopython
import Bio
## Creación de secuencias
from Bio.Seq import Seq
## Lectura de archivos de secuencia SeqIO (sequence Input/Output)
from Bio import SeqIO
## La versión en línea de BLAST
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
### Alineamiento entre dos secuencias
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
## Búsqueda de enzimas de restricción
from Bio.Restriction import *
from Bio.Restriction.PrintFormat import PrintFormat
## Construcción de un árbol filogenético
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO

# Creación de una secuencia para analizar
'''
En este caso se ha consultado la secuencia del gen El gen de la neomicina
fosfotransferasa (nptII), procedente del transposon Tn5, el cual confiere
resistencia a antibióticos aminoglucósidos, como la kanamicina.
La enzima desactiva la sustancia antibiótica mediante fosforilación.
Se utiliza generalmente como marcador de selección.
'''

'''
>ENA|AAB00444|AAB00444.1 Plasmid pRL1063a neomycin phosphotransferase
ATGATTGAACAAGATGGATTGCACGCAGGTTCTCCGGCCGCTTGGGTGGAGAGGCTATTC
GGCTATGACTGGGCACAACAGACAATCGGCTGCTCTGATGCCGCCGTGTTCCGGCTGTCA
GCGCAGGGGCGCCCGGTTCTTTTTGTCAAGACCGACCTGTCCGGTGCCCTGAATGAACTG
CAGGACGAGGCAGCGCGGCTATCGTGGCTGGCCACGACGGGCGTTCCTTGCGCAGCTGTG
CTCGACGTTGTCACTGAAGCGGGAAGGGACTGGCTGCTATTGGGCGAAGTGCCGGGGCAG
GATCTCCTGTCATCTCACCTTGCTCCTGCCGAGAAAGTATCCATCATGGCTGATGCAATG
CGGCGGCTGCATACGCTTGATCCGGCTACCTGCCCATTCGACCACCAAGCGAAACATCGC
ATCGAGCGAGCACGTACTCGGATGGAAGCCGGTCTTGTCGATCAGGATGATCTGGACGAA
GAGCATCAGGGGCTCGCGCCAGCCGAACTGTTCGCCAGGCTCAAGGCGCGCATGCCCGAC
GGCGAGGATCTCGTCGTGACCCATGGCGATGCCTGCTTGCCGAATATCATGGTGGAAAAT
GGCCGCTTTTCTGGATTCATCGACTGTGGCCGGCTGGGTGTGGCGGACCGCTATCAGGAC
ATAGCGTTGGCTACCCGTGATATTGCTGAAGAGCTTGGCGGCGAATGGGCTGACCGCTTC
CTCGTGCTTTACGGTATCGCCGCTCCCGATTCGCAGCGCATCGCCTTCTATCGCCTTCTT
GACGAGTTCTTCTGA
'''
## Nota: si la cadena tiene varias líneas hay que utilizar la triple comilla.

nptII=Seq('''ATGATTGAACAAGATGGATTGCACGCAGGTTCTCCGGCCGCTTGGGTGGAGAGGCTATTC
           GGCTATGACTGGGCACAACAGACAATCGGCTGCTCTGATGCCGCCGTGTTCCGGCTGTCA
           GCGCAGGGGCGCCCGGTTCTTTTTGTCAAGACCGACCTGTCCGGTGCCCTGAATGAACTG
           CAGGACGAGGCAGCGCGGCTATCGTGGCTGGCCACGACGGGCGTTCCTTGCGCAGCTGTG
           CTCGACGTTGTCACTGAAGCGGGAAGGGACTGGCTGCTATTGGGCGAAGTGCCGGGGCAG
           GATCTCCTGTCATCTCACCTTGCTCCTGCCGAGAAAGTATCCATCATGGCTGATGCAATG
           CGGCGGCTGCATACGCTTGATCCGGCTACCTGCCCATTCGACCACCAAGCGAAACATCGC
           ATCGAGCGAGCACGTACTCGGATGGAAGCCGGTCTTGTCGATCAGGATGATCTGGACGAA
           GAGCATCAGGGGCTCGCGCCAGCCGAACTGTTCGCCAGGCTCAAGGCGCGCATGCCCGAC
           GGCGAGGATCTCGTCGTGACCCATGGCGATGCCTGCTTGCCGAATATCATGGTGGAAAAT
           GGCCGCTTTTCTGGATTCATCGACTGTGGCCGGCTGGGTGTGGCGGACCGCTATCAGGAC
           ATAGCGTTGGCTACCCGTGATATTGCTGAAGAGCTTGGCGGCGAATGGGCTGACCGCTTC
           CTCGTGCTTTACGGTATCGCCGCTCCCGATTCGCAGCGCATCGCCTTCTATCGCCTTCTT
           GACGAGTTCTTCTGA''')

nptII

# Estimación de la longitud en nucleótidos del gen nptII
'''nt=nucleótidos'''
print(f'longitud de la secuencia es: {len(nptII)} nt')

# Obtención de la secuencia reversa complementaria

print(f'La secuencia reversa complementaria del gen nptII es: \n {nptII.reverse_complement()}')

# Contar el número de diferentes nucleótidos

nucleotide=['A', 'T', 'C', 'G']
name=['Adenina', 'Timina', 'Citosina', 'Guanina']
for i, nt in enumerate(nucleotide):
    print(f'Número de nucleótidos correspondientes a {name[i]}= {nptII.count(nt)}')

# Búsqueda de aminoácidos en la secuencia
'''
Se puede buscar si existe en la secuencia algún codón de terminación, el cual
no determina ningún aminoácido y su función es simplemente la de acotar el
mensaje cifrado por el ADN. Existen varios codones de terminación "TAG", "TGA"
y "TAA"
'''
codon_stop = ['TAG', 'TGA', 'TAA']

for index in range(0, len(nptII)+1, 3):
    codon=nptII[index:index+3]
    if codon in codon_stop:
        print(codon)
        print(index,index+3)

'''
Se puede observar que el codon stop se corresponde con los últimos tres nucleótidos
de la secuencia del gen.
'''

# Lectura de un archivo de secuencia


for record in SeqIO.parse('AAB00444.1.fasta', 'fasta'):
    print(record.id)
    print(record.seq)
    print(record.description)

'''Se puede observar como utilizando la función SeqIO.parse, se puede obtener toda la información
de la sequencia a través de un iterador que proporciona objetos de SeqRecord'''

# Búsqueda de enzimas de restricción en la secuencia
## Se analiza la secuencia completa para todas las enzimas
a=Analysis(AllEnzymes, nptII)

a.print_that()

## Búsqueda de sitios de restricción de algunas enzimas en la secuencia nptII
restriction_batch = RestrictionBatch(['AccII', 'AfaI', 'RsaI'])
result=restriction_batch.search(nptII)
result

### Visualización de los sitios de restricción

my_map=PrintFormat()
my_map.sequence=nptII
my_map.print_as("map")
my_map.print_that(result)


# Alineamiento de secuencias

## Lectura de la secuencia
record=SeqIO.read('AAB00444.1.fasta', 'fasta')
record

## blast de la secuencia
result=NCBIWWW.qblast('blastn', 'nt', record.format('fasta'))
## Guardar el resultado en un archivo
with open('blast_result.xml', 'w') as out_put:
    out_put.write(result.read())
result.close()

'''
Al abrir el resultado, se puede observar como aparecen muchos plásmido que en su secuencias
contienen este gen de resistencia a antibióticos.
'''
## Lectura del blastn
result=open('blast_result.xml')

blast_result=NCBIXML.read(result)

## Obtención de los alineamientos

print(f'Alineamientos para la secuencia {blast_result.query}')

for alignment in blast_result.alignments:
    print(f'Número de acceso: {alignment.accession}')
    print(f'Nombre del plásmido o secuencia: {alignment.title}')
    print(f'Longitud: {alignment.length}')
    print()

'''
Para analizar el resultado de varias consultas:
'''
result=open('blast_result.xml')
blast_records=NCBIXML.parse(result)
blast_records
for blast_record in blast_records:
    print(f'Resultado del BLAST: {blast_record.query}')
    print(f'Number of alignments: {len(blast_record.alignments)}')
    print()



# Alineamiento global entre dos secuencias
'''
El método "next()" crea un iterador que devuelve uno por uno cada uno de los
elementos
'''

seq1=next(SeqIO.parse('AAB00444.1.fasta', format='fasta'))
seq2=next(SeqIO.parse('CAA23656.1.fasta', format='fasta'))

alignments=pairwise2.align.globalxx(seq1, seq2)

for a in alignments:
    print(format_alignment(*a))



'''
Como se ha visto, se ha realizado un alineamiento global con una puntuación (score) predeterminado.
No obstante, también se puede personalizar el esquema de puntuación
(asignar valores personalizados para coincidencias, discrepancias y brechas).
En este caso, para coincidentes se asignan 2 puntos, se deduce 1 punto por cada carácter no coincidente
y se deducen 0,5 puntos al abrir un hueco y 0,1 puntos al ampliarlo.
'''

alignments=pairwise2.align.globalms(seq1, seq2, 2, -1, -0.5, -0.1)

for a in alignments:
    print(format_alignment(*a))

# Construcción de un árbol filogenético

## Lectura del archivo
aln=AlignIO.read('homo.phy', 'phylip')

## Obtención de la matriz de distancias

calculator=DistanceCalculator('identity')
dm=calculator.get_distance(aln)

## Construcción del árbol filogenético
constructor=DistanceTreeConstructor()
tree=constructor.upgma(dm)

## Representación del árbol filogenético

Phylo.draw_ascii(tree)
'''
Ejemplificación de como, a través de un archivo de frecuencias de distintas especies,
se puede crear un árbol filogenético
'''
