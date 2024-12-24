import turtle
#Celina Sid Abdelkader, anas mranialaoui 
#TP_01, 2024-03-28

"""
Ce projet sert à trancrire, à partir du brin principal d'une séquence d'adn fournie, les gènes existants 
sur le brin principal et complémentaire puis de dessiner les gènes codés par ce brin
"""

##########################################################################################
#Dictionnaire contenant l'équivalant en acide aminé de chaque codon (triplet bases azotée)
##########################################################################################

'''
adn="TCGACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGATACCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATACCCAGCCAGCCAGCCAGCGACG\
GCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGAGTGCCAGCCAGCCAGCCAGCGAACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGAATGCCAGCCAGC\
CAGCCAGCGAAGCCAGCCAGCCGATATTCAGCCAGCCAGCCAGCGAACACTCTTCGACAGCCAGCGAAGCCAGCCAGCCGATATTCAGCCAGCCAGCCAGCGA\
ACTCGACACTCTTCGACAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCATCCCAGCGATACCC\
AGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAG\
CCAGCGAACTCGTCTGCGTTCGACAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGC\
GATTGCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGTATGCCAGCC\
AGCATCCCAGCGA"
'''
#Dictionnaire contenant l'équivalant en mot de chaque codon (triplet bases azotée == mot)
codons_aa = {
    "UUU": "Phénylalanine",
    "UUC": "Phénylalanine",
    "UUA": "Leucine",
    "UUG": "Leucine",
    "CUU": "Leucine",
    "CUC": "Leucine",
    "CUA": "Leucine",
    "CUG": "Leucine",
    "AUU": "Isoleucine",
    "AUC": "Isoleucine",
    "AUA": "Isoleucine",
    "AUG": "Méthionine (Start)",
    "GUU": "Valine",
    "GUC": "Valine",
    "GUA": "Valine",
    "GUG": "Valine",
    "UCU": "Sérine",
    "UCC": "Sérine",
    "UCA": "Sérine",
    "UCG": "Sérine",
    "CCU": "Proline",
    "CCC": "Proline",
    "CCA": "Proline",
    "CCG": "Proline",
    "ACU": "Thrénine",
    "ACC": "Thrénine",
    "ACA": "Thrénine",
    "ACG": "Thrénine",
    "GCU": "Alanine",
    "GCC": "Alanine",
    "GCA": "Alanine",
    "GCG": "Alanine",
    "UAU": "Tyrosine",
    "UAC": "Tyrosine",
    "UAA": "Stop",
    "UAG": "Stop",
    "CAU": "Histidine",
    "CAC": "Histidine",
    "CAA": "Glutamine",
    "CAG": "Glutamine",
    "AAU": "Asparagine",
    "AAC": "Asparagine",
    "AAA": "Lysine",
    "AAG": "Lysine",
    "GAU": "Aspartate",
    "GAC": "Aspartate",
    "GAA": "Glutamate",
    "GAG": "Glutamate",
    "UGU": "Cystéine",
    "UGC": "Cystéine",
    "UGA": "Stop",
    "UGG": "Tryptophane",
    "CGU": "Arginine",
    "CGC": "Arginine",
    "CGA": "Arginine",
    "CGG": "Arginine",
    "AGU": "Sérine",
    "AGC": "Sérine",
    "AGA": "Arginine",
    "AGG": "Arginine",
    "GGU": "Glycine",
    "GGC": "Glycine",
    "GGA": "Glycine",
    "GGG": "Glycine"}

#Dictionnaire contenant l'équivalant en lettre de chaque codon (triplet bases azotée == lettre)
lettreAa = {
    "UUU": "F",
    "UUC": "F",
    "UUA": "L",
    "UUG": "L",
    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",
    "AUU": "I",
    "AUC": "I",
    "AUA": "I",
    "AUG": "M",
    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",
    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",
    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "UAU": "Y",
    "UAC": "Y",
    "UAA": "*",
    "UAG": "*",
    "CAU": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAU": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAU": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "UGU": "C",
    "UGC": "C",
    "UGA": "*",
    "UGG": "W",
    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGU": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G"
}

"""
######################
Sections des fonctions
######################
"""

#/!\ATTENTION/!\ les prints de chaque fonction servent de débogage pour suivre le processus pas-à-pas

#Fonction 1: trouve les debuts de ADN 
def trouveDebut(adn):
    positions_debut = []  # array contenant les position des d/buts
    for i in range(0,len(adn)):                                                 #On itère sur tout ADN
        if adn[i:i+3] == "TAC":                                                 #Si le sous-array i,i+3 == ATT, ATC ou ACT;
            positions_debut.append(i)                                           #On ajoute de sous-array a la liste de positions fin
    #print("debuts existants:")        
    #print(positions_debut)
    return positions_debut


#Fonction 2: trouve les fins de ADN
def trouveFin(adn):
    positions_fin = []                                                          #Initialisation array contenant les position des d/buts
    for i in range(len(adn)):                                                   #On itère sur tout ADN
        if adn[i:i+3] == "ATT" or adn[i:i+3] == "ATC" or adn[i:i+3] == "ACT":   #Si le sous-array i,i+3 == ATT, ATC ou ACT
            positions_fin.append(i)                                             #On ajoute de sous-array à la liste de positions fin
    #print("fins existants:")    
    #print(positions_fin)
    return positions_fin

#Fonction 3: trouve les genes valides
def trouveGene(a,b):
    gene_possible = []                                                          # initialisation array contenant les gènes possibles

    #cette partie sert à crééer tout les gènes possibles (donc des débuts suivis d'une fin
    #avec un multiple de 3 comme distance les séparant) 

    for debut in a:                                                             # On va faire un produit cartésien entre chaque début et chaque fin
        for fin in b:
            if (fin - debut) % 3 == 0 and (fin - debut) > 0:                    #si la fin suit le début et que la différence entre est un multiple de 3;
                gene_possible.append((debut, fin))                              #c'est une gène possible

    #cette partie sert à filter les genes; si une paire2 est entre une paire1, la paire la plus
    #imbriquee est gardee, si les paires partagent un debut ou une fin et que la premiere paire
    #est plus longue, alors la paire1 nest pas valide, sinon on la met dans la chaine de gene

    gene = []
    for paire1 in gene_possible:                                                            #premier comparateur, gene 1 de la liste
        paire1_valide = True                                                                #on le considère bon pour l'instant...
        for paire2 in gene_possible:                                                        #deuxième gène à comparer
            if paire1 != paire2 and ((paire2[0] < paire1[1] and paire2[1] > paire1[0]) or   #si la paire 2 est comprise dans la paire 1
                paire1[0] == paire2[0] or                                                   #si elles partagent le même début
                paire1[1] == paire2[1]):                                                    #si elles partagent la même fin alrs

                if (paire1[1] - paire1[0]) >= (paire2[1] - paire2[0]):                      #on vérifie si la paire 1 est plus longue
                    paire1_valide = False                                                   #si elle l'est alors paire 1 n'est plus valide par chevauchement et longueur
                    break
        if paire1_valide:
            gene.append(paire1)
    #print("genes valides:")
    #print(gene)
    return(gene)

#Fonction 4: transcrire le gene en chaine de bases azotees ("AAGGUUTACXXXACT" == "TACXXXACT")
def transcrireGene(adn,gene):
    chaine_proteine=""
    for paire in gene:
        fragment = adn[paire[0]:paire[1]]
        chaine_proteine = chaine_proteine + fragment
        fragment=""
    #print("chaine proteine:")
    #print(chaine_proteine)
    return chaine_proteine

#Fonction 5 : traduire le gene en arn (ATGC=UACG)
def traduireAdn(chaine_proteine):
    adnA = chaine_proteine.replace("A","U")     # A devient U,
    adnT = adnA.replace("T","A")                # T devient A
    adnG = adnT.replace("G","I")                # G devient une base de transition I
    adnC = adnG.replace("C","G")                # C devient G 
    arn = adnC.replace("I","C")                 # base de transition I devient G, arn transcrit totalement
    #print("traduction en ARN")
    #print(arn)
    return arn

#Fonction 6 : formatte les bases azotees en codons (a,u,g == "aug")
def arnCdons(arn):
    arn_codons = []                             #initialisation du array arn_de codons
    while arn:                                  #tant que arn n'est pas vide (donc == 0)
        arn_codons.append(arn[:3])              #on ajoute les 3 premieres lettres de arn comme 1 élément de arn_codons
        arn = arn[3:]                           #on retire les 3 premieres lettres de arn, on continue jusquà épuisement
    #print("format codons")
    #print(arn_codons)
    return arn_codons

#Fonction 7 : traduit les codons en acides aminees (aug = methionine)
def acidesAmines(arn_codons,codons_aa):         #CODONS
    acides_amine = []                           #initialisation du array acides_amine
    for i in arn_codons:                        #on itère chaque élément du array arn-codons
        if i in codons_aa:                      #si l'élément de l'index i est égal à un élément dans le dictionnaire
            acides_amine.append(codons_aa[i])   #on construit acides_amine avec l'acide aminé du dictionnaire
        else:
            acides_amine.append(i)              #dans le cas d'un adn dont le nbr de bases azotées est pas un multiple de 3 ou 
                                                #le codon n'existe pas dans le dictionnaire, on le laisse comme ca
    #print("format acides amines:")
    #print(acides_amine)
    return acides_amine

#Fonction 8 : formatte l array dacides amines en chaine de string
def arnString(acide_amine):
    arn_string = ""                             #initialisation du array arn_de codons
    for i in range(0, len(acide_amine)):
        if i != (len(acide_amine)-1):
            arn_string=arn_string + acide_amine[i] + "-"
        else:
            arn_string=arn_string + acide_amine[i]
    #print("format chaine de characteres acides amines:")
    #print(arn_string)
    return arn_string

#Fonction 9 : traduit les codons en leur equivalence en lettre (AUG = M)
def sigleLettre(arn_codons,lettreAa):           
    sigle = []                                  #initialisation du array acides_amine
    for i in arn_codons:                        #on itère chaque élément du array arn-codons
        if i in lettreAa:                       #si l'élément de l'index i est égal à un élément dans le dictionnaire
            sigle.append(lettreAa[i])           #on construit acides_amine avec l'acide aminé du dictionnaire
        else:
            sigle.append(i)                     #dans le cas d'un adn dont le nbr de bases azotées est pas un multiple de 3 ou 
                                                #le codon n'existe pas dans le dictionnaire, on le laisse comme ca
    #print(sigle)
    return sigle

#Fonction 10 : utilisation pour le brin complementaire, retourne l adn inverse et complementaire

def antisens(adn):
    adn_inv = adn[::-1]
    adnA = adn_inv.replace("A","I")     # A devient I, une base de transition
    adnT = adnA.replace("T","A")        # T devient A
    adnG = adnT.replace("I","T")        # I devient T
    adnC = adnG.replace("C","I")        # C devient I, une base de transition
    adnI = adnC.replace("G","C")        # G devient C
    adn = adnI.replace("I","G")         # I devient G, arn transcrit totalement
    #print("adn complementaire: ")
    #print(adn)
    return adn

"""
###########################################
Sections des fonctions pour tests unitaires
###########################################
"""

#Tests unitaires de fonction 10
def test_antisens():

    assert antisens("AGTC") == "GACT"
    assert antisens("AAAA") == "TTTT"
    assert antisens("") == ""
    assert antisens("CCCGGG") == "CCCGGG"
    print("Test antisens passé avec succès.")
#test_antisens()
    
#Tests unitaires de fonction 1
def test_trouveDebut():

    assert trouveDebut("TACTACAGT") == [0, 3]
    assert trouveDebut("AGTTAC") == [3]
    assert trouveDebut("TACTACTAC") == [0, 3, 6]
    assert trouveDebut("AAAAGG") == []
    print("Test trouveDebut passé avec succès.")
#test_trouveDebut()

#Test unitaires de fonction 2
def test_trouveFin():

    assert trouveFin("ATTCATCTACT") == [0, 4, 8]
    assert trouveFin("AGTATTGATCATCGACT") == [3, 7, 10, 14]
    assert trouveFin("CCCC") == []
    assert trouveFin("ACTATTATC") == [0, 3, 6]
    print("Test trouveFin passé avec succès.")
#test_trouveFin()
    
#Test unitaires de fonction 3
def test_trouveGene():

    assert trouveGene([0, 3], [6, 9]) == [(3, 6)]
    assert trouveGene([1, 4], [10, 15]) == [(4, 10)]
    assert trouveGene([], []) == []
    assert trouveGene([3, 15], [9, 21]) == [(3, 9), (15, 21)]
    print("Test trouveGene passé avec succès.")
#test_trouveGene()

#Tests unitaires de fonction 4
def test_transcrireGene():

    assert transcrireGene("AGTCAGTCA", [(0, 6)]) == "AGTCAG"
    assert transcrireGene("AGTCAG", [(0, 6)]) == "AGTCAG"
    assert transcrireGene("", []) == ""
    assert transcrireGene("AAAACCCCGGGG", [(4, 8)]) == "CCCC"
    print("Test transcrireGene passé avec succès.")
#test_transcrireGene()
    
#Tests unitaires de fonction 5
    
def test_traduireAdn():

    assert traduireAdn("ATG") == "UAC"
    assert traduireAdn("TTT") == "AAA"
    assert traduireAdn("") == ""
    assert traduireAdn("AGTC") == "UCAG"
    print("Test traduireAdn passé avec succès.")
#test_traduireAdn()

#Tests unitaires de fonction 6
def test_arnCdons():

    assert arnCdons("AUGUAC") == ["AUG", "UAC"]
    assert arnCdons("AUGUACAAA") == ["AUG", "UAC", "AAA"]
    assert arnCdons("") == []
    assert arnCdons("AUG") == ["AUG"]
    print("Test arnCdons passé avec succès.")
#test_arnCdons()

#Tests unitaires de fonction 7
def test_acidesAmines():

    codons_aa = {"AUG": "Méthionine", "UAC": "Tyrosine", "AAA": "Lysine"}
    assert acidesAmines(["AUG", "UAC"], codons_aa) == ["Méthionine", "Tyrosine"]
    assert acidesAmines(["AAA"], codons_aa) == ["Lysine"]
    assert acidesAmines([], codons_aa) == []
    assert acidesAmines(["AUG", "AAA", "UAC"], codons_aa) == ["Méthionine", "Lysine", "Tyrosine"]
    print("Test acidesAmines passé avec succès.")
#test_acidesAmines()
#Tests unitaires de fonction 8
def test_arnString():

    assert arnString(["Isoleucine", "Sérine"]) == "Isoleucine-Sérine"
    assert arnString(["Lysine"]) == "Lysine"
    assert arnString([]) == ""
    assert arnString(["Méthionine", "Lysine", "Tyrosine"]) == "Méthionine-Lysine-Tyrosine"
    print("Test arnString passé avec succès.")
#test_arnString

#Tests unitaires de fonction 9
def test_sigleLettre():

    lettreAa = {"AUG": "M", "UAC": "Y", "AAA": "K"}
    assert sigleLettre(["AUG", "UAC"], lettreAa) == ["M", "Y"]
    assert sigleLettre(["AAA"], lettreAa) == ["K"]
    assert sigleLettre([], lettreAa) == []
    assert sigleLettre(["AUG", "AAA", "UAC"], lettreAa) == ["M", "K", "Y"]
    print("Test sigleLettre passé avec succès.")

"""
Brin d'ADN utilisé pour le projet
"""
adn="TCGACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGATACCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATACCCAGCCAGCCAGCCAGCGACGGCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGAGTGCCAGCCAGCCAGCCAGCGAACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGAATGCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATATTCAGCCAGCCAGCCAGCGAACACTCTTCGACAGCCAGCGAAGCCAGCCAGCCGATATTCAGCCAGCCAGCCAGCGAACTCGACACTCTTCGACAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCATCCCAGCGATACCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAACTCGTCTGCGTTCGACAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGATTGCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGTATGCCAGCCAGCATCCCAGCGA"


#Pour afficher et lire en string les acides aminés du brin régulier
print("---------------------------------------")
print("Acides aminés du brin principal :")
print(arnString(acidesAmines(arnCdons(traduireAdn(transcrireGene(adn,trouveGene(trouveDebut(adn),trouveFin(adn))))),codons_aa)))

#Pour afficher et lire en string les acides aminés du brin régulier
print("---------------------------------------")
print("Acides aminés du brin complémentaire :")
print(arnString(acidesAmines(arnCdons(traduireAdn(transcrireGene(antisens(adn),trouveGene(trouveDebut(antisens(adn)),trouveFin(antisens(adn)))))),codons_aa)))

proteine_adn = sigleLettre(arnCdons(traduireAdn(transcrireGene(adn,trouveGene(trouveDebut(adn),trouveFin(adn))))),lettreAa)                              
proteine_complementaire = sigleLettre(arnCdons(traduireAdn(transcrireGene(antisens(adn),trouveGene(trouveDebut(antisens(adn)),trouveFin(antisens(adn)))))),lettreAa)

proteines=proteine_adn+proteine_complementaire      #array contenant le brin adn et adn complémentaire en sigle (format lettres)


"""
#####################################
Sections des fonctions pour le dessin
#####################################
"""

clear(800,600)

#Fonction de dessin de la lettre symbolique
def dessinLettre(acideA):
    pu();fd(10);lt(90);fd(10);rt(90);pd()       #positionnement de la lettre, on l'aligne au centre du carré quon dessinera
    if acideA == 'M':                           #on vérifie tout les cas de correspondance de lettre existants (pris du dictionnaire)
        turtle.write('M')
    elif acideA == 'F':
        turtle.write('F')
    elif acideA == 'L':
        turtle.write('L')
    elif acideA == 'I':
        turtle.write('I')
    elif acideA == 'V':
        turtle.write('V')
    elif acideA == 'S':
        turtle.write('S')
    elif acideA == 'P':
        turtle.write('P')
    elif acideA == 'T':
        turtle.write('t')
    elif acideA == 'A':
        turtle.write('A')
    elif acideA == 'Y':
        turtle.write('Y')
    elif acideA == 'H':
        turtle.write('H')
    elif acideA == 'Q':
        turtle.write('Q')
    elif acideA == 'N':
        turtle.write('N')
    elif acideA == 'K':
        turtle.write('K')
    elif acideA == 'D':
        turtle.write('D')
    elif acideA == 'E':
        turtle.write('E')
    elif acideA == 'C':
        turtle.write('C')
    elif acideA == 'W':
        turtle.write('W')
    elif acideA == 'R':
        turtle.write('R')
    elif acideA == 'G':
        turtle.write('G')

#Fonction de dessin de carrés

def carre(proteines):
    positionX=-170                                  #on comence le dessin en -170,300
    positionY=300
    goto(positionX,positionY)
    compte=0

    for acideA in proteines:                        #nbr de carres à dessiner
        if compte%15 == 0 or acideA=='M':           #si le compte de lettres pour un gène depasse 15 ou on trouve un nouveau début
            positionX=-170                          #on retourne en x à -170, mais on descend d'un carré
            positionY-=20
            pu()
            goto(positionX,positionY)
            pd()
            compte=0
        else:
            positionX+=20                           #sinon, on dessine le prochain carré à droite
            pu()
            goto(positionX,positionY)
            pd()
        for j in range(0,4):                        #fct qui dessine les carres
            fd(20)                                  #loop qui dessine les 4 cotés du carré
            lt(90)
        dessinLettre(acideA)   
        compte+=1

carre(proteines) #appel de fonction pour dessiner la protéines (tout les gènes présents dans la séquence adn)