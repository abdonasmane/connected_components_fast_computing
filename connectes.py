#!/usr/bin/env python3
"""
compute sizes of all connected components.
sort and display.
"""
from sys import argv


k = 2 #dimension des points, donc on construira des 2DTree

class KDTree:
    """
    KDtree, où on va stocker les points
    """
    def __init__(self, points):
        def construction_KDTree(points, profondeur):
            """
            On construit recursivement la KDTree
            """
            if not points:
                return None
            # on calcule la dimension qui sera considérée pour diviser le plan à la profondeur profondeur
            axe = profondeur % k

            # on cherche le point optimale pour diviser le plan : la median
            points.sort(key=lambda p: p.coordinates[axe])

            indice_mediane = len(points) // 2
            # On construit les sous arbres gauche et droite récursivement
            return (
                points[indice_mediane],
                construction_KDTree(points[:indice_mediane], profondeur+1),
                construction_KDTree(points[indice_mediane+1:], profondeur+1)
            )

        self.root = construction_KDTree(points, profondeur=0)


    def recherche_voisins(self, point, rayon):
        """
        cherche les point qui se trouvent dans le disque de centre "point" et de rayon "rayon"
        """
        ray_2 = rayon**2 # pour éviter de calculer le carrée à chaque fois
        def recherche_sous_arbre(noeud, centre, ray_2, profondeur=0):
            if not noeud:
                return []

            voisins_trouvés = []
            if (centre.coordinates[0] - noeud[0].coordinates[0])**2 + (centre.coordinates[1] - noeud[0].coordinates[1])**2 <= ray_2:
                voisins_trouvés.append(noeud[0])

            axe = profondeur % k
            if noeud[1] and centre.coordinates[axe] - rayon <= noeud[0].coordinates[axe]:
                voisins_trouvés+=(recherche_sous_arbre(noeud[1], centre, ray_2, profondeur+1))

            if noeud[2] and centre.coordinates[axe] + rayon >= noeud[0].coordinates[axe]:
                voisins_trouvés+=(recherche_sous_arbre(noeud[2], centre, ray_2, profondeur+1))

            return voisins_trouvés

        return recherche_sous_arbre(self.root, point, ray_2)


class mon_Point:
    """ 
    represente un point et sa collection
    """
    def __init__(self, coordonnees):
        self.coordinates = coordonnees
        self.collection = Collection_de_point(self)


class Collection_de_point:
    """ 
    represente une collection de points
    """
    def __init__(self, racine):
        self.collection = [racine]

    def add_collection(self, nouv_collection):
        """
        fusionne 2 collections de pts
        """
        if self == nouv_collection:
            return
        for pt in nouv_collection.collection:
            self.collection.append(pt)
            pt.collection = self

    def __len__(self):
        return len(self.collection)


def load_instance(filename):
    """
    loads .pts file.
    returns distance limit and points.
    """
    with open(filename, "r") as instance_file:
        lines = iter(instance_file)
        distance = float(next(lines))
        points = [mon_Point([float(f) for f in l.split(",")]) for l in lines]
    return distance, points


def Partition_avec_matrice(points, distance):
    """
    partitionne le plan sous forme de carrés et calcule les composantes connexes
    """
    cote, dis2 = distance/(2**0.5), distance**2
    nb_carres_par_ligne = int(1/cote)+1
    # initialisation de la matrice
    carres = [[[] for _ in range(nb_carres_par_ligne)] for _ in range(nb_carres_par_ligne)]
    # insertion des points dans la mattrice
    for p_t in points:
        carres[int(p_t.coordinates[0]//cote)][int(p_t.coordinates[1]//cote)].append(p_t)
    # traitement des points de chaque carré
    for carres_x in carres:
        for carres_y in carres_x:
            for p_t in carres_y[1:]:
                carres_y[0].collection.add_collection(p_t.collection)

    mvt = [(0, 1), (1, 1), (1, 0), (1, -1), (0, 2), (1, 2), (2, 1), (2, 0), (2, -1), (1, -2)]
    # recherche parmis les carrés voisins lesquelles entreront dans la composante connexe
    for i in range(len(carres)):
        for j in range(len(carres)):
            if len(carres[i][j]) == 0:
                continue

            for mv in mvt:
                k = i + mv[0]
                l = j + mv[1]
                if k < 0 or k >= len(carres):
                    continue
                if l < 0 or l >= len(carres):
                    continue

                if conn_exists(carres[i][j], carres[k][l], dis2):
                    carres[i][j][0].collection.add_collection(carres[k][l][0].collection)


def Partition_avec_dictionnaire(pts, distance):
    """
    partitionne le plan sous forme de carrés et calcule les composantes connexes
    """
    cote, dis2 = distance/(2**0.5), distance**2
    nb_carres = int(1/cote)+1
    # initialisation du dictionnaire
    cases_points = {}
    for p_t in pts:
        key = int((p_t.coordinates[0]//cote)+nb_carres*(p_t.coordinates[1]//cote))
        try:
            cases_points[key].append(p_t)
        except KeyError:
            cases_points[key] = [p_t]
    
    # traitement des points de chaque carré
    for carre1 in cases_points:
        head = cases_points[carre1][0]
        for point in cases_points[carre1][1:]:
            head.collection.add_collection(point.collection)

    mvt = [(0, 1), (1, 1), (1, 0), (1, -1), (0, 2), (1, 2), (2, 1), (2, 0), (2, -1), (1, -2)]
    # recherche parmis les carrés voisins lesquelles entreront dans la composante connexe
    for ind_case in cases_points:
        j, i = ind_case//nb_carres, ind_case%nb_carres
        for mv in mvt:
                ind_candidat = i + mv[0] + nb_carres * (j+mv[1])
                if ind_candidat in cases_points and conn_exists(cases_points[ind_case], cases_points[ind_candidat], dis2):
                    cases_points[ind_case][0].collection.add_collection(cases_points[ind_candidat][0].collection)


def kdtree(points, dist_limit):
    """
    affichage des tailles triees de chaque composante
    """
    # on initialise notre KDTree 
    tree = KDTree(points)
    # ici on marque les points visités
    visited = set()
    # ici on marque les points non visités
    unvisited = set(points)
    # c'est là où on stockera les tailles des composantes connexes
    comp_sizes = []
    while unvisited:
        if len(visited) == len(points):
            break
        curr_point = unvisited.pop()
        visited.add(curr_point)
        comp_points = set([curr_point])
        queue = [curr_point]
        while queue:
            curr_point = queue.pop()
            voisins = KDTree.recherche_voisins(tree, curr_point, dist_limit)
            for voisin in set(voisins) - visited:
                queue.append(voisin)
                visited.add(voisin)
                comp_points.add(voisin)
            if len(visited) == len(points):
                break
        # on supprime tous les points qui appartiennent à la composante connexe trouvée du unvisited
        unvisited.difference_update(comp_points)
        comp_sizes.append(len(comp_points))
    comp_sizes.sort(reverse=True)
    print(comp_sizes)


def conn_exists(carre1, carre2, dis2):
    """ vérifie si 2 carres appartiennent à la même composante """
    for pt in carre1:
        for pt2 in carre2:
            if (pt.coordinates[0]-pt2.coordinates[0])**2 + (pt.coordinates[1]-pt2.coordinates[1])**2 <= dis2:
                return True

    return False


def print_components_sizes(funct, distance, points):
    """
    affichage des tailles triees de chaque composante
    """
    if funct == kdtree:
        funct(points, distance)
    else:
        funct(points, distance)
        graphs = list(set([pt.collection for pt in points]))
        res = [len(composante) for composante in graphs]
        print(sorted(res, reverse=1))


def main():
    for instance in argv[1:]:
        distance, points = load_instance(instance)
        if 0.009 <= distance <= 1:
            print_components_sizes(Partition_avec_matrice, distance, points)
        elif 0.005 < distance < 0.009:
            if len(points) <= 200_000:
                print_components_sizes(kdtree, distance, points)
            else:
                print_components_sizes(Partition_avec_matrice, distance, points)
        elif 0.002 < distance <= 0.005:
            if len(points) <= 160_000:
                print_components_sizes(Partition_avec_matrice, distance, points)
            else:
                print_components_sizes(kdtree, distance, points)
        elif 0.001 < distance <= 0.002:
            if len(points) <= 1_000_000:
                print_components_sizes(Partition_avec_dictionnaire, distance, points)
            elif 1_000_000 < len(points) <= 1_350_000:
                print_components_sizes(Partition_avec_matrice, distance, points)
            else:
                print_components_sizes(kdtree, distance, points)
        else:
                print_components_sizes(Partition_avec_dictionnaire, distance, points)


main()
