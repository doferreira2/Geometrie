#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <iostream>
#include <fstream>
#include <time.h>



typedef struct s_color
{
	float R;
	float V;
	float B;
} color;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

typedef Polyhedron::Facet_const_iterator Facet_iterator;
typedef Polyhedron::Vertex_const_iterator Vertex_iterator;
typedef Polyhedron::Halfedge_const_iterator Halfedge_iterator;
typedef Polyhedron::Halfedge_around_facet_const_circulator Halfedge_facet_circulator;
typedef Polyhedron::Facet_handle Facet_handle;

typedef std::map<Polyhedron::Facet_const_handle, double> Facet_double_map;
typedef std::map<Polyhedron::Facet_const_handle, int> Facet_int_map;
typedef std::map<Polyhedron::Facet_const_handle, bool> Facet_bool_map;

// variable global
int N_Class;
int Nb_Color;

void initColor(color *colorPalette)
{
	for (int i = 0; i < Nb_Color; i++)
	{
		colorPalette[i].R = (rand() * i % 1000) / 1000.;
		colorPalette[i].V = (rand() * (50 - i) % 1000) / 1000.;
		colorPalette[i].B = (rand() * (100 - i) % 1000) / 1000.;
	}
}

void computBorne(Facet_double_map &facetMap, double &min, double &max)
{
	max = facetMap.begin()->second;
	min = facetMap.begin()->second;

	for (const auto &elem : facetMap)
	{
		if (elem.second > max)
		{
			max = elem.second;
		}
		if (elem.second < min)
		{
			min = elem.second;
		}
	}
}

/// @brief map all the values from [min, max] to [0, 1]
/// @param facetMap non-const reference to the map (it is an in/out parameter)
void normalizeMap(Facet_double_map &facetMap)
{
	double maxValue, minValue;

	computBorne(facetMap, minValue, maxValue);

	for (auto &elem : facetMap)
	{
		elem.second -= minValue;
		elem.second /= (maxValue - minValue);
	}
}

/**
 * @brief Permet de calculer un tableau de avec la classe de chaque face
 *
 * @param facetMap
 * @return Facet_int_map
 */
Facet_int_map simpleSegmentation(Facet_double_map &facetMap)
{
	double avgValue = 0;
	int nb = 0, em = 0;
	Facet_iterator i = 0;
	Facet_int_map ClassInt;

	// Initialize the min and max perimeters
	double maxValue, minValue;

	computBorne(facetMap, minValue, maxValue);

	float intervalSize = (maxValue - minValue) / N_Class;

	for (const auto &elem : facetMap)
	{
		int intervalIndex = (elem.second - minValue) / intervalSize;

		if (intervalIndex < 0)
		{
			intervalIndex = 0;
		}
		if (intervalIndex >= N_Class)
		{
			intervalIndex = N_Class - 1;
		}

		ClassInt[elem.first] = intervalIndex;
	}
	Nb_Color = N_Class;
	return ClassInt;
}

/// @brief Generate in a .off file a colored mesh according to a value map (green to red shades)
/// @param mesh the input mesh
/// @param facetMap map of values between 0 and 1 (see "normalize()") for each facet of mesh
/// @param filePath path to the colored .off file to be generated
void writeOFFfromValueMap(const Polyhedron &mesh, const Facet_double_map &facetMap, Facet_int_map &facetIntM, std::string filePath)
{
	std::ofstream in_myfile;
	in_myfile.open(filePath);

	CGAL::set_ascii_mode(in_myfile);

	in_myfile << "COFF" << std::endl // "COFF" makes the file support color informations
			  << mesh.size_of_vertices() << ' '
			  << mesh.size_of_facets() << " 0" << std::endl;
	// nb of vertices, faces and edges (the latter is optional, thus 0)

	std::copy(mesh.points_begin(), mesh.points_end(),
			  std::ostream_iterator<Kernel::Point_3>(in_myfile, "\n"));

	color colorPalette[Nb_Color];
	initColor(colorPalette);

	for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
	{
		Halfedge_facet_circulator j = i->facet_begin();

		CGAL_assertion(CGAL::circulator_size(j) >= 3);

		in_myfile << CGAL::circulator_size(j) << ' ';
		do
		{
			in_myfile << ' ' << std::distance(mesh.vertices_begin(), j->vertex());

		} while (++j != i->facet_begin());

		in_myfile << std::setprecision(5) << std::fixed; // set the format of floats to X.XXXXX

		auto redValue = colorPalette[facetIntM.at(i)].R;   // low values will be closer to red
		auto greenValue = colorPalette[facetIntM.at(i)].V; // high values will be closer to green
		auto blueValue = colorPalette[facetIntM.at(i)].B;

		in_myfile << " " << redValue << " " << greenValue << " " << blueValue;

		in_myfile << std::endl;
	}

	in_myfile.close();

	std::cout << "Le résultat a été exporté dans " << filePath << " !" << std::endl;
}

/**
 * @brief Cree un tableau avec les perimaitre de chaque face
 *
 * @param mesh
 * @return Facet_double_map
 */
Facet_double_map computePerimMap(const Polyhedron &mesh)
{
	Facet_double_map out;

	for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
	{
		double current_perimeter = 0.;
		Halfedge_facet_circulator j = i->facet_begin();
		do
		{
			current_perimeter += std::sqrt(CGAL::squared_distance(j->vertex()->point(), j->opposite()->vertex()->point()));
		} while (++j != i->facet_begin());

		std::cout << "perim(" << std::distance(mesh.facets_begin(), i) << ")=" << current_perimeter << std::endl;

		out[i] = current_perimeter;
	}

	return out;
}
/**
 * @brief Cree un tableau avec les aires de chaque faces
 *
 * @param mesh
 * @return Facet_double_map
 */
Facet_double_map computeAreaMap(const Polyhedron &mesh)
{
	Facet_double_map out;

	for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
	{

		Halfedge_facet_circulator j = i->facet_begin();

		Polyhedron::Vertex_const_handle firstVertex = j->vertex();

		double current_area = 0;
		// a facet is not necessarily a triangle, so we decompose one facet into multiple triangles,
		// and sum up all their areas. Only works for convex faces.
		// (illustration: http://mathbitsnotebook.com/JuniorMath/Polygons/polygons3g.jpg)
		do
		{
			current_area += CGAL::squared_area(firstVertex->point(), j->vertex()->point(), j->opposite()->vertex()->point());
		} while (++j != i->facet_begin());

		std::cout << "area(" << std::distance(mesh.facets_begin(), i) << ")=" << current_area << std::endl;

		out[i] = current_area;
	}

	return out;
}

/**
 * @brief Fonction de parcour des voisin d'un face recursivement
 *
 * @param mesh
 * @param facet
 * @param visit
 * @param segOut
 * @param segmentation
 * @param Nclass
 */
void visit_facet(Polyhedron &mesh, Facet_iterator facet, Facet_bool_map &visitMap, Facet_int_map &segOut, Facet_int_map &segIni, int Nclass)
{
	// Mark the facet as visited
	visitMap[facet] = true;
	segOut[facet] = Nclass;

	Halfedge_facet_circulator he = facet->facet_begin();

	do
	{
		Facet_iterator face_voisin = he->opposite()->facet();

		if (face_voisin != NULL)
		{
			if (!visitMap[face_voisin] && segIni[face_voisin] == segIni[facet])
			{
				visit_facet(mesh, face_voisin, visitMap, segOut, segIni, Nclass);
			}
		}

	} while (++he != facet->facet_begin());
}
/**
 * @brief Permet de faire un changement de class en fonction des face adajcente
 *
 * @param mesh
 * @param segmentation
 * @return Facet_int_map
 *
 */
Facet_int_map segmentationParCC(Polyhedron &mesh, Facet_int_map &segmentation)
{
	Facet_int_map segOut = segmentation;
	Facet_bool_map is_visit;

	int nClass = 0;

	//  tableau de visit
	for (Facet_iterator f = mesh.facets_begin(); f != mesh.facets_end(); ++f)
	{
		is_visit[f] = false;
	}

	// parcour du mesh
	for (Facet_iterator f = mesh.facets_begin(); f != mesh.facets_end(); ++f)
	{
		if (!is_visit[f])
		{
			visit_facet(mesh, f, is_visit, segOut, segmentation, nClass);
			nClass++;
		}
	}

	Nb_Color = nClass;
	return segOut;
}

void RecVisit(Polyhedron &mesh, Facet_iterator facet, Facet_bool_map &visit, Facet_int_map &seg, Facet_double_map &data, int nbClass, float seuil)
{
	visit[facet] = true;
	seg[facet] = nbClass;
	Halfedge_facet_circulator he = facet->facet_begin();

	do
	{

		Facet_iterator neighbor_facet = he->opposite()->facet();

		// Check if neighbor_facet is in the same segmentation class and has not been visited yet
		if (neighbor_facet != NULL)
		{
			if (!visit[neighbor_facet] && (fabs(data[neighbor_facet] - data[facet]) <= seuil))
			{
				RecVisit(mesh, neighbor_facet, visit, seg, data, nbClass, seuil);
			}
		}

	} while (++he != facet->facet_begin());
}

Facet_int_map complexSegmentation(Polyhedron &mesh, Facet_double_map &facetMap)
{
	Facet_int_map clxSeg;
	Facet_bool_map is_visit;
	int nClass = 0;

	double maxValue = facetMap.begin()->second;
	double minValue = facetMap.begin()->second;

	computBorne(facetMap, minValue, maxValue);

	float seuil = (maxValue - minValue) / N_Class;

	for (Facet_iterator f = mesh.facets_begin(); f != mesh.facets_end(); ++f)
	{
		is_visit[f] = false;
		clxSeg[f] = 0;
	}

	for (Facet_iterator f = mesh.facets_begin(); f != mesh.facets_end(); ++f)
	{
		if (!is_visit[f])
		{
			RecVisit(mesh, f, is_visit, clxSeg, facetMap, nClass, seuil);
			nClass++;
		}
	}
	Nb_Color = nClass;

	return clxSeg;
}

// ---------------------------------------------------
int main(int argc, char *argv[])
{

	if (argc < 2)
	{
		std::cerr << "Il manque un paramètre au programme. Veuillez lui donner en entrée un nom de fichier au format off." << std::endl;
		return 1;
	}

	auto mode = (argc >= 3 ? argv[2] : "A");
	N_Class = (argc >= 4 ? atoi(argv[3]) : 5);

	Polyhedron mesh;

	std::ifstream input(argv[1]);

	if (!input || !(input >> mesh) || mesh.is_empty())
	{
		std::cerr << "Le fichier donné n'est pas un fichier off valide." << std::endl;
		return 1;
	}
	srand(time(NULL));

	Facet_double_map mapPerim;

	if (strcmp(mode, "P") == 0)
	{
		mapPerim = computePerimMap(mesh);
	}
	else
	{
		mapPerim = computeAreaMap(mesh);
	}

	normalizeMap(mapPerim);

	auto PerimInt = simpleSegmentation(mapPerim);
	std::cout << "Nb Class avec segmentation simple : " << Nb_Color << std::endl;
	writeOFFfromValueMap(mesh, mapPerim, PerimInt, "result.off");

	auto Em = segmentationParCC(mesh, PerimInt);
	std::cout << "Nb Class avec segmentation CC : " << Nb_Color << std::endl;
	writeOFFfromValueMap(mesh, mapPerim, Em, "resultCC.off");

	auto cpS = complexSegmentation(mesh, mapPerim);
	std::cout << "Nb Class avec segmentation complex : " << Nb_Color << std::endl;
	writeOFFfromValueMap(mesh, mapPerim, cpS, "resultcompl.off");


	std::ofstream file;
	file.open("test.txt");
	file << "simple, CC, complex" << std::endl;
	for (Facet_iterator f = mesh.facets_begin(); f != mesh.facets_end(); ++f)
	{
		file << PerimInt[f] << "," << Em[f] << "," << cpS[f] << std::endl;
	}
	file.close();

	return 0;
}
