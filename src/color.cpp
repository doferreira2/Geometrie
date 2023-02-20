#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <iostream>
#include <fstream>
#include <vector>

#define N_Class 5

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

// Define the color palette
color colorPalette[N_Class] = { { 1.0, 0.0, 0.0 }, // red
                               { 1.0, 0.5, 0.0 }, // orange
                               { 1.0, 1.0, 0.0 }, // yellow
                               { 0.0, 1.0, 0.0 }, // green
                               { 0.0, 0.0, 1.0 }  // blue
                             };

/// @brief map all the values from [min, max] to [0, 1]
/// @param facetMap non-const reference to the map (it is an in/out parameter)
void normalizeMap(Facet_double_map &facetMap)
{
	double maxValue = facetMap.begin()->second;
	double minValue = facetMap.begin()->second;

	// look for min and max value in the map
	for (const auto &elem : facetMap)
	{
		if (elem.second > maxValue)
		{
			maxValue = elem.second;
		}
		if (elem.second < minValue)
		{
			minValue = elem.second;
		}
	}

	for (auto &elem : facetMap)
	{
		elem.second -= minValue;
		elem.second /= (maxValue - minValue);
	}
}

Facet_int_map simpleThershold(Facet_double_map &facetMap)
{
	double avgValue = 0;
	int nb = 0, em = 0;
	Facet_iterator i = 0;
	Facet_int_map PerimInt;

	// Initialize the min and max perimeters
	double maxValue = facetMap.begin()->second;
	double minValue = facetMap.begin()->second;

	for (const auto &elem : facetMap)
	{
		if (elem.second > maxValue)
		{
			maxValue = elem.second;
		}
		if (elem.second < minValue)
		{
			minValue = elem.second;
		}
	}
	double perimeterThreshold = (minValue + maxValue) / 2.0;
	float intervalSize = (maxValue - minValue) / N_Class;


	for (const auto &elem : facetMap)
	{
		int intervalIndex = (elem.second - minValue) / intervalSize;

		if (intervalIndex < 0) {
			intervalIndex = 0;
		}
		if (intervalIndex >= N_Class) {
			intervalIndex = N_Class - 1;
		}

		PerimInt[elem.first] = intervalIndex;
	}
	return PerimInt;
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

		auto redValue = colorPalette[facetIntM.at(i)].R; // low values will be closer to red
		auto greenValue = colorPalette[facetIntM.at(i)].V; 				   // high values will be closer to green
		auto blueValue = colorPalette[facetIntM.at(i)].B;

		in_myfile << " " << redValue << " " << greenValue << " " << blueValue;

		in_myfile << std::endl;
	}

	in_myfile.close();

	std::cout << "Le résultat a été exporté dans " << filePath << " !" << std::endl;
}

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

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		std::cerr << "Il manque un paramètre au programme. Veuillez lui donner en entrée un nom de fichier au format off." << std::endl;
		return 1;
	}

	Polyhedron mesh;

	std::ifstream input(argv[1]);

	if (!input || !(input >> mesh) || mesh.is_empty())
	{
		std::cerr << "Le fichier donné n'est pas un fichier off valide." << std::endl;
		return 1;
	}

	auto mapPerim = computePerimMap(mesh);

	normalizeMap(mapPerim);

	std::cout << "test1" << std::endl;
	auto PerimInt = simpleThershold(mapPerim);
	std::cout << "test2" << std::endl;

	writeOFFfromValueMap(mesh, mapPerim, PerimInt, argc >= 3 ? argv[2] : "result.off");

	return 0;
}
