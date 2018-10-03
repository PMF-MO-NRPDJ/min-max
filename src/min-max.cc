#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <limits>

#include "dune/common/parallel/mpihelper.hh"
#include <dune/common/exceptions.hh> 

#include <dune/grid/uggrid.hh>  
#include <dune/grid/common/gridinfo.hh> 
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

// Izračunaj kut između p1-p2 i p3-p2
template <typename Point>
double kut(Point const & p1, Point const & p2, Point const & p3)
{
   // računanje kuta (VAŠ KOD) ide ovdje    
   // Point će biti Dune::FieldVector<double, dim>. Norma se računa 
   // pomoću funkcije članice klase two_norm(), a skalarni produkt 
   // pomoću funkcije članice klase dot(). Pogledati dokumentaciju klase.
    return kut;
}

int main(int argc, char** argv)
{
    const int dim = 2;
    typedef Dune::UGGrid<dim> GridType;
    GridType * p_grid =nullptr;

    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

    if(argc < 2){
      std::cout << "Need msh file name!" << std::endl;
      std::exit(-1);
    } 

    // Kreiramo dvodimenzionalnu mrežu čitanjem gmsh datoteke
    p_grid = Dune::GmshReader<GridType>::read(argv[1]);

    // Ako je mreža paralelna rasporedi je ravnomjerno po procesorima
    p_grid->loadBalance();
    if(helper.rank() == 0) Dune::gridinfo(*p_grid);
   
    // Uzmi GridView 
    typedef typename GridType::LeafGridView LeafGridView;
    LeafGridView gridView = p_grid->leafGridView(); 
    
    double min = std::numeric_limits<double>::max();  // najveći double
    double max = std::numeric_limits<double>::lowest(); // najveći negativni double
    int count = 0;
    for(auto const & element : elements(gridView))
    {

/*     VAŠ KOD dolazi ovdje.
 *     alpha, beta i gamma su kutevi trokuta na kojem se nalazimo
 *     loc_min i loc_max su minimalni i maksimalni kut trokuta
       double alpha = ...
       double beta  = ...
       double gamma = M_PI - alpha -beta;
       double loc_min =  ...
       double loc_max = ...
*/
       std::cout <<"Element " << count << " min kut = " <<  180*loc_min/M_PI 
                 << ", max kut = " << 180*loc_max/M_PI << "\n";

       min = std::min(loc_min, min);
       max = std::max(loc_max, max);
       count++;
    } 

    std::cout << "Broj (leaf) elemenata = " << count  << std::endl;
    std::cout << "Minimalni kut  = " << 180*min/M_PI 
              << ", maksimalni kut = " << 180*max/M_PI  << std::endl;


    // Ispis mreže u VTK formatu (u datoteku poluvijenac.vtu)
    Dune::VTKWriter<LeafGridView> vtkwriter(gridView);
    vtkwriter.write("poluvijenac");

    delete p_grid;
 
    return 0;
}
