/*This is a template file for use with 3D finite elements (vector field).
  The portions of the code you need to fill in are marked with the comment "//EDIT".

  Do not change the name of any existing functions, but feel free
  to create additional functions, variables, and constants.
  It uses the deal.II FEM library.*/

//Include files
//Data structures and solvers
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
//Mesh related classes
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
//Finite element implementation classes
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
//Standard C++ libraries
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace dealii;

//Define the order of the basis functions (Lagrange polynomials)
//and the order of the quadrature rule globally
const unsigned int order = 1;
const unsigned int quadRule = 2;

template <int dim>
class FEM
{
 public:
  //Class functions
  FEM (); // Class constructor 
  ~FEM(); //Class destructor

  //Function to calculate components of the elasticity tensor
  double C(unsigned int i,unsigned int j,unsigned int k,unsigned int l);

  //Solution steps
  void generate_mesh(std::vector<unsigned int> numberOfElements);
  void define_boundary_conds();
  void setup_system();
  void assemble_system();
  void solve();
  void output_results();

  //Class objects
  Triangulation<dim> triangulation; //mesh
  FESystem<dim>      fe;            //FE element
  DoFHandler<dim>    dof_handler;   //Connectivity matrices

  //NEW - deal.II quadrature
  QGauss<dim>   quadrature_formula;      //Quadrature
  QGauss<dim-1> face_quadrature_formula; //Face Quadrature

  //Data structures
  SparsityPattern      sparsity_pattern; //Sparse matrix pattern
  SparseMatrix<double> K;                //Global stiffness matrix - Sparse matrix - used in the solver
  Vector<double>       D, F;             //Global vectors - Solution vector (D) and Global force vector (F)

  Table<2,double>	        dofLocation;	 //Table of the coordinates of dofs by global dof number
  std::map<unsigned int,double> boundary_values; //Map of dirichlet boundary conditions

  //solution name array
  std::vector<std::string> nodal_solution_names;
  std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation;
};

// Class constructor for a vector field
template <int dim>
FEM<dim>::FEM ()
:
fe (FE_Q<dim>(order), dim),
  dof_handler (triangulation),
  quadrature_formula(quadRule),
  face_quadrature_formula(quadRule)
{	
	
  //Nodal Solution names - this is for writing the output file
  for (unsigned int i=0; i<dim; ++i){
    nodal_solution_names.push_back("u");
    nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  }
}

//Class destructor
template <int dim>
FEM<dim>::~FEM (){dof_handler.clear ();}

//Function to calculate the components of the 4th order elasticity tensor
template <int dim>
double FEM<dim>::C(unsigned int i,unsigned int j,unsigned int k,unsigned int l){

  //Define the material parameters of Young's modulus and Poisson's ratio
  double E=2.0e11 ,  //EDIT
    nu=0.3 ; //EDIT
  double lambda=(E*nu)/((1.+nu)*(1.-2.*nu)),
    mu=E/(2.*(1.+nu));

  return lambda*(i==j)*(k==l) + mu*((i==k)*(j==l) + (i==l)*(j==k));

}

//Define the problem domain and generate the mesh
template <int dim>
void FEM<dim>::generate_mesh(std::vector<unsigned int> numberOfElements){

  //Define the limits of your domain
  double x_min = 0.0, //EDIT
    x_max = 1.0, //EDIT
    y_min = 0.0, //EDIT
    y_max = 1.0, //EDIT
    z_min = 0.0, //EDIT
    z_max = 1.0; //EDIT

  Point<dim,double> min(x_min,y_min,z_min),
    max(x_max,y_max,z_max);
  GridGenerator::subdivided_hyper_rectangle (triangulation, numberOfElements, min, max);
}

//Specify the Dirichlet boundary conditions
template <int dim>
void FEM<dim>::define_boundary_conds(){

  //EDIT - Define the Dirichlet boundary conditions.
	
  /*Note: this will be very similiar to the define_boundary_conds function
    in the HW2 template. You will loop over all degrees of freedom and use "dofLocation"
    to check if the dof is on the boundary with a Dirichlet condition.

    (Aside: Since we	now have more than 1 dof per node, it is possible to apply a different Dirichlet
    condition on each of the 3 dofs on the same node. Note that this is NOT the case for the current 
    assignment. But if you wanted to do it, you would need to check
    the nodal dof (0, 1, or 2) as well as the location. For example, if you wanted to fix displacements
    only in the x-direction at x=0, you would have an condition such as this:
    (dofLocation[globalDOF][0] == 0 && nodalDOF == 0)
    Note that you can get the nodal dof from the global dof using the "modulo" operator,
    i.e nodalDOF = globalDOF%dim. Here "%" gives you the remainder when dividing globalDOF by dim.)

    Add the global dof number and the specified value (displacement in this problem)
    to the boundary values map, something like this:

    boundary_values[globalDOFIndex] = dirichletDisplacementValue

    Note that "dofLocation" is now a Table of degree of freedom locations, not just node locations,
    so, for example, dofs 0, 1, and 2 correspond to the same node, so that have the same coordinates.
    The row index is the global dof number; the column index refers to the x, y, or z component (0, 1, or 2 for 3D).
    e.g. dofLocation[7][2] is the z-coordinate of global dof 7*/

  const unsigned int totalDOFs = dof_handler.n_dofs(); //Total number of degrees of freedom

for(unsigned int globalDof=0; globalDof<totalDOFs; globalDof++)
{
 if(dofLocation[globalDof][2]==0)
{
 boundary_values[globalDof] = 0.0;
}
}

}

//Setup data structures (sparse matrix, vectors)
template <int dim>
void FEM<dim>::setup_system(){

  //Let deal.II organize degrees of freedom
  dof_handler.distribute_dofs (fe);

  //Get a vector of global degree-of-freedom x-coordinates
  MappingQ1<dim,dim> mapping;
  std::vector< Point<dim,double> > dof_coords(dof_handler.n_dofs());
  dofLocation.reinit(dof_handler.n_dofs(),dim);
  DoFTools::map_dofs_to_support_points<dim,dim>(mapping,dof_handler,dof_coords);
  for(unsigned int i=0; i<dof_coords.size(); i++){
    for(unsigned int j=0; j<dim; j++){
      dofLocation[i][j] = dof_coords[i][j];
    }
  }

  //Specify boundary condtions (call the function)
  define_boundary_conds();

  //Define the size of the global matrices and vectors
  sparsity_pattern.reinit (dof_handler.n_dofs(), 
			   dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();
  K.reinit (sparsity_pattern);
  D.reinit(dof_handler.n_dofs());
  F.reinit(dof_handler.n_dofs());

  //Just some notes...
  std::cout << "   Number of active elems:       " << triangulation.n_active_cells() << std::endl;
  std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;   
}

//Form elmental vectors and matrices and assemble to the global vector (F) and matrix (K)
template <int dim>
void FEM<dim>::assemble_system(){

  /*NEW - deal.II basis functions, etc. The third input values (after fe and quadrature_formula) 
    specify what information we	want to be updated. For fe_values, we need the basis function values,
    basis function gradients,	and det(Jacobian) times the quadrature weights (JxW). For fe_face_values,
    we need the basis function values, the value of x at the quadrature points, and JxW.*/

  //For volume integration/quadrature points
  FEValues<dim> fe_values (fe,
			   quadrature_formula, 
			   update_values | 
			   update_gradients | 
			   update_JxW_values);

  //For surface integration/quadrature points
  FEFaceValues<dim> fe_face_values (fe,
				    face_quadrature_formula, 
				    update_values | 
				    update_quadrature_points | 
				    update_JxW_values);

  K=0; F=0;

  const unsigned int dofs_per_elem = fe.dofs_per_cell;                      //This gives you dofs per element
  const unsigned int nodes_per_elem = GeometryInfo<dim>::vertices_per_cell;
  const unsigned int num_quad_pts = quadrature_formula.size();              //Total number of quad points in the element
  const unsigned int num_face_quad_pts = face_quadrature_formula.size();    //Total number of quad points in the face
  const unsigned int faces_per_elem = GeometryInfo<dim>::faces_per_cell;
  FullMatrix<double> Klocal (dofs_per_elem, dofs_per_elem);
  Vector<double>     Flocal (dofs_per_elem);

  std::vector<unsigned int> local_dof_indices (dofs_per_elem);              //This relates local dof numbering to global dof numbering

  //loop over elements  
  typename DoFHandler<dim>::active_cell_iterator elem = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (;elem!=endc; ++elem){
    //Update fe_values for the current element
    fe_values.reinit(elem);

    //Retrieve the effective "connectivity matrix" for this element
    elem->get_dof_indices (local_dof_indices);
		
    /*Global Assembly - 
      Get the current Flocal and Klocal from the functions you wrote above 
      and populate F_assembly and K_assembly using local_dof_indices to relate local and global DOFs.*/
    Klocal = 0.;

    //Loop over local DOFs and quadrature points to populate Klocal
    //Note that all quadrature points are included in this single loop
    for (unsigned int q=0; q<num_quad_pts; ++q){
      //evaluate elemental stiffness matrix, K^{AB}_{ik} = \integral N^A_{,j}*C_{ijkl}*N^B_{,l} dV 
      for (unsigned int A=0; A<nodes_per_elem; A++) { //Loop over nodes
	for(unsigned int i=0; i<dim; i++){ //Loop over nodal dofs
	  for (unsigned int B=0; B<nodes_per_elem; B++) {
	    for(unsigned int k=0; k<dim; k++){
	      for (unsigned int j = 0; j<dim; j++){
		for (unsigned int l = 0; l<dim; l++){
		  /*//EDIT - You need to define Klocal here. Note that the indices of Klocal are the element dof numbers (0 through 23),
		    which you can caluclate from the element node numbers (0 through 8) and the nodal dofs (0 through 2).
		    You'll need the following information:
		    basis gradient vector: fe_values.shape_grad(elementDOF,q), where elementDOF is dim*A+i or dim*B+k
		    NOTE: this is the gradient with respect to the real domain (not the bi-unit domain)
		    elasticity tensor: use the function C(i,j,k,l)
		    det(J) times the total quadrature weight: fe_values.JxW(q)*/
		Klocal[3*A+i][3*B+k] += fe_values.shape_grad(3*A+i,q)[j]*C(i,j,k,l)* fe_values.shape_grad(3*B+k,q)[l]*fe_values.JxW(q);
		}
	      }
	    }
	  }
	}
      }
    }

    Flocal = 0.;

    //Loop over faces (for Neumann BCs), local DOFs and quadrature points to populate Flocal.

    //If you had a forcing function, you would need to use FEValues here as well to integrate over the volume.

    //Add Neumann boundary conditions here in Flocal by integrating over the appropriate surface
    Vector<double> h(dim); h=0.;
    for (unsigned int f=0; f < faces_per_elem; f++){
      //Update fe_face_values from current element and face
      fe_face_values.reinit (elem, f);
      /*elem->face(f)->center() gives a position vector (in the real domain) of the center point on face f
	of the current element. We can use it to see if we are at the Neumann boundary, x_3 = 1.*/
      if(elem->face(f)->center()[2] == 1){
	//To integrate over this face, loop over all face quadrature points with this single loop
	for (unsigned int q=0; q<num_face_quad_pts; ++q){
	  double x = fe_face_values.quadrature_point(q)[0]; //x-coordinate at the current surface quad. point
	  //EDIT - define the value of the traction vector, h
	 h[2]  = 1.0e9*x;
	  for (unsigned int A=0; A<nodes_per_elem; A++){ //loop over all element nodes
	    for(unsigned int i=0; i<dim; i++){ //loop over nodal dofs
	      /*//EDIT - define Flocal. Again, the indices of Flocal are the element dof numbers (0 through 23).
		Evaluate the basis functions using the elementDOF: fe_face_values.shape_value(elementDOF,q)

		Note that we are looping over all element dofs, not just those on the Neumann face. However,
		the face quadrature points are only on the Neumann face, so we are indeed doing a surface integral.

		For det(J) times the total quadrature weight: fe_face_values.JxW(q)*/
		Flocal[3*A+i] += fe_face_values.shape_value(3*A+i,q)*fe_face_values.JxW(q)*h[i];
	    }
	  }
	}
      }
    }

    //Assemble local K and F into global K and F
    for(unsigned int i=0; i<dofs_per_elem; i++){
      //EDIT - Assemble F from Flocal (you can look at HW2)
      F[local_dof_indices[i]]  += Flocal[i];
      for(unsigned int j=0; j<dofs_per_elem; j++){
	//EDIT - Assemble K from Klocal (you can look at HW2)
        K.add(local_dof_indices[i],local_dof_indices[j],Klocal[i][j]);
      }
    }
  }

  //Let deal.II apply Dirichlet conditions WITHOUT modifying the size of K and F global
  MatrixTools::apply_boundary_values (boundary_values, K, D, F, false);
}

//Solve for D in KD=F
template <int dim>
void FEM<dim>::solve(){

  //Solve for D
  SparseDirectUMFPACK  A;
  A.initialize(K);
  A.vmult (D, F); //D=K^{-1}*F

}

//Output results
template <int dim>
void FEM<dim>::output_results (){

  //Write results to VTK file
  std::ofstream output1 ("solution.vtk");
  DataOut<dim> data_out; data_out.attach_dof_handler (dof_handler);

  //Add nodal DOF data
  data_out.add_data_vector (D,
			    nodal_solution_names,
			    DataOut<dim>::type_dof_data,
			    nodal_data_component_interpretation);
  data_out.build_patches();
  data_out.write_vtk(output1);
  output1.close();
}
