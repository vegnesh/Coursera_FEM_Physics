/*This is a template file for use with 1D finite elements.
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

template <int dim>
class FEM
{
 public:
  //Class functions
  FEM (unsigned int order,unsigned int problem); // Class constructor 
  ~FEM(); //Class destructor

  //Function to find the value of xi at the given node (using deal.II node numbering)
  double xi_at_node(unsigned int dealNode);

  //Define your 1D basis functions and derivatives
  double basis_function(unsigned int node, double xi);
  double basis_gradient(unsigned int node, double xi);

  //Solution steps
  void generate_mesh(unsigned int numberOfElements);
  void define_boundary_conds();
  void setup_system();
  void assemble_system();
  void solve();
  void output_results();

  //Function to calculate the l2 norm of the error in the finite element sol'n vs. the exact solution
  double l2norm_of_error();

  //Class objects
  Triangulation<dim>   triangulation; //mesh
  FESystem<dim>        fe; 	      //FE element
  DoFHandler<dim>      dof_handler;   //Connectivity matrices

  //Gaussian quadrature - These will be defined in setup_system()
  unsigned int	        quadRule;    //quadrature rule, i.e. number of quadrature points
  std::vector<double>	quad_points; //vector of Gauss quadrature points
  std::vector<double>	quad_weight; //vector of the quadrature point weights
    
  //Data structures
  SparsityPattern       sparsity_pattern; //Sparse matrix pattern
  SparseMatrix<double>  K;		 //Global stiffness (sparse) matrix
  Vector<double>        D, F; 		 //Global vectors - Solution vector (D) and Global force vector (F)
  std::vector<double>   nodeLocation;	 //Vector of the x-coordinate of nodes by global dof number
  std::map<unsigned int,double> boundary_values;	//Map of dirichlet boundary conditions
  double                basisFunctionOrder, prob, L, g1, g2;    

  //solution name array
  std::vector<std::string> nodal_solution_names;
  std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation;
};

// Class constructor for a vector field
template <int dim>
FEM<dim>::FEM(unsigned int order,unsigned int problem)
:
fe (FE_Q<dim>(order), dim), 
  dof_handler (triangulation)
{
  basisFunctionOrder = order;
  if(problem == 1 || problem == 2){
    prob = problem;
  }
  else{
    std::cout << "Error: problem number should be 1 or 2.\n";
    exit(0);
  }

  //Nodal Solution names - this is for writing the output file
  for (unsigned int i=0; i<dim; ++i){
    nodal_solution_names.push_back("u");
    nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  }
}

//Class destructor
template <int dim>
FEM<dim>::~FEM(){
  dof_handler.clear();
}

//Find the value of xi at the given node (using deal.II node numbering)
template <int dim>
double FEM<dim>::xi_at_node(unsigned int dealNode){
  double xi;

  if(dealNode == 0){
    xi = -1.;
  }
  else if(dealNode == 1){
    xi = 1.;
  }
  else if(dealNode <= basisFunctionOrder){
    xi = -1. + 2.*(dealNode-1.)/basisFunctionOrder;
  }
  else{
    std::cout << "Error: you input node number "
	      << dealNode << " but there are only " 
	      << basisFunctionOrder + 1 << " nodes in an element.\n";
    exit(0);
  }

  return xi;
}

//Define basis functions
template <int dim>
double FEM<dim>::basis_function(unsigned int node, double xi){
  /*"basisFunctionOrder" defines the polynomial order of the basis function,
    "node" specifies which node the basis function corresponds to, 
    "xi" is the point (in the bi-unit domain) where the function is being evaluated.
    You need to calculate the value of the specified basis function and order at the given quadrature pt.*/

  double value = 1.; //Store the value of the basis function in this variable

  /*You can use the function "xi_at_node" (defined above) to get the value of xi (in the bi-unit domain)
    at any node in the element - using deal.II's element node numbering pattern.*/

 //FINISH EDIT
  double xi_b = xi_at_node(node);
  for(unsigned int index = 0; index <= basisFunctionOrder; index++)
 {
  if(index != node)
  value *= ((xi - xi_at_node(index))/(xi_b-xi_at_node(index)));
 }

  return value;
}

//Define basis function gradient
template <int dim>
double FEM<dim>::basis_gradient(unsigned int node, double xi){
  /*"basisFunctionOrder" defines the polynomial order of the basis function,
    "node" specifies which node the basis function corresponds to, 
    "xi" is the point (in the bi-unit domain) where the function is being evaluated.
    You need to calculate the value of the derivative of the specified basis function and order at the given quadrature pt.
    Note that this is the derivative with respect to xi (not x)*/

  double value = 0.; //Store the value of the gradient of the basis function in this variable

  /*You can use the function "xi_at_node" (defined above) to get the value of xi (in the bi-unit domain)
    at any node in the element - using deal.II's element node numbering pattern.*/

  //FINISH EDIT

  double xi_b = xi_at_node(node);
  
  for(unsigned int outer =0;outer<= basisFunctionOrder; outer++)
 {
  if(outer != node )
 {
  value += (basis_function(node,xi)/(xi - xi_at_node(outer)));

 }
  

 }
  return value;
}

//Define the problem domain and generate the mesh
template <int dim>
void FEM<dim>::generate_mesh(unsigned int numberOfElements){

  //Define the limits of your domain
  L = 0.1; //FINISH EDIT
  double x_min = 0.;
  double x_max = L;

  Point<dim,double> min(x_min),
    max(x_max);
  std::vector<unsigned int> meshDimensions (dim,numberOfElements);
  GridGenerator::subdivided_hyper_rectangle (triangulation, meshDimensions, min, max);
}

//Specify the Dirichlet boundary conditions
template <int dim>
void FEM<dim>::define_boundary_conds(){
  const unsigned int totalNodes = dof_handler.n_dofs(); //Total number of nodes

  //Identify dirichlet boundary nodes and specify their values.
  //This function is called from within "setup_system"

  /*The vector "nodeLocation" gives the x-coordinate in the real domain of each node,
    organized by the global node number.*/

  /*This loops through all nodes in the system and checks to see if they are
    at one of the boundaries. If at a Dirichlet boundary, it stores the node number
    and the applied displacement value in the std::map "boundary_values". Deal.II
    will use this information later to apply Dirichlet boundary conditions.
    Neumann boundary conditions are applied when constructing Flocal in "assembly"*/
  for(unsigned int globalNode=0; globalNode<totalNodes; globalNode++){
    if(nodeLocation[globalNode] == 0){
      boundary_values[globalNode] = g1;
    }
    if(nodeLocation[globalNode] == L){
      if(prob == 1){
	boundary_values[globalNode] = g2;
      }
    }
  }
			
}

//Setup data structures (sparse matrix, vectors)
template <int dim>
void FEM<dim>::setup_system(){

  //Define constants for problem (Dirichlet boundary values)
  g1 = 0.; g2 = 0.001; //FINISH EDIT

  //Let deal.II organize degrees of freedom
  dof_handler.distribute_dofs (fe);

  //Enter the global node x-coordinates into the vector "nodeLocation"
  MappingQ1<dim,dim> mapping;
  std::vector< Point<dim,double> > dof_coords(dof_handler.n_dofs());
  nodeLocation.resize(dof_handler.n_dofs());
  DoFTools::map_dofs_to_support_points<dim,dim>(mapping,dof_handler,dof_coords);
  for(unsigned int i=0; i<dof_coords.size(); i++){
    nodeLocation[i] = dof_coords[i][0];
  }

  //Specify boundary condtions (call the function)
  define_boundary_conds();

  //Define the size of the global matrices and vectors
  sparsity_pattern.reinit (dof_handler.n_dofs(), dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();
  K.reinit (sparsity_pattern);
  F.reinit (dof_handler.n_dofs());
  D.reinit (dof_handler.n_dofs());

  //Define quadrature rule
  /*A quad rule of 2 is included here as an example. You will need to decide
    what quadrature rule is needed for the given problems*/

  quadRule = 2; // FINISH EDIT - Number of quadrature points along one dimension
  quad_points.resize(quadRule); quad_weight.resize(quadRule);

  quad_points[0] = -sqrt(1./3.); //FINISH EDIT
  quad_points[1] = sqrt(1./3.); //FINISH EDIT

  quad_weight[0] = 1.; //FINISH EDIT
  quad_weight[1] = 1.; //FINISH EDIT

	

  quadRule = 3; //FINISH EDIT - Number of quadrature points along one dimension
  quad_points.resize(quadRule); quad_weight.resize(quadRule);

  quad_points[0] = 0.; //FINISH EDIT
  quad_points[1] = -sqrt(3./5.); //FINISH EDIT
  quad_points[2] = sqrt(3./5.); //FINISH EDIT

  quad_weight[0] = 8./9.; //FINISH EDIT
  quad_weight[1] = 5./9.; //FINISH EDIT
  quad_weight[2] = 5./9.; //FINISH EDIT


  std::cout<<"	Quad Order:	"<<quadRule<<std::endl;
  //Just some notes...
  std::cout << "   Number of active elems:       " << triangulation.n_active_cells() << std::endl;
  std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;   
}

//Form elmental vectors and matrices and assemble to the global vector (F) and matrix (K)
template <int dim>
void FEM<dim>::assemble_system(){

  K=0; F=0;

  const unsigned int   			dofs_per_elem = fe.dofs_per_cell; //This gives you number of degrees of freedom per element
  FullMatrix<double> 				Klocal (dofs_per_elem, dofs_per_elem);
  Vector<double>      			Flocal (dofs_per_elem);
  std::vector<unsigned int> local_dof_indices (dofs_per_elem);
  double										h_e, x, f;

  //loop over elements  
  typename DoFHandler<dim>::active_cell_iterator elem = dof_handler.begin_active(), 
    endc = dof_handler.end();
  for (;elem!=endc; ++elem){

    /*Retrieve the effective "connectivity matrix" for this element
      "local_dof_indices" relates local dofs to global dofs,
      i.e. local_dof_indices[i] gives the global dof number for local dof i.*/
    elem->get_dof_indices (local_dof_indices);

    /*We find the element length by subtracting the x-coordinates of the two end nodes
      of the element. Remember that the vector "nodeLocation" holds the x-coordinates, indexed
      by the global node number. "local_dof_indices" gives us the global node number indexed by
      the element node number.*/
    h_e = nodeLocation[local_dof_indices[1]] - nodeLocation[local_dof_indices[0]];

    //Loop over local DOFs and quadrature points to populate Flocal and Klocal.
    double Area = 1.e-4;
    Flocal = 0.;
    f = 1.e11;
    for(unsigned int A=0; A<dofs_per_elem; A++){
      for(unsigned int q=0; q<quadRule; q++){
	x = 0;
	//Interpolate the x-coordinates at the nodes to find the x-coordinate at the quad pt.
	for(unsigned int B=0; B<dofs_per_elem; B++){
	  x += nodeLocation[local_dof_indices[B]]*basis_function(B,quad_points[q]);
	}
	//FINISH EDIT - Define Flocal.
	Flocal[A] += basis_function(A,quad_points[q])*quad_weight[q]*h_e*0.5*f*Area*x; 

      }
//    std::cout<<"\nFlocal: "<<f*Area*h_e*0.5 <<"\t"<< Flocal[A];
    }
    //Add nonzero Neumann condition, if applicable
    double hval = 1.e6;
    if(prob == 2){ 
      if(nodeLocation[local_dof_indices[1]] == L){
	//FINISH EDIT VERIFY - Modify Flocal to include the traction on the right boundary.
//	std::cout<<"\n"<<Flocal[1]<<" This is Fla1\n";

	Flocal[1] += basis_function(1,1.0)*hval;
//	std::cout<<"\n"<<Flocal[1]<<" This is Fla2\n";
      }
    }

    //Loop over local DOFs and quadrature points to populate Klocal
    Klocal = 0;
    double Eval = 1.e11;
    for(unsigned int A=0; A<dofs_per_elem; A++){
      for(unsigned int B=0; B<dofs_per_elem; B++){
	for(unsigned int q=0; q<quadRule; q++){
	  //FINISH EDIT - Define Klocal.
	  Klocal[A][B] += (2*Eval/h_e)*basis_gradient(A,quad_points[q])*basis_gradient(B,quad_points[q])*quad_weight[q]*Area;
//	  std::cout<<"\nBasis Gradient:"<<basis_gradient(A,quad_points[q])<<"\t"<<basis_gradient(B,quad_points[q]);
	}
      }
  //  std::cout<<"\n"<<Eval*Area/h_e<<"\t"<<Klocal[A][0]<<"\t"<<Klocal[A][1];
     }

    //Assemble local K and F into global K and F
    //You will need to used local_dof_indices[A]
    for(unsigned int A=0; A<dofs_per_elem; A++){
      //FINISH EDIT - add component A of Flocal to the correct location in F
      /*Remember, local_dof_indices[A] is the global degree-of-freedom number
	corresponding to element node number A*/
	F[local_dof_indices[A]]  += Flocal[A];
//	F[0] = 100000.0;
//	if(nodeLocation[local_dof_indices[1]] == L){
//	 F[local_dof_indices[1]] = 100000.0;
//	}

//        std::cout<<"F:\t"<<local_dof_indices[A]<<"\t"<<F[local_dof_indices[A]]<<"\t"<<Flocal[A]<<"\n";
      for(unsigned int B=0; B<dofs_per_elem; B++){
	//FINISH EDIT - add component A,B of Klocal to the correct location in K (using local_dof_indices)
	/*Note: K is a sparse matrix, so you need to use the function "add".
	  For example, to add the variable C to K[i][j], you would use:
	  K.add(i,j,C);*/
	  K.add(local_dof_indices[A],local_dof_indices[B],Klocal[A][B]);
//	  std::cout<<"K:\t"<<local_dof_indices[A]<<"\t"<<local_dof_indices[B]<<"\t"<<Klocal[A][B]<<"\n";
      }
    }

  }

  //Apply Dirichlet boundary conditions
  /*deal.II applies Dirichlet boundary conditions (using the boundary_values map we
    defined in the function "define_boundary_conds") without resizing K or F*/
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
  std::ofstream output1("solution.vtk");
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  //Add nodal DOF data
  data_out.add_data_vector(D, nodal_solution_names, DataOut<dim>::type_dof_data,
			   nodal_data_component_interpretation);
  data_out.build_patches();
  data_out.write_vtk(output1);
  output1.close();
}

template <int dim>
double FEM<dim>::l2norm_of_error(){
	
  double l2norm = 0.;

  //Find the l2 norm of the error between the finite element sol'n and the exact sol'n
  const unsigned int   			dofs_per_elem = fe.dofs_per_cell; //This gives you dofs per element
  std::vector<unsigned int> local_dof_indices (dofs_per_elem);
  double u_exact, u_h, x, h_e;

  //loop over elements  
  typename DoFHandler<dim>::active_cell_iterator elem = dof_handler.begin_active (), 
    endc = dof_handler.end();
  for (;elem!=endc; ++elem){

    //Retrieve the effective "connectivity matrix" for this element
    elem->get_dof_indices (local_dof_indices);

    //Find the element length
    h_e = nodeLocation[local_dof_indices[1]] - nodeLocation[local_dof_indices[0]];
    double cval, L = 0.1;
    if(prob==1)
{
//    cval = 1./100. + L*L/6;
      cval  = 7.*L*L/6.;
}
else
{
    cval = 1./10. + L*L/2.;
}
    for(unsigned int q=0; q<quadRule; q++){
      x = 0.; u_h = 0.;
      //Find the values of x and u_h (the finite element solution) at the quadrature points
      for(unsigned int B=0; B<dofs_per_elem; B++){
	x += nodeLocation[local_dof_indices[B]]*basis_function(B,quad_points[q]);
	u_h += D[local_dof_indices[B]]*basis_function(B,quad_points[q]);
      }
      //EDIT - Find the l2-norm of the error through numerical integration.
      /*This includes evaluating the exact solution at the quadrature points*/
	u_exact = -x*x*x/6. + cval*x;
	l2norm += (u_exact - u_h)*(u_exact-u_h)*h_e*0.5*quad_weight[q];				
    }
  }

  return sqrt(l2norm);
}
