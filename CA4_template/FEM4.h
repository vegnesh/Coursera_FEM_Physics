/*This is a template file for use with 3D finite elements.
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
  FEM (double Alpha); // Class constructor 
  ~FEM(); //Class destructor

  //Solution steps
  void generate_mesh(std::vector<unsigned int> numberOfElements);
  void define_boundary_conds();
  void setup_system();
  void assemble_system();

  void solve_steady();
  void apply_initial_conditions();
  void solve_trans();
  void output_steady_results();
  void output_trans_results(unsigned int index);

  //Calculate the l2norm of the difference between the steady state and transient solution
  double l2norm();

  //Class objects
  Triangulation<dim> triangulation; //mesh
  FESystem<dim>      fe;            //FE element
  DoFHandler<dim>    dof_handler;   // Connectivity matrices

  QGauss<dim>  	     quadrature_formula; //Quadrature

  //Data structures
  SparsityPattern      sparsity_pattern;                   //Sparse matrix pattern
  SparseMatrix<double> M, K, system_matrix;                //Global stiffness matrix - Sparse matrix - used in the solver
  Vector<double>       D_steady, D_trans, V_trans, F, RHS; //Global vectors - Solution vector (D) and Global force vector (F)

  Table<2,double>	        nodeLocation;	      //Table of the coordinates of nodes by global dof number
  std::map<unsigned int,double> boundary_values_of_D; //Map of dirichlet boundary conditions for the temperature
  std::map<unsigned int,double> boundary_values_of_V; //Map of dirichlet boundary conditions for the time derivative of temperature

  std::vector<double> l2norm_results; //A vector to store the l2norms calculated in the time loop in solve_trans()
  double	      alpha; 	      //Specifies the Euler method, 0 <= alpha <= 1
    
  //solution name array
  std::vector<std::string> nodal_solution_names;
  std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation;
};

// Class constructor for a scalar field
template <int dim>
FEM<dim>::FEM (double Alpha)
:
fe (FE_Q<dim>(order), 1),
  dof_handler (triangulation),
  quadrature_formula(quadRule)
{
  alpha = Alpha;

  nodal_solution_names.push_back("D");
  nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
}

//Class destructor
template <int dim>
FEM<dim>::~FEM (){dof_handler.clear ();}

//Define the problem domain and generate the mesh
template <int dim>
void FEM<dim>::generate_mesh(std::vector<unsigned int> numberOfElements){

  //Define the limits of your domain
  double x_min = 0.0, //EDIT - define the left limit of the domain, etc.
    x_max = 1.0, //EDIT
    y_min = 0.0, //EDIT
    y_max = 1.0, //EDIT
    z_min = 0.0, //EDIT
    z_max = 0.1; //EDIT

  Point<dim,double> min(x_min,y_min,z_min),
    max(x_max,y_max,z_max);
  GridGenerator::subdivided_hyper_rectangle (triangulation, numberOfElements, min, max);
}

//Specify the Dirichlet boundary conditions
template <int dim>
void FEM<dim>::define_boundary_conds(){

  //EDIT - Define the Dirichlet boundary conditions.
	
  /*Note: this will be very similiar to the define_boundary_conds function
    in the HW2 template. You will loop over all nodes and use "nodeLocations"
    to check if the node is on the boundary with a Dirichlet condition. If it is,
    then add the node number and the specified value
    to the boundary values map, something like this:

    boundary_values_of_D[globalNodeIndex] = dirichletTemperatureValue
    boundary_values_of_V[globalNodeIndex] = dirichletTemperatureTimeDerivative

    Note that there are boundary conditions on temperature and on the time derivative of temperature.
    "nodeLocation" is a Table. The row index is
    the global node number; the column index refers to the x, y, or z component (0, 1, or 2 for 3D).
    e.g. nodeLocation[7][2] is the z coordinate of global node 7*/

  const unsigned int totalNodes = dof_handler.n_dofs(); //Total number of nodes
for(unsigned int globalNode=0; globalNode<totalNodes; globalNode++)
{
 if(nodeLocation[globalNode][0]==0.0)
{
 boundary_values_of_D[globalNode] = 300.0; 
 boundary_values_of_V[globalNode] = 0.0; 

}
 if(nodeLocation[globalNode][0]==1.0)
{
 boundary_values_of_D[globalNode] = 310.0;
 boundary_values_of_V[globalNode] = 0.0; 

}
}
}

//Setup data structures (sparse matrix, vectors)
template <int dim>
void FEM<dim>::setup_system(){

  //Let deal.II organize degrees of freedom
  dof_handler.distribute_dofs (fe);

  //Fill in the Table "nodeLocations" with the x, y, and z coordinates of each node by its global index
  MappingQ1<dim,dim> mapping;
  std::vector< Point<dim,double> > dof_coords(dof_handler.n_dofs());
  nodeLocation.reinit(dof_handler.n_dofs(),dim);
  DoFTools::map_dofs_to_support_points<dim,dim>(mapping,dof_handler,dof_coords);
  for(unsigned int i=0; i<dof_coords.size(); i++){
    for(unsigned int j=0; j<dim; j++){
      nodeLocation[i][j] = dof_coords[i][j];
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
  M.reinit (sparsity_pattern);
  system_matrix.reinit (sparsity_pattern);
  D_steady.reinit(dof_handler.n_dofs());
  D_trans.reinit(dof_handler.n_dofs());
  V_trans.reinit(dof_handler.n_dofs());
  RHS.reinit(dof_handler.n_dofs());
  F.reinit(dof_handler.n_dofs());

  //Just some notes...
  std::cout << "   Number of active elems:       " << triangulation.n_active_cells() << std::endl;
  std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;   
}

//Form elmental vectors and matrices and assemble to the global vector (F) and matrix (K)
template <int dim>
void FEM<dim>::assemble_system(){

  M=0; K=0; F=0;

  FEValues<dim> fe_values(fe,
			  quadrature_formula,
			  update_values | 
			  update_gradients | 
			  update_JxW_values);

  const unsigned int dofs_per_elem = fe.dofs_per_cell;         //This gives you dofs per element
  unsigned int 	     num_quad_pts = quadrature_formula.size(); //Total number of quad points in the element
  FullMatrix<double> Mlocal (dofs_per_elem, dofs_per_elem);
  FullMatrix<double> Klocal (dofs_per_elem, dofs_per_elem);
  Vector<double>     Flocal (dofs_per_elem);

  std::vector<unsigned int> local_dof_indices (dofs_per_elem); //This relates local dof numbering to global dof numbering
  double		    rho = 3.8151e6 ;                            //EDIT - specify the specific heat per unit volume

  //loop over elements  
  typename DoFHandler<dim>::active_cell_iterator elem = dof_handler.begin_active (),
    endc = dof_handler.end();
  for (;elem!=endc; ++elem){

    //Retrieve the effective "connectivity matrix" for this element
    elem->get_dof_indices (local_dof_indices);

    fe_values.reinit(elem); //Retrieve values from current element
    /*FEValues will give you the basis functions, gradients, jacobian, etc. that you need at each quad point,
      just as in HW4:
      Basis function value with: fe_values.shape_value(element_dof_number,quad_point_number);
      Basis function gradient with: fe_values.shape_grad(element_dof_number,quad_point_number)[component];
      (again, this is the gradient with respect to the real (not bi-unit) domain)
      Determinant of the Jacobian times the quadrature weight with: fe_values.JxW(quad_point_number);*/

    elem->get_dof_indices (local_dof_indices);

    //Loop over local DOFs and quadrature points to populate Mlocal
    Mlocal = 0.;
    for(unsigned int q=0; q<num_quad_pts; q++){
      for(unsigned int A=0; A<fe.dofs_per_cell; A++){
	for(unsigned int B=0; B<fe.dofs_per_cell; B++){
	  //EDIT - define Mlocal[A][B]
	Mlocal[A][B] += fe_values.shape_value(A,q)*fe_values.shape_value(B,q)*fe_values.JxW(q)*rho;
	}
      }
    }

    FullMatrix<double> kappa(dim,dim);
    kappa[0][0] = 385.;
    kappa[1][1] = 385.;
    kappa[2][2] = 385.;

    //Loop over local DOFs and quadrature points to populate Klocal
    Klocal = 0.;
    for(unsigned int A=0; A<fe.dofs_per_cell; A++){
      for(unsigned int B=0; B<fe.dofs_per_cell; B++){
	for(unsigned int q=0; q<num_quad_pts; q++){
	  for(unsigned int i=0; i<dim; i++){
	    for(unsigned int j=0; j<dim; j++){
	      //EDIT - define Klocal[A][B]
		Klocal[A][B] += fe_values.shape_grad(A,q)[i]*kappa[i][j]*fe_values.shape_grad(B,q)[j]*fe_values.JxW(q);
	    }
	  }
	}
      }
    }

    for (unsigned int i=0; i<dofs_per_elem; ++i){
      for (unsigned int j=0; j<dofs_per_elem; ++j){
	//EDIT - assemble K and M from Klocal and Mlocal
	M.add(local_dof_indices[i],local_dof_indices[j],Mlocal[i][j]);
	K.add(local_dof_indices[i],local_dof_indices[j],Klocal[i][j]);
      }
    }
  }
}

//Solve for D_steady in KD=F
template <int dim>
void FEM<dim>::solve_steady(){

  //Note that I'm applying boundary conditions within each solve function
  //Let deal.II apply Dirichlet conditions
  MatrixTools::apply_boundary_values (boundary_values_of_D, K, D_steady, F, false);

  //Solve for D
  SparseDirectUMFPACK  A;
  A.initialize(K);
  A.vmult (D_steady, F); //D=K^{-1}*F

  //Call output results directly from the solve functions
  output_steady_results();
}

//Apply initial conditions for the transient problem
template <int dim>
void FEM<dim>::apply_initial_conditions(){

  /*Loop over global nodes. Use nodeLocation to determine the position of the node
    and specify the correct initial temperature value in D_trans.*/

  const unsigned int totalNodes = dof_handler.n_dofs(); //Total number of nodes

  for(unsigned int i=0; i<totalNodes; i++){
    if(nodeLocation[i][0] < 0.5){
      D_trans[i] = 300.0 ; //EDIT
    }
    else{
      D_trans[i] = 300.0 + 20.0*(nodeLocation[i][0]-0.5); //EDIT
    }
  }

  //Find V_0 = M^{-1}*(F_0 - K*D_0)
  system_matrix.copy_from(M); //system_matrix = M

  //We need to define the right-hand-side vector (RHS = F_0 - K*D_0) step by step...
  K.vmult(RHS,D_trans); //RHS = K*D_trans
  RHS *= -1.; 		//RHS = -1.*RHS = -K*D_trans
  RHS.add(1.,F); 	//RHS = RHS + 1.*F = F - K*D_trans
  MatrixTools::apply_boundary_values (boundary_values_of_V, system_matrix, V_trans, RHS, false);

  //Solve for solution in system_matrix*V_trans = RHS
  SparseDirectUMFPACK  A;
  A.initialize(system_matrix);
  A.vmult (V_trans, RHS); //V_trans=system_matrix^{-1}*RHS

	//Output initial state
  output_trans_results(0);

  double current_l2norm = l2norm();
  l2norm_results.push_back(current_l2norm);

}

//Solve for D_transient
template <int dim>
void FEM<dim>::solve_trans(){

  //Call the function to initialize D_trans and V_trans as D_0 and V_0
  apply_initial_conditions();

  //Define delta_t (leave at 1. when you turn in your assignment)
  const double delta_t = 1.;

  const unsigned int totalNodes = dof_handler.n_dofs(); //Total number of nodes
  Vector<double>     D_tilde(totalNodes);

  //Loop over time steps. For each time step, update D_transient from D_n to D_{n+1} using the V method
  for(unsigned int t_step=1; t_step<3001; t_step++){
	
    /*SparseMatrix and Vector operations in deal.II:
      To add a matrix with a coefficient to system matrix, use system_matrix.add(coeffient,matrix)
      To add a vector with a coefficient to RHS, use RHS.add(coefficient,vector)
      To multiply a sparse matrix and a vector, for example RHS = K*D_tilde, use
      K.vmult(RHS,D_tilde)
      To set system_matrix equal to M, use system_matrix.copy_from(M)
      To set RHS equal to F, use RHS = F
      To set RHS equal to -RHS, use RHS *= -1
      For some examples, look at the apply_initial_conditions() function*/

    //Find D_tilde. Remember, at this point D_trans = D_n and V_trans = V_n
    //EDIT
	D_tilde = D_trans;
	D_tilde.add(delta_t*(1.0-alpha),V_trans);
    /*Use D_tilde to update V_trans from V_n to V_{n+1}. This involves solving 
      a matrix/vector system: system_matrix*Vtrans = RHS. You need to define
      system_matrix and RHS to correctly solve for V_trans = V_{n+1} = system_matrix^{-1}*RHS*/
    //EDIT
       system_matrix.copy_from(M);
       system_matrix.add(alpha*delta_t,K);
       K.vmult(RHS,D_tilde); //RHS = K*D_trans
       RHS *= -1.;           //RHS = -1.*RHS = -K*D_trans
       RHS.add(1.,F);

    //Apply boundary conditions on V_trans before solving the matrix/vector system
    MatrixTools::apply_boundary_values (boundary_values_of_V, system_matrix, V_trans, RHS, false);

    //Solve for V_trans (V_{n+1}) in system_matrix*solution = RHS
    SparseDirectUMFPACK  A;
    A.initialize(system_matrix);
    A.vmult (V_trans, RHS); //V_trans=system_matrix^{-1}*RHS

    /*To clarify, .vmult with a SparseDirectUMFPACK matrix (like A)
      multiplies a vector by the matrix inverse,
      but .vmult with a SparseMatrix (like K or M)
      multiplies a vector by the matrix itself (non-inverted)*/

    //Update D_trans to D_{n+1} using D_tilde and V_trans (V_{n+1})
    //EDIT
    D_trans = D_tilde;
    D_trans.add(alpha*delta_t,V_trans);

    //Output the results every 100 seconds
    if(t_step%100 == 0){
      output_trans_results(t_step);

      double current_l2norm = l2norm();
      l2norm_results.push_back(current_l2norm);
    }
  }
}

//Output transient results for a given time step
template <int dim>
void FEM<dim>::output_steady_results (){
  //Write results to VTK file
  std::ofstream output1 ("solution.vtk");
  DataOut<dim> data_out; data_out.attach_dof_handler (dof_handler);

  //Add nodal DOF data
  data_out.add_data_vector (D_steady, nodal_solution_names, DataOut<dim>::type_dof_data, nodal_data_component_interpretation);
  data_out.build_patches (); data_out.write_vtk (output1); output1.close();
}

//Output transient results for a given time step
template <int dim>
void FEM<dim>::output_trans_results (unsigned int index){
  //This adds an index to your filename so that you can distinguish between time steps

  //Write results to VTK file
  char filename[100];
  snprintf(filename, 100, "solution_%d.vtk", index);
  std::ofstream output1 (filename);
  DataOut<dim> data_out; data_out.attach_dof_handler (dof_handler);

  //Add nodal DOF data
  data_out.add_data_vector (D_trans, nodal_solution_names, DataOut<dim>::type_dof_data, nodal_data_component_interpretation);
  data_out.build_patches (); data_out.write_vtk (output1); output1.close();
}

//Function to calculate the l2norm of the difference between the current and steady state solutions.
template <int dim>
double FEM<dim>::l2norm(){
  double l2norm = 0.;

  FEValues<dim> fe_values(fe,
			  quadrature_formula,
			  update_values |
			  update_JxW_values);

  const unsigned int 				dofs_per_elem = fe.dofs_per_cell; //This gives you dofs per element
  std::vector<unsigned int> local_dof_indices (dofs_per_elem);
  const unsigned int 				num_quad_pts = quadrature_formula.size(); //Total number of quad points in the element
  double 										u_steady, u_trans;

  //loop over elements  
  typename DoFHandler<dim>::active_cell_iterator elem = dof_handler.begin_active (), 
    endc = dof_handler.end();
  for (;elem!=endc; ++elem){
    elem->get_dof_indices (local_dof_indices);
    fe_values.reinit(elem);

    for(unsigned int q=0; q<num_quad_pts; q++){
      u_steady = 0.; u_trans = 0.;
      for(unsigned int A=0; A<dofs_per_elem; A++){
	/*//EDIT - interpolate the steady state solution (u_steady) and transient solution (u_trans)
	  at the current quadrature point using D_steady and D_trans. Similar to finding u_h in HW2*/
//	x += nodeLocation[local_dof_indices[A]]*fe_values.shape_value(A,q);
        u_steady += D_steady[local_dof_indices[A]]*fe_values.shape_value(A,q);
	u_trans += D_trans[local_dof_indices[A]]*fe_values.shape_value(A,q); 
      }
      //EDIT - define the l2norm of the difference between u_steady and u_trans
    l2norm += (u_steady-u_trans)*(u_steady-u_trans)*fe_values.JxW(q);
    }

  }
	
  //The l2norm is the square root of the integral of (u_steady - u_trans)^2, using D_steady and D_trans

  return sqrt(l2norm);
}
