/**
 * @author Thomas Rückert
 * @date 2017
 */

#include <grpc++/grpc++.h>
#include <vector>
#include <math.h>

#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
//#include "vis.c"

#include "hypre.pb.h"
#include "hypre.grpc.pb.h"

#include "../../../hypre/src/utilities/_hypre_utilities.h"
#include "../../../hypre/src/krylov/HYPRE_krylov.h"
#include "../../../hypre/src/HYPRE.h"
#include "../../../hypre/src/parcsr_ls/HYPRE_parcsr_ls.h"
#include "../../../hypre/src/examples/vis.c"

using grpc::Status;
using grpc::ServerBuilder;
using grpc::Server;

using namespace ::rpc_hypre;

/**
 * think about this approach to store all parameters received
 */
class HypreContext {

  ::rpc_hypre::RPC_HYPRE_Solver solver;

  HYPRE_Solver hypreSolver;

  ::rpc_hypre::RPC_HYPRE_IJMatrix matrix;

  HYPRE_IJMatrix hypreMatrix;

//  HYPRE_ParCSRMatrix parcsr_A;
  std::vector<HYPRE_ParCSRMatrix> matrices;

//  HYPRE_ParVector par_b;
//  HYPRE_ParVector par_x;
  std::vector<HYPRE_ParVector> vectors;

};

/**
 * the rpc hypre service class
 */
class HypreService final : public rpc_hypre::hypreSrv::Service {

private:
protected:

  std::vector<HypreContext> contextList;

  std::vector<::rpc_hypre::RPC_HYPRE_Solver> solversList;
  std::vector<HYPRE_Solver> hypreSolversList;

  std::vector<::rpc_hypre::RPC_HYPRE_IJMatrix> matrixList;
  std::vector<HYPRE_IJMatrix> hypreMatrixList;

  std::vector<::rpc_hypre::RPC_HYPRE_ParCSRMatrix> parMatrixList;
  std::vector<HYPRE_ParCSRMatrix> hypreParMatrixList;

  /**
   * @return int
   */
  int getSolverIdentifier() {

    return hypreSolversList.size()-1;

  }

  /**
   * @return int
   */
  int getMatrixIdentifier() {

    return hypreMatrixList.size()-1;

  }

  /**
   * @return int
   */
  int getParMatrixIdentifier() {

    return hypreParMatrixList.size()-1;

  }

public:

  /* for MPI */
  int argc;
  char** argv;
  int myid, num_procs;

  /**
   * construct
   */
  HypreService() {
    argc = 0;
    argv = 0;
  }

  int* repeatedToArray(int size, const ::google::protobuf::RepeatedField< ::google::protobuf::int32> &repeated) {

    int* intArray = new int[size];
    for (int i = 0; i < size; i++) {

      intArray[i] = repeated[i];

    }
    return intArray;

  }

  double* repeatedToArray(int size, const ::google::protobuf::RepeatedField<double> &repeated) {

    double* doubleArray = new double[size];
    for (int i = 0; i < size; i++) {

      doubleArray[i] = repeated[i];

    }
    return doubleArray;

  }

  /**
   *
   * @param context
   * @param request
   * @param matrix
   * @return Status
   */
  Status RPC_HYPRE_IJMatrixCreate(::grpc::ServerContext* context,
                                  const ::rpc_hypre::RPC_HYPRE_IJMatrixCreateMessage* request,
                                  ::rpc_hypre::RPC_HYPRE_IJMatrix* matrix) override {

    HYPRE_IJMatrix hypreMatrix;
    HYPRE_IJMatrixCreate(
        MPI_COMM_WORLD, request->ilower(), request->iupper(), request->jlower(), request->jupper(), &hypreMatrix
    );

    //save the matrix
    matrixList.push_back(*matrix);
    hypreMatrixList.push_back(hypreMatrix);

    //get the id
    matrix->set_identifier(getMatrixIdentifier());

    return Status::OK;

  }

  /**
   *
   * @param context
   * @param request
   * @param response
   * @return Status
   */
  Status RPC_HYPRE_IJMatrixSetObjectType(::grpc::ServerContext* context,
                                         const ::rpc_hypre::RPC_HYPRE_GenericMatrixIntegerParam* request,
                                         ::rpc_hypre::RPC_HYPRE_IJMatrix* response) override {

    HYPRE_IJMatrix hypreMatrix = hypreMatrixList[request->matrix().identifier()];
    HYPRE_IJMatrixSetObjectType(hypreMatrix, request->value());

    return Status::OK;

  }

  /**
   *
   * @param context
   * @param request
   * @param response
   * @return Status
   */
  Status RPC_HYPRE_IJMatrixInitialize(::grpc::ServerContext* context,
                                      const ::rpc_hypre::RPC_HYPRE_IJMatrix* request,
                                      ::rpc_hypre::RPC_HYPRE_IJMatrix* response) override {

    HYPRE_IJMatrix hypreMatrix = hypreMatrixList[request->identifier()];
    HYPRE_IJMatrixInitialize(hypreMatrix);

    return Status::OK;

  }

  /**
   *
   * @param context
   * @param request
   * @param response
   * @return Status
   */
  Status RPC_HYPRE_IJMatrixSetValues(::grpc::ServerContext* context,
                                     const ::rpc_hypre::RPC_HYPRE_IJMatrixSetValuesMessage* request,
                                     ::rpc_hypre::RPC_HYPRE_IJMatrix* response) override {

    HYPRE_IJMatrix hypreMatrix = hypreMatrixList[request->matrix().identifier()];

    HYPRE_IJMatrixSetValues(
        hypreMatrix,
        request->nrows(),
        repeatedToArray(request->ncols_size(), request->ncols()),
        repeatedToArray(request->rows_size(), request->rows()),
        repeatedToArray(request->cols_size(), request->cols()),
        repeatedToArray(request->values_size(), request->values())
    );

    return Status::OK;

  }

  /**
   *
   * @param context
   * @param request
   * @param response
   * @return Status
   */
  Status RPC_HYPRE_IJMatrixAssemble(::grpc::ServerContext* context,
                                    const ::rpc_hypre::RPC_HYPRE_IJMatrix* request,
                                    ::rpc_hypre::RPC_HYPRE_IJMatrix* response) override {

    HYPRE_IJMatrix hypreMatrix = hypreMatrixList[request->identifier()];
    HYPRE_IJMatrixAssemble(hypreMatrix);

    return Status::OK;

  }

  /**
   *
   * @param context
   * @param request
   * @param response
   * @return Status
   */
  Status RPC_HYPRE_IJMatrixGetObject(::grpc::ServerContext* context,
                                     const ::rpc_hypre::RPC_HYPRE_GenericMatrixIntegerParam* request,
                                     ::rpc_hypre::RPC_HYPRE_IJMatrix* response) override {

    HYPRE_IJMatrix hypreMatrix = hypreMatrixList[request->matrix().identifier()];

    RPC_HYPRE_ParCSRMatrix parcsr_A;
    HYPRE_ParCSRMatrix hypreParcsr_A;

    HYPRE_IJMatrixGetObject(hypreMatrix, (void**) &hypreParcsr_A);

    parMatrixList.push_back(parcsr_A);
    hypreParMatrixList.push_back(hypreParcsr_A);

    return Status::OK;

  }

  /**
   *
   * @param context
   * @param request
   * @param response
   * @return Status
   */
  Status RPC_HYPRE_IJMatrixDestroy(::grpc::ServerContext* context,
                                   const ::rpc_hypre::RPC_HYPRE_IJMatrix* request,
                                   ::rpc_hypre::Empty* response) override {

    HYPRE_IJMatrix hypreMatrix = hypreMatrixList[request->identifier()];
    HYPRE_IJMatrixDestroy(hypreMatrix);

    return Status::OK;

  }


  /**
   *
   * @param context
   * @param request
   * @param response
   * @return
   */
  Status RPC_HYPRE_BoomerAMGCreate(::grpc::ServerContext* context,
                                   const ::rpc_hypre::Empty* request,
                                   ::rpc_hypre::RPC_HYPRE_Solver* solver) override {

    std::cout << "Running RPC_HYPRE_BoomerAMGCreate()." << std::endl;

//    ::rpc_hypre::RPC_HYPRE_Solver* solver = new ::rpc_hypre::RPC_HYPRE_Solver();
    HYPRE_Solver hypreSolver;

    //create the hypre solver
    HYPRE_BoomerAMGCreate(&hypreSolver);

    //save the solver
    solversList.push_back(*solver);
    hypreSolversList.push_back(hypreSolver);

    //get the id
    solver->set_identifier(getSolverIdentifier());

    std::cout << "Finished running RPC_HYPRE_BoomerAMGCreate()." << std::endl;

    return Status::OK;
  }

  /**
   *
   * @param context
   * @param request
   * @param response
   * @return
   */
  Status RPC_HYPRE_BoomerAMGSetup(::grpc::ServerContext* context,
                                  const ::rpc_hypre::RPC_HYPRE_BoomerAMGSolveMessage* request,
                                  ::rpc_hypre::Empty* response) {

    std::cout << "Running RPC_HYPRE_BoomerAMGSetup()." << std::endl;
    std::cout << "Using solver:" << request->solver().identifier() << std::endl;

//    ::rpc_hypre::RPC_HYPRE_ParVector* par_b = request->mutable_par_b();
//    ::rpc_hypre::RPC_HYPRE_ParVector* par_x = request->mutable_par_x();
//    ::rpc_hypre::RPC_HYPRE_ParCSRMatrix* parcsr_A = request->mutable_parcsr_a();

    //vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv just for testing vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    int i;
    int N, n;

    int ilower, iupper;
    int local_size, extra;

//    int solver_id;
//    int print_system;

    double h, h2;

    HYPRE_IJMatrix A;
    HYPRE_ParCSRMatrix parcsr_A;
    HYPRE_IJVector b;
    HYPRE_ParVector par_b;
    HYPRE_IJVector x;
    HYPRE_ParVector par_x;

    /* Default problem parameters */
    n = 33;

    /* Preliminaries: want at least one processor per row */
    if (n*n < num_procs) n = sqrt(num_procs) + 1;
    N = n*n; /* global number of rows */
    h = 1.0/(n+1); /* mesh size*/
    h2 = h*h;

    std::cout << "global number of rows: " << N << std::endl;
    std::cout << "mesh size: " << h << std::endl;

    /* Each processor knows only of its own rows - the range is denoted by ilower
       and upper.  Here we partition the rows. We account for the fact that
       N may not divide evenly by the number of processors. */
    local_size = N/num_procs;
    extra = N - local_size*num_procs;

    std::cout << "num procs: " << num_procs << std::endl;
    std::cout << "local size: " << local_size << std::endl;

    ilower = local_size*myid;
    ilower += hypre_min(myid, extra);

    iupper = local_size*(myid+1);
    iupper += hypre_min(myid+1, extra);
    iupper = iupper - 1;

    /* How many rows do I have? */
    local_size = iupper - ilower + 1;

    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);

    /* Choose a parallel csr format storage (see the User's Manual) */
    HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);

    /* Initialize before setting coefficients */
    HYPRE_IJMatrixInitialize(A);

    /* Now go through my local rows and set the matrix entries.
       Each row has at most 5 entries. For example, if n=3:

       A = [M -I 0; -I M -I; 0 -I M]
       M = [4 -1 0; -1 4 -1; 0 -1 4]

       Note that here we are setting one row at a time, though
       one could set all the rows together (see the User's Manual).
    */
    {
      int nnz;
      double values[5];
      int cols[5];

      std::cout << "ilower: " << ilower << std::endl;
      std::cout << "iupper: " << iupper << std::endl;

      for (i = ilower; i <= iupper; i++)
      {
        nnz = 0;

        /* The left identity block:position i-n */
        if ((i-n)>=0)
        {
          cols[nnz] = i-n;
          values[nnz] = -1.0;
          nnz++;
        }

        /* The left -1: position i-1 */
        if (i%n)
        {
          cols[nnz] = i-1;
          values[nnz] = -1.0;
          nnz++;
        }

        /* Set the diagonal: position i */
        cols[nnz] = i;
        values[nnz] = 4.0;
        nnz++;

        /* The right -1: position i+1 */
        if ((i+1)%n)
        {
          cols[nnz] = i+1;
          values[nnz] = -1.0;
          nnz++;
        }

        /* The right identity block:position i+n */
        if ((i+n)< N)
        {
          cols[nnz] = i+n;
          values[nnz] = -1.0;
          nnz++;
        }

        /* Set the values for row i */
        HYPRE_IJMatrixSetValues(A, 1, &nnz, &i, cols, values);
      }
    }

    /* Assemble after setting the coefficients */
    HYPRE_IJMatrixAssemble(A);

    /* Note: for the testing of small problems, one may wish to read
       in a matrix in IJ format (for the format, see the output files
       from the -print_system option).
       In this case, one would use the following routine:
       HYPRE_IJMatrixRead( <filename>, MPI_COMM_WORLD,
                           HYPRE_PARCSR, &A );
       <filename>  = IJ.A.out to read in what has been printed out
       by -print_system (processor numbers are omitted).
       A call to HYPRE_IJMatrixRead is an *alternative* to the
       following sequence of HYPRE_IJMatrix calls:
       Create, SetObjectType, Initialize, SetValues, and Assemble
    */


    /* Get the parcsr matrix object to use */
    HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);


    /* Create the rhs and solution */
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper,&b);
    HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(b);

    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper,&x);
    HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(x);

    /* Set the rhs values to h^2 and the solution to zero */
    {
      double *rhs_values, *x_values;
      int    *rows;

      rhs_values =  (double*) calloc(local_size, sizeof(double));
      x_values =  (double*) calloc(local_size, sizeof(double));
      rows = (int*) calloc(local_size, sizeof(int));

      for (i=0; i<local_size; i++)
      {
        rhs_values[i] = h2;
        x_values[i] = 0.0;
        rows[i] = ilower + i;
      }

      HYPRE_IJVectorSetValues(b, local_size, rows, rhs_values);
      HYPRE_IJVectorSetValues(x, local_size, rows, x_values);

      free(x_values);
      free(rhs_values);
      free(rows);
    }


    HYPRE_IJVectorAssemble(b);
    /*  As with the matrix, for testing purposes, one may wish to read in a rhs:
        HYPRE_IJVectorRead( <filename>, MPI_COMM_WORLD,
                                  HYPRE_PARCSR, &b );
        as an alternative to the
        following sequence of HYPRE_IJVectors calls:
        Create, SetObjectType, Initialize, SetValues, and Assemble
    */
    HYPRE_IJVectorGetObject(b, (void **) &par_b);

    HYPRE_IJVectorAssemble(x);
    HYPRE_IJVectorGetObject(x, (void **) &par_x);
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ just for testing ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ::rpc_hypre::RPC_HYPRE_Solver solver = request->solver();

    //there are no identifiers < 0
    if ((unsigned)solver.identifier() >= hypreSolversList.size()) {
      return Status::CANCELLED;
    }

    HYPRE_Solver hypreSolver = hypreSolversList[solver.identifier()];

    HYPRE_BoomerAMGSetup(hypreSolver, parcsr_A, par_b, par_x);

    //vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv just for testing vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

    HYPRE_BoomerAMGSolve(hypreSolver, parcsr_A, par_b, par_x);

    int num_iterations;
    double final_res_norm;

    /* Run info - needed logging turned on */
    HYPRE_BoomerAMGGetNumIterations(hypreSolver, &num_iterations);
    HYPRE_BoomerAMGGetFinalRelativeResidualNorm(hypreSolver, &final_res_norm);

    if (myid == 0)
    {
      printf("\n");
      printf("Iterations = %d\n", num_iterations);
      printf("Final Relative Residual Norm = %e\n", final_res_norm);
      printf("\n");
    }

    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ just for testing ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    std::cout << "Finished running RPC_HYPRE_BoomerAMGSetup()." << std::endl;

    return Status::OK;
  }

  /**
   *
   * @param context
   * @param request
   * @param response
   * @return
   */
  Status RPC_HYPRE_BoomerAMGSolve(::grpc::ServerContext* context,
                                  const ::rpc_hypre::RPC_HYPRE_BoomerAMGSolveMessage* request,
                                  ::rpc_hypre::RPC_HYPRE_ParVector* response) {

    std::cout << "Using solver:" << request->solver().identifier() << std::endl;

    ::rpc_hypre::RPC_HYPRE_Solver solver = request->solver();

    //there are no identifiers < 0
    if ((unsigned)solver.identifier() >= hypreSolversList.size()) {
      return Status::CANCELLED;
    }

    HYPRE_Solver hypreSolver = hypreSolversList[solver.identifier()];

//    HYPRE_BoomerAMGSolve(hypreSolver, parcsr_A, par_b, par_x);
//
//    int num_iterations;
//    double final_res_norm;
//
//    /* Run info - needed logging turned on */
//    HYPRE_BoomerAMGGetNumIterations(hypreSolver, &num_iterations);
//    HYPRE_BoomerAMGGetFinalRelativeResidualNorm(hypreSolver, &final_res_norm);
//
//    if (myid == 0)
//    {
//      printf("\n");
//      printf("Iterations = %d\n", num_iterations);
//      printf("Final Relative Residual Norm = %e\n", final_res_norm);
//      printf("\n");
//    }

    return Status::OK;

  }

};

void RunServer(int argc, char** argv) {

  int myid, num_procs;

  /* Initialize MPI */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  if (myid == 0)
  {

    std::string server_address("0.0.0.0:50051");
    HypreService service;

    service.argc = argc;
    service.argv = argv;
    service.myid = myid;
    service.num_procs = num_procs;

    ServerBuilder builder;
    // Listen on the given address without any authentication mechanism.
    builder.AddListeningPort(server_address, grpc::InsecureServerCredentials());
    // Register "service" as the instance through which we'll communicate with
    // clients. In this case it corresponds to an *synchronous* service.
    builder.RegisterService(&service);
    // Finally assemble the server.
    std::unique_ptr<Server> server(builder.BuildAndStart());
    std::cout << "Server listening on " << server_address << std::endl;

    // Wait for the server to shutdown. Note that some other thread must be
    // responsible for shutting down the server for this call to ever return.
    server->Wait();

  }

  MPI_Finalize();
}

int main(int argc, char** argv) {

  RunServer(argc, argv);

  return 0;
}
