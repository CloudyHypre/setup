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

using grpc::Channel;
using grpc::ClientContext;
using grpc::Status;
using rpc_hypre::hypreSrv;

class HypreClient {

protected:

  int solverIdentifier;
  int matrixIdentifier;

public:

  /**
   * construct
   * @param channel
   */
  HypreClient(std::shared_ptr<Channel> channel) : stub_(hypreSrv::NewStub(channel)) {}

//  ~HypreClient() {}

  Status createSolver() {

    //create

    ClientContext context;
    ::rpc_hypre::Empty emptyRequest;
    rpc_hypre::RPC_HYPRE_Solver solverMessage;

    std::cout << "Running RPC_HYPRE_BoomerAMGCreate()." << std::endl;
    Status status = stub_->RPC_HYPRE_BoomerAMGCreate(&context, emptyRequest, &solverMessage);
    std::cout << "Finished running RPC_HYPRE_BoomerAMGCreate()." << std::endl;

    std::cout << "Got solver: " << solverMessage.identifier() << std::endl;

    solverIdentifier = solverMessage.identifier();

    return status;

  }

  Status setupSolver() {

    //set it up

    rpc_hypre::RPC_HYPRE_BoomerAMGSolveMessage* request = new rpc_hypre::RPC_HYPRE_BoomerAMGSolveMessage;

    rpc_hypre::RPC_HYPRE_ParCSRMatrix* matrixA = new rpc_hypre::RPC_HYPRE_ParCSRMatrix;
    rpc_hypre::RPC_HYPRE_ParVector* par_b = new rpc_hypre::RPC_HYPRE_ParVector;
    rpc_hypre::RPC_HYPRE_ParVector* par_x = new rpc_hypre::RPC_HYPRE_ParVector;
    rpc_hypre::RPC_HYPRE_Solver* solverMessage = new rpc_hypre::RPC_HYPRE_Solver;

    solverMessage->set_identifier(solverIdentifier);
    matrixA->set_identifier(matrixIdentifier);

    request->set_allocated_solver(solverMessage);
    request->set_allocated_parcsr_a(matrixA);
    request->set_allocated_par_b(par_b);
    request->set_allocated_par_x(par_x);

    ::rpc_hypre::Empty emptyResponse;

    ClientContext context;
    std::cout << "Running RPC_HYPRE_BoomerAMGSetup()." << std::endl;
    Status status = stub_->RPC_HYPRE_BoomerAMGSetup(&context, *request, &emptyResponse);
    std::cout << "Finished running RPC_HYPRE_BoomerAMGSetup()." << std::endl;

    delete request;

    //results in seg fault:
//    delete matrixA;
//    delete par_b;
//    delete par_x;
//    delete solver2;

    return status;

  }

  Status solve() {

    //solve

    ClientContext context;
    rpc_hypre::RPC_HYPRE_BoomerAMGSolveMessage* request = new rpc_hypre::RPC_HYPRE_BoomerAMGSolveMessage;

    rpc_hypre::RPC_HYPRE_Solver* solver = new rpc_hypre::RPC_HYPRE_Solver;
    rpc_hypre::RPC_HYPRE_ParCSRMatrix* parMatrix = new rpc_hypre::RPC_HYPRE_ParCSRMatrix;
    rpc_hypre::RPC_HYPRE_ParVector* par_b = new rpc_hypre::RPC_HYPRE_ParVector;
    rpc_hypre::RPC_HYPRE_ParVector* par_x = new rpc_hypre::RPC_HYPRE_ParVector;

    solver->set_identifier(solverIdentifier);
    parMatrix->set_identifier(matrixIdentifier);

    request->set_allocated_solver(solver);
    request->set_allocated_parcsr_a(parMatrix);
    request->set_allocated_par_b(par_b);
    request->set_allocated_par_x(par_x);

    rpc_hypre::RPC_HYPRE_ParVector response;

    std::cout << "Running RPC_HYPRE_BoomerAMGSolve()." << std::endl;
    Status status = stub_->RPC_HYPRE_BoomerAMGSolve(&context, *request, &response);
    std::cout << "Finished running RPC_HYPRE_BoomerAMGSolve()." << std::endl;

    return status;

  }

  Status setupMatrix() {

    std::cout << "Running matrix()" << std::endl;

    int i;
    int N, n;

    int ilower, iupper;
    int local_size, extra;

    int num_procs = 1;
    int myid = 1;

//    int solver_id;
//    int print_system;

//    double h, h2;

//    HYPRE_IJMatrix A;
//    HYPRE_ParCSRMatrix parcsr_A;
//    HYPRE_IJVector b;
//    HYPRE_ParVector par_b;
//    HYPRE_IJVector x;
//    HYPRE_ParVector par_x;

    /* Default problem parameters */
    n = 33;

    /* Preliminaries: want at least one processor per row */
    if (n*n < num_procs) n = sqrt(num_procs) + 1;
    N = n*n; /* global number of rows */
//    h = 1.0/(n+1); /* mesh size*/
//    h2 = h*h;

//    std::cout << "global number of rows: " << N << std::endl;
//    std::cout << "mesh size: " << h << std::endl;

    /* Each processor knows only of its own rows - the range is denoted by ilower
       and upper.  Here we partition the rows. We account for the fact that
       N may not divide evenly by the number of processors. */
    local_size = N/num_procs;
    extra = N - local_size*num_procs;

//    std::cout << "num procs: " << num_procs << std::endl;
//    std::cout << "local size: " << local_size << std::endl;

    ilower = local_size*myid;
    ilower += hypre_min(myid, extra);

    iupper = local_size*(myid+1);
    iupper += hypre_min(myid+1, extra);
    iupper = iupper - 1;

    /* How many rows do I have? */
    local_size = iupper - ilower + 1;


    /// ----------------------------------- CREATE -----------------------------------
    //    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);

    ClientContext createContext;
    ::rpc_hypre::RPC_HYPRE_IJMatrixCreateMessage* createRequest = new ::rpc_hypre::RPC_HYPRE_IJMatrixCreateMessage;
    createRequest->set_ilower(ilower);
    createRequest->set_iupper(iupper);
    createRequest->set_jlower(ilower);
    createRequest->set_jupper(iupper);
    ::rpc_hypre::RPC_HYPRE_IJMatrix matrix;
    Status status = stub_->RPC_HYPRE_IJMatrixCreate(&createContext, *createRequest, &matrix);

    if (!status.ok()) {
      return status;
    }

    /* Choose a parallel csr format storage (see the User's Manual) */
    /// ----------------------------------- SET TYPE -----------------------------------

//    HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);

    ClientContext setTypeContext;
    ::rpc_hypre::RPC_HYPRE_GenericMatrixIntegerParam* setTypeRequest =
        new ::rpc_hypre::RPC_HYPRE_GenericMatrixIntegerParam;
    setTypeRequest->set_value(HYPRE_PARCSR);
    setTypeRequest->set_allocated_matrix(&matrix);
    status = stub_->RPC_HYPRE_IJMatrixSetObjectType(&setTypeContext, *setTypeRequest, &matrix);

    if (!status.ok()) {
      return status;
    }

    /* Initialize before setting coefficients */

    /// ----------------------------------- INIT -----------------------------------

//    HYPRE_IJMatrixInitialize(A);

    ClientContext initializeContext;
    ::rpc_hypre::RPC_HYPRE_IJMatrix* initRequest = new ::rpc_hypre::RPC_HYPRE_IJMatrix;
    initRequest->set_identifier(matrix.identifier());

    status = stub_->RPC_HYPRE_IJMatrixInitialize(&initializeContext, *initRequest, &matrix);

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

//      std::cout << "ilower: " << ilower << std::endl;
//      std::cout << "iupper: " << iupper << std::endl;

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

        /// ----------------------------------- SET -----------------------------------
//        HYPRE_IJMatrixSetValues(A, 1, &nnz, &i, cols, values);

        ClientContext setValuesContext;
        ::rpc_hypre::RPC_HYPRE_IJMatrixSetValuesMessage* setValuesRequest =
            new ::rpc_hypre::RPC_HYPRE_IJMatrixSetValuesMessage;

        setValuesRequest->set_allocated_matrix(&matrix);
        setValuesRequest->set_nrows(1);
        setValuesRequest->add_ncols(nnz);
        setValuesRequest->add_rows(i);

        for (int m = 0; m < nnz; ++m) {

          setValuesRequest->add_cols(cols[m]);
          setValuesRequest->add_values(values[m]);

        }

        stub_->RPC_HYPRE_IJMatrixSetValues(&setValuesContext, *setValuesRequest, &matrix);

      }
    }

    /* Assemble after setting the coefficients */

    /// ----------------------------------- ASSEMBLE -----------------------------------
//    HYPRE_IJMatrixAssemble(A);

    ClientContext assembleContext;
    stub_->RPC_HYPRE_IJMatrixAssemble(&assembleContext, matrix, &matrix);

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
    /// ----------------------------------- GET OBJECT -----------------------------------
//    HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);

    ClientContext getObjectContext;
    ::rpc_hypre::RPC_HYPRE_GenericMatrixIntegerParam* getObjectRequest =
        new ::rpc_hypre::RPC_HYPRE_GenericMatrixIntegerParam;
    stub_->RPC_HYPRE_IJMatrixGetObject(&getObjectContext, *getObjectRequest, &matrix);
    matrixIdentifier = matrix.identifier();

    /* Create the rhs and solution */

    /// ----------------------------------- VECTOR -----------------------------------

//    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper,&b);
//    HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
//    HYPRE_IJVectorInitialize(b);
//
//    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper,&x);
//    HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
//    HYPRE_IJVectorInitialize(x);
//
//    /* Set the rhs values to h^2 and the solution to zero */
//    {
//      double *rhs_values, *x_values;
//      int    *rows;
//
//      rhs_values =  (double*) calloc(local_size, sizeof(double));
//      x_values =  (double*) calloc(local_size, sizeof(double));
//      rows = (int*) calloc(local_size, sizeof(int));
//
//      for (i=0; i<local_size; i++)
//      {
//        rhs_values[i] = h2;
//        x_values[i] = 0.0;
//        rows[i] = ilower + i;
//      }
//
//      HYPRE_IJVectorSetValues(b, local_size, rows, rhs_values);
//      HYPRE_IJVectorSetValues(x, local_size, rows, x_values);
//
//      free(x_values);
//      free(rhs_values);
//      free(rows);
//    }
//
//
//    HYPRE_IJVectorAssemble(b);
//    /*  As with the matrix, for testing purposes, one may wish to read in a rhs:
//        HYPRE_IJVectorRead( <filename>, MPI_COMM_WORLD,
//                                  HYPRE_PARCSR, &b );
//        as an alternative to the
//        following sequence of HYPRE_IJVectors calls:
//        Create, SetObjectType, Initialize, SetValues, and Assemble
//    */
//    HYPRE_IJVectorGetObject(b, (void **) &par_b);
//
//    HYPRE_IJVectorAssemble(x);
//    HYPRE_IJVectorGetObject(x, (void **) &par_x);

    std::cout << "Finished matrix()" << std::endl;

    return Status::OK;

  }

  /**
   * execute ex5
   *
   * @param user
   * @return
   */
  int executeExample5() {

    if (!checkStatus(createSolver())) return 0;

    if (!checkStatus(setupMatrix())) return 0;

    if (!checkStatus(setupSolver())) return 0;

    if (!checkStatus(solve())) return 0;

    return 1;
  }

  int checkStatus(Status status) {
    if (!status.ok()) {
      std::cout << "Call was not successful. Error code: " << status.error_code() << ": " << status.error_message()
                << std::endl;
      return 0;
    }
    return 1;
  }

private:

//  hypreSrv::Stub* stub_;
  std::unique_ptr<hypreSrv::Stub> stub_;

};

int main(int argc, char** argv) {

  std::cout << "Starting Client." << std::endl;

  HypreClient client(grpc::CreateChannel("localhost:50051", grpc::InsecureChannelCredentials()));

  std::cout << "Started Client." << std::endl;

  std::cout << "Starting Example 5 execution." << std::endl;
  int success = client.executeExample5();
  std::cout << "Finished Example 5 execution." << std::endl;

  if (success) {
    std::cout << "Success!" << std::endl;
  } else {
    std::cout << "Failed!" << std::endl;
  }

  std::cout << "Closing Client." << std::endl;

  return 0;
}
