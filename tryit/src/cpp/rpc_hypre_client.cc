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
  int vectorBIdentifier;
  int vectorXIdentifier;

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
    rpc_hypre::RPC_HYPRE_ParCSRMatrix* matrixA = new rpc_hypre::RPC_HYPRE_ParCSRMatrix;
    rpc_hypre::RPC_HYPRE_ParVector* par_b = new rpc_hypre::RPC_HYPRE_ParVector;
    rpc_hypre::RPC_HYPRE_ParVector* par_x = new rpc_hypre::RPC_HYPRE_ParVector;
    rpc_hypre::RPC_HYPRE_Solver* solverMessage = new rpc_hypre::RPC_HYPRE_Solver;
    solverMessage->set_identifier(solverIdentifier);
    matrixA->set_identifier(matrixIdentifier);

    rpc_hypre::RPC_HYPRE_BoomerAMGSolveMessage request;
    request.set_allocated_solver(solverMessage);
    request.set_allocated_parcsr_a(matrixA);
    request.set_allocated_par_b(par_b);
    request.set_allocated_par_x(par_x);

    Status status = stub_->RPC_HYPRE_BoomerAMGSetup((new ClientContext), request, (new ::rpc_hypre::Empty));

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
    int myid = 0;

    /* Default problem parameters */
    n = 33;

    /* Preliminaries: want at least one processor per row */
    if (n*n < num_procs) n = sqrt(num_procs) + 1;
    N = n*n; /* global number of rows */

    /* Each processor knows only of its own rows - the range is denoted by ilower
       and upper.  Here we partition the rows. We account for the fact that
       N may not divide evenly by the number of processors. */
    local_size = N/num_procs;
    extra = N - local_size*num_procs;

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
    if (!status.ok()) {
      return status;
    }

    {
      int nnz;
      double values[5];
      int cols[5];

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

        status = stub_->RPC_HYPRE_IJMatrixSetValues(&setValuesContext, *setValuesRequest, &matrix);
        if (!status.ok()) {
          return status;
        }

      }
    }

    /* Assemble after setting the coefficients */

    /// ----------------------------------- ASSEMBLE -----------------------------------
//    HYPRE_IJMatrixAssemble(A);

    ClientContext assembleContext;
    status = stub_->RPC_HYPRE_IJMatrixAssemble(&assembleContext, matrix, &matrix);
    if (!status.ok()) {
      return status;
    }

    /* Get the parcsr matrix object to use */
    /// ----------------------------------- GET OBJECT -----------------------------------
//    HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);

    ClientContext getObjectContext;
    ::rpc_hypre::RPC_HYPRE_GenericMatrixIntegerParam* getObjectRequest =
        new ::rpc_hypre::RPC_HYPRE_GenericMatrixIntegerParam;

    getObjectRequest->set_allocated_matrix(&matrix);
    status = stub_->RPC_HYPRE_IJMatrixGetObject(&getObjectContext, *getObjectRequest, &matrix);
    if (!status.ok()) {
      return status;
    }

    matrixIdentifier = matrix.identifier();

    std::cout << "Finished matrix()" << std::endl;

    return Status::OK;

  }

  Status setupVectors() {

    int i;
    int N, n;

    int ilower, iupper;
    int local_size, extra;

    int num_procs = 1;
    int myid = 0;

    double h, h2;

    Status status;
    ::rpc_hypre::RPC_HYPRE_IJVector vectorB;
    ::rpc_hypre::RPC_HYPRE_IJVector vectorX;

    /* Default problem parameters */
    n = 33;

    /* Preliminaries: want at least one processor per row */
    if (n*n < num_procs) n = sqrt(num_procs) + 1;
    N = n*n; /* global number of rows */
    h = 1.0/(n+1); /* mesh size*/
    h2 = h*h;


    /* Each processor knows only of its own rows - the range is denoted by ilower
       and upper.  Here we partition the rows. We account for the fact that
       N may not divide evenly by the number of processors. */
    local_size = N/num_procs;
    extra = N - local_size*num_procs;

    ilower = local_size*myid;
    ilower += hypre_min(myid, extra);

    iupper = local_size*(myid+1);
    iupper += hypre_min(myid+1, extra);
    iupper = iupper - 1;

    /* How many rows do I have? */
    local_size = iupper - ilower + 1;

    /* Create the rhs and solution */
    /// ----------------------------------- CREATE -----------------------------------
//    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper,&b);

    ::rpc_hypre::RPC_HYPRE_IJVectorCreateMessage createBRequest;
    createBRequest.set_jlower(ilower);
    createBRequest.set_jupper(iupper);

    status = stub_->RPC_HYPRE_IJVectorCreate((new ClientContext), createBRequest, &vectorB);
    if (!status.ok()) {
      return status;
    }

    /// ----------------------------------- SET OBJECT -----------------------------------
//    HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);

    ::rpc_hypre::RPC_HYPRE_GenericVectorxIntegerParam setBTypeRequest;
    setBTypeRequest.set_value(HYPRE_PARCSR);
    setBTypeRequest.set_allocated_vector(&vectorB);
    status = stub_->RPC_HYPRE_IJVectorSetObjectType((new ClientContext), setBTypeRequest, &vectorB);
    if (!status.ok()) {
      return status;
    }

    //we want to reuse and avoid double free of vector
    setBTypeRequest.release_vector();

        /// ----------------------------------- INIT -----------------------------------
//    HYPRE_IJVectorInitialize(b);

    status = stub_->RPC_HYPRE_IJVectorInitialize((new ClientContext), vectorB, &vectorB);
    if (!status.ok()) {
      return status;
    }

    /// ----------------------------------- CREATE -----------------------------------
//    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper,&x);
    ::rpc_hypre::RPC_HYPRE_IJVectorCreateMessage createXRequest;
    createXRequest.set_jlower(ilower);
    createXRequest.set_jupper(iupper);

    status = stub_->RPC_HYPRE_IJVectorCreate((new ClientContext), createXRequest, &vectorX);
    if (!status.ok()) {
      return status;
    }

    /// ----------------------------------- SET OBJECT -----------------------------------
//    HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
    ::rpc_hypre::RPC_HYPRE_GenericVectorxIntegerParam setXTypeRequest;
    setXTypeRequest.set_value(HYPRE_PARCSR);
    setXTypeRequest.set_allocated_vector(&vectorX);

    status = stub_->RPC_HYPRE_IJVectorSetObjectType((new ClientContext), setXTypeRequest, &vectorX);
    //we want to reuse and avoid double free of vector
    setXTypeRequest.release_vector();

    if (!status.ok()) {
      return status;
    }


    /// ----------------------------------- INIT -----------------------------------
//    HYPRE_IJVectorInitialize(x);
    status = stub_->RPC_HYPRE_IJVectorInitialize((new ClientContext), vectorX, &vectorX);
    if (!status.ok()) {
      return status;
    }

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

      /// ----------------------------------- SET VAL -----------------------------------
//      HYPRE_IJVectorSetValues(b, local_size, rows, rhs_values);

      ::rpc_hypre::RPC_HYPRE_IJVectorSetValuesMessage setBValuesRequest;
      setBValuesRequest.set_allocated_vector(&vectorB);
      setBValuesRequest.set_nvalues(local_size);
      for (i=0; i<local_size; i++)
      {
        setBValuesRequest.add_indices(rows[i]);
        setBValuesRequest.add_values(rhs_values[i]);
      }

      status = stub_->RPC_HYPRE_IJVectorSetValues((new ClientContext), setBValuesRequest, &vectorB);
      setBValuesRequest.release_vector();
      if (!status.ok()) {
        return status;
      }


      /// ----------------------------------- SET VAL -----------------------------------
//      HYPRE_IJVectorSetValues(x, local_size, rows, x_values);

      ::rpc_hypre::RPC_HYPRE_IJVectorSetValuesMessage setXValuesRequest;
      setXValuesRequest.set_allocated_vector(&vectorX);
      setXValuesRequest.set_nvalues(local_size);
      for (i=0; i<local_size; i++)
      {
        setXValuesRequest.add_indices(rows[i]);
        setXValuesRequest.add_values(x_values[i]);
      }

      status = stub_->RPC_HYPRE_IJVectorSetValues((new ClientContext), setXValuesRequest, &vectorX);
      setXValuesRequest.release_vector();
      if (!status.ok()) {
        return status;
      }

      free(x_values);
      free(rhs_values);
      free(rows);
    }


    /// ----------------------------------- ASSEMBLE -----------------------------------
//    HYPRE_IJVectorAssemble(b);
    status = stub_->RPC_HYPRE_IJVectorAssemble((new ClientContext), vectorB, &vectorB);
    if (!status.ok()) {
      return status;
    }

    /// ----------------------------------- GET OBJ -----------------------------------
//    HYPRE_IJVectorGetObject(b, (void **) &par_b);

    ::rpc_hypre::RPC_HYPRE_GenericVectorxIntegerParam getBObjRequest;
    getBObjRequest.set_allocated_vector(&vectorB);

    status = stub_->RPC_HYPRE_IJVectorGetObject((new ClientContext), getBObjRequest, &vectorB);
    getBObjRequest.release_vector();
    if (!status.ok()) {
      return status;
    }

    /// ----------------------------------- ASSEMBLE -----------------------------------
//    HYPRE_IJVectorAssemble(x);
    status = stub_->RPC_HYPRE_IJVectorAssemble((new ClientContext), vectorX, &vectorX);
    if (!status.ok()) {
      return status;
    }

    /// ----------------------------------- GET OBJ -----------------------------------
//    HYPRE_IJVectorGetObject(x, (void **) &par_x);

    ::rpc_hypre::RPC_HYPRE_GenericVectorxIntegerParam getXObjRequest;
    getXObjRequest.set_allocated_vector(&vectorX);

    status = stub_->RPC_HYPRE_IJVectorGetObject((new ClientContext), getXObjRequest, &vectorX);
    getXObjRequest.release_vector();
    if (!status.ok()) {
      return status;
    }

    vectorBIdentifier = vectorB.identifier();
    vectorXIdentifier = vectorX.identifier();

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
    if (!checkStatus(setupVectors())) return 0;
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
