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
#include "../../../hypre/src/IJ_mv/_hypre_IJ_mv.h"
#include "../../../hypre/src/parcsr_mv/HYPRE_parcsr_mv.h"

using grpc::Status;
using grpc::ServerBuilder;
using grpc::Server;

using namespace ::rpc_hypre;

/**
 * think about this approach to store all parameters received
 */
class HypreContext {

  //some security related stuff could be stored here, too

  ::rpc_hypre::RPC_HYPRE_Solver solver;
  HYPRE_Solver hypreSolver;

  ::rpc_hypre::RPC_HYPRE_IJMatrix matrixMessageObject;
  HYPRE_IJMatrix hypreIJMatrix;
  HYPRE_ParCSRMatrix parcsr_A;

  ::rpc_hypre::RPC_HYPRE_IJVector vectorMessageObject;
  HYPRE_IJVector hypreIJVector;
  HYPRE_ParVector par_b;
  HYPRE_ParVector par_x;

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

  std::vector<::rpc_hypre::RPC_HYPRE_IJVector> vectorList;
  std::vector<HYPRE_IJVector> hypreVectorList;
  std::vector<::rpc_hypre::RPC_HYPRE_ParVector> parVectorList;
  std::vector<HYPRE_ParVector> hypreParVectorList;

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

  int getVectorIdentifier() {

    return hypreVectorList.size()-1;

  }

  int getParVectorIdentifier() {

    return hypreParVectorList.size()-1;

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

  /// ----------------------------------- MATRIX -----------------------------------

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

    std::cout << "Create Matrix..." << std::endl;

    HYPRE_IJMatrix hypreMatrix;
    HYPRE_IJMatrixCreate(
        MPI_COMM_WORLD, request->ilower(), request->iupper(), request->jlower(), request->jupper(), &hypreMatrix
    );

    //save the matrix
    matrixList.push_back(*matrix);
    hypreMatrixList.push_back(hypreMatrix);

    int id = getMatrixIdentifier();

    std::cout << "Created Matrix:" << id << std::endl;

    //get the id
    matrix->set_identifier(id);
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

    std::cout << "Matrix:" << request->matrix().identifier()
              << " set object type"
              << std::endl;

    HYPRE_IJMatrix hypreMatrix = hypreMatrixList[request->matrix().identifier()];
    HYPRE_IJMatrixSetObjectType(hypreMatrix, request->value());

    response->set_identifier(request->matrix().identifier());
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

    std::cout << "Matrix:" << request->identifier()
              << " init"
              << std::endl;

    HYPRE_IJMatrix hypreMatrix = hypreMatrixList[request->identifier()];
    HYPRE_IJMatrixInitialize(hypreMatrix);

    response->set_identifier(request->identifier());
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

    response->set_identifier(request->matrix().identifier());
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

    std::cout << "Matrix:" << request->identifier()
              << " assemble"
              << std::endl;

    HYPRE_IJMatrix hypreMatrix = hypreMatrixList[request->identifier()];
    HYPRE_IJMatrixAssemble(hypreMatrix);

    response->set_identifier(request->identifier());
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

    std::cout << "Matrix:" << request->matrix().identifier()
              << " get object"
              << std::endl;

    HYPRE_IJMatrix hypreMatrix = hypreMatrixList[request->matrix().identifier()];
    RPC_HYPRE_ParCSRMatrix parcsr_A;
    HYPRE_ParCSRMatrix hypreParcsr_A;

    HYPRE_IJMatrixGetObject(hypreMatrix, (void**) &hypreParcsr_A);
    parMatrixList.push_back(parcsr_A);
    hypreParMatrixList.push_back(hypreParcsr_A);

    int id = getParMatrixIdentifier();
    parcsr_A.set_identifier(id);
    response->set_identifier(id);
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

    std::cout << "Matrix:" << request->identifier()
              << " destroy"
              << std::endl;

    HYPRE_IJMatrix hypreMatrix = hypreMatrixList[request->identifier()];
    HYPRE_IJMatrixDestroy(hypreMatrix);

    return Status::OK;

  }

  /// ----------------------------------- VECTOR -----------------------------------

  /**
   * @param context
   * @param request
   * @param vector
   * @return Status
   */
  Status RPC_HYPRE_IJVectorCreate(::grpc::ServerContext* context,
                                  const ::rpc_hypre::RPC_HYPRE_IJVectorCreateMessage* request,
                                  ::rpc_hypre::RPC_HYPRE_IJVector* vector) override {

    HYPRE_IJVector hypreVector;
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, request->jlower(), request->jupper(), &hypreVector);

    vectorList.push_back(*vector);
    hypreVectorList.push_back(hypreVector);

    int id = getVectorIdentifier();
    std::cout << "Created Vector: " << id << std::endl;

    vector->set_identifier(id);
    return Status::OK;

  }

  /**
   * @param context
   * @param request
   * @param response
   * @return Status
   */
  Status RPC_HYPRE_IJVectorSetObjectType(::grpc::ServerContext* context,
                                         const ::rpc_hypre::RPC_HYPRE_GenericVectorxIntegerParam* request,
                                         ::rpc_hypre::RPC_HYPRE_IJVector* response) override {

    HYPRE_IJVector hypreIjVector = hypreVectorList[request->vector().identifier()];
    HYPRE_IJVectorSetObjectType(hypreIjVector, request->value());

    response->set_identifier(request->vector().identifier());
    return Status::OK;

  }

  /**
   * @param context
   * @param request
   * @param response
   * @return Status
   */
  Status RPC_HYPRE_IJVectorInitialize(::grpc::ServerContext* context,
                                      const ::rpc_hypre::RPC_HYPRE_IJVector* request,
                                      ::rpc_hypre::RPC_HYPRE_IJVector* response) override {

    HYPRE_IJVector hypreIjVector = hypreVectorList[request->identifier()];
    HYPRE_IJVectorInitialize(hypreIjVector);

    response->set_identifier(request->identifier());
    return Status::OK;

  }

  /**
   * @param context
   * @param request
   * @param response
   * @return Status
   */
  Status RPC_HYPRE_IJVectorSetValues(::grpc::ServerContext* context,
                                     const ::rpc_hypre::RPC_HYPRE_IJVectorSetValuesMessage* request,
                                     ::rpc_hypre::RPC_HYPRE_IJVector* response) override {

    HYPRE_IJVector hypreIjVector = hypreVectorList[request->vector().identifier()];
    HYPRE_IJVectorSetValues(
        hypreIjVector,
        request->nvalues(),
        repeatedToArray(request->indices_size(), request->indices()),
        repeatedToArray(request->values_size(), request->values())
    );

    response->set_identifier(request->vector().identifier());
    return Status::OK;

  }

  /**
   * @param context
   * @param request
   * @param response
   * @return Status
   */
  Status RPC_HYPRE_IJVectorAssemble(::grpc::ServerContext* context,
                                    const ::rpc_hypre::RPC_HYPRE_IJVector* request,
                                    ::rpc_hypre::RPC_HYPRE_IJVector* response) override {

    HYPRE_IJVector hypreIjVector = hypreVectorList[request->identifier()];
    HYPRE_IJVectorAssemble(hypreIjVector);

    response->set_identifier(request->identifier());
    return Status::OK;

  }

  /**
   * @param context
   * @param request
   * @param response
   * @return Status
   */
  Status RPC_HYPRE_IJVectorGetObject(::grpc::ServerContext* context,
                                     const ::rpc_hypre::RPC_HYPRE_GenericVectorxIntegerParam* request,
                                     ::rpc_hypre::RPC_HYPRE_IJVector* response) override {

    HYPRE_IJVector hypreIjVector = hypreVectorList[request->vector().identifier()];
    RPC_HYPRE_ParVector parVector;
    HYPRE_ParVector hypreParVector;

    HYPRE_IJVectorGetObject(hypreIjVector, (void**) &hypreParVector);
    parVectorList.push_back(parVector);
    hypreParVectorList.push_back(hypreParVector);

    int id = getParVectorIdentifier();
    std::cout << "par vector created: " << id << std::endl;

    parVector.set_identifier(id);
    response->set_identifier(id);
    return Status::OK;

  }

  /**
   * @param context
   * @param request
   * @param response
   * @return Status
   */
  Status RPC_HYPRE_IJVectorDestroy(::grpc::ServerContext* context,
                                   const ::rpc_hypre::RPC_HYPRE_IJVector* request,
                                   ::rpc_hypre::Empty* response) override {

    HYPRE_IJVector hypreIjVector = hypreVectorList[request->identifier()];
    HYPRE_IJVectorDestroy(hypreIjVector);

    return Status::OK;

  }


  /// ----------------------------------- BOOMER -----------------------------------

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

    ///solver
    ::rpc_hypre::RPC_HYPRE_Solver solver = request->solver();
    //there are no identifiers < 0
    if ((unsigned)solver.identifier() >= hypreSolversList.size()) {
      return Status::CANCELLED;
    }
    HYPRE_Solver hypreSolver = hypreSolversList[solver.identifier()];

    ///matrix
    ::rpc_hypre::RPC_HYPRE_ParCSRMatrix parMatrix = request->parcsr_a();
    if ((unsigned)parMatrix.identifier() >= hypreParMatrixList.size()) {
      return Status::CANCELLED;
    }
    HYPRE_ParCSRMatrix hypreParMatrix = hypreParMatrixList[parMatrix.identifier()];

    ///vector b
    ::rpc_hypre::RPC_HYPRE_ParVector parVectorB = request->par_b();
    if ((unsigned)parVectorB.identifier() >= hypreParVectorList.size()) {
      return Status::CANCELLED;
    }
    HYPRE_ParVector hypreParVectorB = hypreParVectorList[parVectorB.identifier()];

    ///vector x
    ::rpc_hypre::RPC_HYPRE_ParVector parVectorX = request->par_x();
    if ((unsigned)parVectorX.identifier() >= hypreParVectorList.size()) {
      return Status::CANCELLED;
    }
    HYPRE_ParVector hypreParVectorX = hypreParVectorList[parVectorX.identifier()];

    std::cout << "setup" << std::endl;
    HYPRE_BoomerAMGSetup(hypreSolver, hypreParMatrix, hypreParVectorB, hypreParVectorX);
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

    ///solver
    ::rpc_hypre::RPC_HYPRE_Solver solver = request->solver();
    //there are no identifiers < 0
    if ((unsigned)solver.identifier() >= hypreSolversList.size()) {
      return Status::CANCELLED;
    }
    HYPRE_Solver hypreSolver = hypreSolversList[solver.identifier()];

    ///matrix
    ::rpc_hypre::RPC_HYPRE_ParCSRMatrix parMatrix = request->parcsr_a();
    if ((unsigned)parMatrix.identifier() >= hypreParMatrixList.size()) {
      return Status::CANCELLED;
    }
    HYPRE_ParCSRMatrix hypreParMatrix = hypreParMatrixList[parMatrix.identifier()];

    ///vector b
    ::rpc_hypre::RPC_HYPRE_ParVector parVectorB = request->par_b();
    if ((unsigned)parVectorB.identifier() >= hypreParVectorList.size()) {
      return Status::CANCELLED;
    }
    HYPRE_ParVector hypreParVectorB = hypreParVectorList[parVectorB.identifier()];

    ///vector x
    ::rpc_hypre::RPC_HYPRE_ParVector parVectorX = request->par_x();
    if ((unsigned)parVectorX.identifier() >= hypreParVectorList.size()) {
      return Status::CANCELLED;
    }
    HYPRE_ParVector hypreParVectorX = hypreParVectorList[parVectorX.identifier()];

    std::cout << "solve" << std::endl;
    HYPRE_BoomerAMGSolve(hypreSolver, hypreParMatrix, hypreParVectorB, hypreParVectorX);

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
