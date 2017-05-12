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

//#include "../../../hypre/src/utilities/_hypre_utilities.h"
//#include "../../../hypre/src/krylov/HYPRE_krylov.h"
//#include "../../../hypre/src/HYPRE.h"
//#include "../../../hypre/src/parcsr_ls/HYPRE_parcsr_ls.h"
//#include "../../../hypre/src/examples/vis.c"

using grpc::Channel;
using grpc::ClientContext;
using grpc::Status;
using rpc_hypre::hypreSrv;

class HypreClient {

protected:

  int solverIdentifier;

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
    rpc_hypre::RPC_HYPRE_BoomerAMGSolveMessage request;
    rpc_hypre::RPC_HYPRE_ParVector response;

    std::cout << "Running RPC_HYPRE_BoomerAMGSolve()." << std::endl;
    Status status = stub_->RPC_HYPRE_BoomerAMGSolve(&context, request, &response);
    std::cout << "Finished running RPC_HYPRE_BoomerAMGSolve()." << std::endl;

    return status;

  }

  /**
   * execute ex5
   *
   * @param user
   * @return
   */
  int executeExample5() {

    Status status = createSolver();

    if (!status.ok()) {
      std::cout << "Call was not successful. Error code: " << status.error_code() << ": " << status.error_message()
                << std::endl;
      return 0;
    }

    Status status2 = setupSolver();

    if (!status2.ok()) {
      std::cout << "Call was not successful. Error code: " << status2.error_code() << ": " << status2.error_message()
                << std::endl;
      return 0;
    }

    Status status3 = solve();

    if (!status3.ok()) {
      std::cout << "Call was not successful. Error code: " << status3.error_code() << ": " << status3.error_message()
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
