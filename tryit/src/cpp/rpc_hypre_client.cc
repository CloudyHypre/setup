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

public:

  /**
   * construct
   * @param channel
   */
  HypreClient(std::shared_ptr<Channel> channel) : stub_(hypreSrv::NewStub(channel)) {}

  /**
   * execute ex5
   *
   * @param user
   * @return
   */
  std::string executeExample5() {

    //create

    std::cout << "Starting Example 5 execution." << std::endl;

    ClientContext context;
    ::google::protobuf::Empty emptyRequest;
    rpc_hypre::RPC_HYPRE_Solver solver;

    std::cout << "Running RPC_HYPRE_BoomerAMGCreate()." << std::endl;
    Status status = stub_->RPC_HYPRE_BoomerAMGCreate(&context, emptyRequest, &solver);
    std::cout << "Finished running RPC_HYPRE_BoomerAMGCreate()." << std::endl;

    if (!status.ok()) {
      std::cout << status.error_code() << ": " << status.error_message()
                << std::endl;
      return 0;
    }

    //set it up

    rpc_hypre::RPC_HYPRE_BoomerAMGSolveMessage request;
    rpc_hypre::RPC_HYPRE_ParCSRMatrix matrixA;
    rpc_hypre::RPC_HYPRE_ParVector par_b;
    rpc_hypre::RPC_HYPRE_ParVector par_x;

    request.set_allocated_solver(&solver);
    request.set_allocated_parcsr_a(&matrixA);
    request.set_allocated_par_b(&par_b);
    request.set_allocated_par_x(&par_x);

    ::google::protobuf::Empty emptyResponse;

    ClientContext context2;
    std::cout << "Running RPC_HYPRE_BoomerAMGSetup()." << std::endl;
    stub_->RPC_HYPRE_BoomerAMGSetup(&context2, request, &emptyResponse);
    std::cout << "Finished running RPC_HYPRE_BoomerAMGSetup()." << std::endl;

    if (!status.ok()) {
      std::cout << status.error_code() << ": " << status.error_message()
                << std::endl;
      return 0;
    }

    //solve (later)

//    stub_->RPC_HYPRE_BoomerAMGSolve();

//    if (!status.ok()) {
//      std::cout << status.error_code() << ": " << status.error_message()
//                << std::endl;
//      return 0;
//    }

    std::cout << "Finished Example 5 execution." << std::endl;


  }

private:

  std::unique_ptr<hypreSrv::Stub> stub_;
//  hypreSrv::Stub* stub_;

};

int main(int argc, char** argv) {

  std::cout << "Starting Client." << std::endl;

  HypreClient client(grpc::CreateChannel(
      "localhost:50051", grpc::InsecureChannelCredentials()));

  std::cout << "Started Client." << std::endl;

  std::string reply = client.executeExample5();
  std::cout << "Received: " << reply << std::endl;

  std::cout << "Closing Client." << std::endl;

  return 0;
}
