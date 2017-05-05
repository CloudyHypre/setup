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
#include "vis.c"

#include "hypre.pb.h"
#include "hypre.grpc.pb.h"

#include "../../../hypre/src/utilities/_hypre_utilities.h"
#include "../../../hypre/src/krylov/HYPRE_krylov.h"
#include "../../../hypre/src/HYPRE.h"
#include "../../../hypre/src/parcsr_ls/HYPRE_parcsr_ls.h"
#include "../../../hypre/src/examples/vis.c"

using grpc::Channel;
using grpc::ClientContext;
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

    ClientContext context;
    ::google::protobuf::Empty emptyRequest;
    rpc_hypre::RPC_HYPRE_Solver solver;
    Status status = stub_->RPC_HYPRE_BoomerAMGCreate(&context, emptyRequest, &solver);

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
    request.set_allocated_solver(solver);
    request.set_allocated_parcsr_a(matrixA);
    request.set_allocated_par_b(par_b);
    request.set_allocated_par_x(par_x);
    ::google::protobuf::Empty emptyResponse;
    stub_->RPC_HYPRE_BoomerAMGSetup(&context, request, emptyResponse);

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

  }

private:

//  std::unique_ptr<hypreSrv::Stub> stub_;
  hypreSrv::Stub* stub_;

};

int main(int argc, char** argv) {

  HypreClient client(grpc::CreateChannel(
      "localhost:50051", grpc::InsecureChannelCredentials()));

  std::string reply = client.executeExample5();
  std::cout << "Received: " << reply << std::endl;

  return 0;
}
