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

  }

private:
  std::unique_ptr<hypreSrv::Stub> stub_;

};

int main(int argc, char** argv) {

  HypreClient client(grpc::CreateChannel(
      "localhost:50051", grpc::InsecureChannelCredentials()));

  std::string reply = client.executeExample5();
  std::cout << "Received: " << reply << std::endl;

  return 0;
}
