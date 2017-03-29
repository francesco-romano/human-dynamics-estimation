
#include <thrifts/XsensSegmentsFrame.h>

#include <human-ik/InverseKinematics.h>
#include <yarp/os/Network.h>
#include <yarp/os/BufferedPort.h>
#include <iDynTree/ModelIO/ModelLoader.h>
#include <iDynTree/Model/Model.h>
#include <iDynTree/KinDynComputations.h>
#include <yarp/os/ResourceFinder.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <iDynTree/Core/TestUtils.h>
#include <Eigen/Core>

int main(int argc, char ** argv) {

    using namespace yarp::os;
    using namespace iDynTree;

    Network net;

    ResourceFinder rf = ResourceFinder::getResourceFinderSingleton();
    rf.configure(argc, argv);

    human::InverseKinematics solver;
    solver.setVerbosityLevel(6);

    std::string parentFrame = "RightUpperLeg", childFrame = "Pelvis";
    ModelLoader loader;
    loader.loadModelFromFile(rf.findFile("urdf"));
    if (!loader.isValid()) {
        return -1;
    }

    KinDynComputations kinDyn;
    if (!solver.setModel(loader.model(), parentFrame, childFrame)) {
        return -1;
    }
    //Now configure the kinDynComputation objects
    std::vector<std::string> consideredJoints;
    solver.getConsideredJoints(consideredJoints);
    loader.loadReducedModelFromFile(rf.findFile("urdf"), consideredJoints);
    kinDyn.loadRobotModel(loader.model());

    FrameIndex parentFrameModelIndex = loader.model().getFrameIndex(parentFrame);
    FrameIndex childFrameModelIndex = loader.model().getFrameIndex(childFrame);

    size_t xsensParentIndex = 15;
    size_t xsensChildIndex = 0;

#if 0
    BufferedPort<xsens::XsensSegmentsFrame> port;

    if (!port.open("/human-state-provider-test/state:i")) {
        return -1;
    }

    if (!Network::connect("/xsens/frames:o", port.getName())) {
        return -1;
    }


    while (true) {

        xsens::XsensSegmentsFrame* frames = port.read();
        if (!frames) break;

        //read frames into transform
        xsens::XsensSegmentData rawFrame = frames->segmentsData[xsensParentIndex];

        Position position;
        position(0) = rawFrame.position.x;
        position(1) = rawFrame.position.y;
        position(2) = rawFrame.position.z;

        Vector4 quaternion;
        quaternion(0) = rawFrame.orientation.w;
        quaternion(1) = rawFrame.orientation.imaginary.x;
        quaternion(2) = rawFrame.orientation.imaginary.y;
        quaternion(3) = rawFrame.orientation.imaginary.z;

        Transform tParent(Rotation::RotationFromQuaternion(quaternion), position);

        rawFrame = frames->segmentsData[xsensChildIndex];
        position(0) = rawFrame.position.x;
        position(1) = rawFrame.position.y;
        position(2) = rawFrame.position.z;
        quaternion(0) = rawFrame.orientation.w;
        quaternion(1) = rawFrame.orientation.imaginary.x;
        quaternion(2) = rawFrame.orientation.imaginary.y;
        quaternion(3) = rawFrame.orientation.imaginary.z;

        Transform tChild(Rotation::RotationFromQuaternion(quaternion), position);

        Transform relative = tParent.inverse() * tChild;

        solver.setDesiredParentFrameAndEndEffectorTransformations(tParent, tChild);
        int resultCode = solver.runIK();
        VectorDynSize result = solver.getLastSolution();

        kinDyn.setJointPos(result);
        Transform ikRelative = kinDyn.getRelativeTransform(parentFrameModelIndex, childFrameModelIndex);

        Transform error = ikRelative.inverse() * relative;

        std::cerr << "XSENS\n" << toEigen(relative.asHomogeneousTransform()) << std::endl << std::endl;
        std::cerr << "\t\t\n" << toEigen(relative.getRotation().asRPY()).transpose() << std::endl << std::endl;
        std::cerr << "IK\n" << toEigen(ikRelative.asHomogeneousTransform()) << std::endl << std::endl;
        std::cerr << "ERR\n" << toEigen(error.asHomogeneousTransform()) << std::endl << std::endl;


        ASSERT_EQUAL_MATRIX_TOL(error.getRotation(), Rotation::Identity(), 1e-3);
        ASSERT_EQUAL_VECTOR_TOL(error.getPosition(), Position::Zero(), 1e-2);

    }
#else
    Eigen::Matrix<double, 4, 4, Eigen::RowMajor> errorRelativeTransform;
//    errorRelativeTransform << -0.0403048,  -0.299502,   0.953244, -0.0248676,
//                              -0.314952,    0.909194,   0.272345,  0.0745783,
//                              -0.948252,   -0.289249,  -0.130974, -0.0239621,
//                                      0,           0,          0,          1;

    errorRelativeTransform << 0.763888,  -0.151962,   0.627202, -0.0123844,
    0.0324943,   0.979705,   0.197792 , 0.0810412,
    -0.64453 , -0.130711 ,  0.753323, -0.0104554,
    0         , 0,          0,          1;

    Rotation rot;
    toEigen(rot) = errorRelativeTransform.block<3, 3>(0, 0);
    Position pos;
    pos(0) = errorRelativeTransform.coeffRef(0, 3);
    pos(1) = errorRelativeTransform.coeffRef(1, 3);
    pos(2) = errorRelativeTransform.coeffRef(2, 3);

    Transform desiredTransform = Transform(rot, pos);
    std::cerr << desiredTransform.inverse().getRotation().asRPY().toString() << std::endl;

    Vector3 weights;
    weights.zero();
    weights(1) = 1;
    weights(2) = 0;
    solver.setWeights(weights);
    VectorDynSize guess(3);
    toEigen(guess) = toEigen(desiredTransform.getRotation().asRPY());
    solver.setGuess(guess);
    solver.setDesiredTransformation(desiredTransform);
    solver.setDesiredJointPositions(guess);

    int resultCode = solver.runIK();
    VectorDynSize result = solver.getLastSolution();

    kinDyn.setJointPos(result);
    Transform ikRelative = kinDyn.getRelativeTransform(parentFrameModelIndex, childFrameModelIndex);

    Transform error = ikRelative.inverse() * desiredTransform;

    std::cerr << "Solution\n" << toEigen(result).transpose() << std::endl ;
    std::cerr << "XSENS\n" << toEigen(desiredTransform.asHomogeneousTransform()) << std::endl << std::endl;
    std::cerr << "IK\n" << toEigen(ikRelative.asHomogeneousTransform()) << std::endl << std::endl;
    std::cerr << "ERR\n" << toEigen(error.asHomogeneousTransform()) << std::endl << std::endl;
    std::cerr << "Rotation in RPY\n" << (180.0 / M_PI * toEigen(error.getRotation().asRPY())) << std::endl;


    ASSERT_EQUAL_MATRIX_TOL(error.getRotation(), Rotation::Identity(), 1e-3);
    ASSERT_EQUAL_VECTOR_TOL(error.getPosition(), Position::Zero(), 1e-2);
#endif
    std::cerr << "Finished reading frames" << std::endl;





}
