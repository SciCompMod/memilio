#include "epidemiology/secir/contact_matrix.h"
#include "matchers.h"
#include "gtest/gtest.h"

TEST(TestContactMatrix, initZero)
{
    epi::ContactMatrix cm(Eigen::Index(3));
    EXPECT_EQ(print_wrap(cm.get_matrix_at(-1e5)), print_wrap(Eigen::MatrixXd::Zero(3, 3)));
    EXPECT_EQ(print_wrap(cm.get_matrix_at(0)), print_wrap(Eigen::MatrixXd::Zero(3, 3)));
    EXPECT_EQ(print_wrap(cm.get_matrix_at(1e-32)), print_wrap(Eigen::MatrixXd::Zero(3, 3)));
    EXPECT_EQ(print_wrap(cm.get_matrix_at(1e5)), print_wrap(Eigen::MatrixXd::Zero(3, 3)));
}

TEST(TestContactMatrix, initBaseAndMin)
{
    auto B = (Eigen::MatrixXd(2, 2) << 1, 2, 3, 4).finished();
    auto M = Eigen::MatrixXd::Constant(2, 2, 0.1);
    epi::ContactMatrix cm(B, M);
    EXPECT_EQ(print_wrap(cm.get_matrix_at(-1e5)), print_wrap(B));
    EXPECT_EQ(print_wrap(cm.get_matrix_at(0)), print_wrap(B));
    EXPECT_EQ(print_wrap(cm.get_matrix_at(1e-32)), print_wrap(B));
    EXPECT_EQ(print_wrap(cm.get_matrix_at(1e5)), print_wrap(B));
}

TEST(TestContactMatrix, dampings)
{
    auto B = (Eigen::MatrixXd(2, 2) << 1, 2, 3, 4).finished();
    auto M = Eigen::MatrixXd::Constant(2, 2, 0.1);
    epi::ContactMatrix cm(B, M);
    auto D1 = 0.25;
    auto D2 = (Eigen::MatrixXd(2, 2) << 0.0, 0.75, 0.5, 0.25).finished();
    cm.add_damping(D1, epi::DampingLevel(7), epi::DampingType(3), epi::SimulationTime(0.5));
    cm.add_damping(D2, epi::DampingLevel(7), epi::DampingType(2), epi::SimulationTime(2.0));

    EXPECT_EQ(cm.get_dampings().size(), 2);
    EXPECT_THAT(
        cm.get_dampings(),
        testing::ElementsAre(epi::Damping(2, D1, epi::DampingLevel(7), epi::DampingType(3), epi::SimulationTime(0.5)),
                             epi::Damping(D2, epi::DampingLevel(7), epi::DampingType(2), epi::SimulationTime(2.0))));

    EXPECT_EQ(print_wrap(cm.get_matrix_at(-1e5)), print_wrap(B));
    EXPECT_EQ(print_wrap(cm.get_matrix_at(-0.5)), print_wrap(B));
    EXPECT_THAT(print_wrap(cm.get_matrix_at(0.5 + 1e-32)), MatrixNear(B - D1 * (B - M)));
    EXPECT_THAT(print_wrap(cm.get_matrix_at(1e5)), MatrixNear(B - ((D1 + D2.array()) * (B - M).array()).matrix()));
}

TEST(TestContactMatrixGroup, sum)
{
    epi::ContactMatrixGroup cmg(2, 3);
    cmg[0] = epi::ContactMatrix(Eigen::MatrixXd::Constant(3, 3, 1.0));
    cmg[1] = epi::ContactMatrix(Eigen::MatrixXd::Constant(3, 3, 2.0));
    cmg[2] = epi::ContactMatrix(Eigen::MatrixXd::Constant(3, 3, 3.0));
    cmg.add_damping(0.5, epi::DampingLevel(3), epi::DampingType(1), epi::SimulationTime(1.0));

    EXPECT_THAT(print_wrap(cmg.get_matrix_at(0.0)), MatrixNear(Eigen::MatrixXd::Constant(3, 3, 6.0)));
    EXPECT_THAT(print_wrap(cmg.get_matrix_at(1.0)), MatrixNear(Eigen::MatrixXd::Constant(3, 3, 3.0)));
}