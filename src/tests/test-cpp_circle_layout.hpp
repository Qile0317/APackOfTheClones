#include <Rcpp.h>
#include "CircleNode.h"
#include "CirclePacker.h"

// function for testing of two nodes' numeric parameters x, y, rad.
bool xyr_are_equal(CircleNode& c1, CircleNode c2, double thr = 5e-5) {
    return approx_equal(c1.x, c2.x, thr) &&
        approx_equal(c1.y, c2.y, thr) &&
        approx_equal(c1.rad, c2.rad, thr);
}

bool nodes_are_equal(CircleNode& c1, CircleNode& c2, double thr = 5e-5) {
    return xyr_are_equal(c1, c2, thr) && (c1.nxt == c2.nxt) && (c1.prv == c2.prv);
}

context("C++ | cpp_circle_layout") {
    test_that("init_boundary works") {
        std::vector<double> v = {1,2,3,4};
        CirclePacker nodes = CirclePacker(v);

        nodes.init_boundary();
        expect_true(nodes.data[0].nxt == 1);
        expect_true(nodes.data[1].nxt == 2);
        expect_true(nodes.data[2].nxt == 0);

        expect_true(nodes.data[0].prv == 2);
        expect_true(nodes.data[1].prv == 0);
        expect_true(nodes.data[2].prv == 1);

        expect_true(nodes.data[3].nxt == -1);
        expect_true(nodes.data[3].prv == -1);
    }

    test_that("insert_circle works") {
        std::vector<double> v = {1,2,3,4};
        CirclePacker nodes = CirclePacker(v);
        nodes.init_boundary();

        int c1 = 1, c2 = 2, cn = 3;
        nodes.insert_circle(c1, c2, cn);

        expect_true(nodes.data[c1].nxt == cn);
        expect_true(nodes.data[cn].prv == c1);
        expect_true(nodes.data[c2].prv == cn);
        expect_true(nodes.data[cn].nxt == c2);
    }

    test_that("fwd_dist works") {
        int c1 = 0, c2 = 1, c3 = 2;
        std::vector<double> v = {1,2,3};
        CirclePacker nodes = CirclePacker(v);
        nodes.init_boundary();

        expect_true(nodes.fwd_dist(c3, c1) == 1);
        expect_true(nodes.fwd_dist(c1, c2) == 1);
        expect_true(nodes.fwd_dist(c2, c1) == 2);
    }

    test_that("fwd_remove works") {
        int c1 = 0, c2 = 1, c5 = 4;
        std::vector<double> v = {1,2,3,4,5};
        CirclePacker nodes = CirclePacker(v);
        nodes.init_boundary(c1, c5);

        nodes.fwd_remove(c2, c5);
        expect_true(nodes.data[c2].nxt == c5);
        expect_true(nodes.data[c2].prv == c1);
        expect_true(nodes.data[c5].prv == c2);
        expect_true(nodes.data[c5].nxt == c1);
    }

    test_that("centre_dist works") {
        CircleNode c1 = CircleNode(3, 4, 1);
        expect_true(centre_dist(c1) == 5);

        CircleNode c2 = CircleNode(-12, -35, 1);
        expect_true(centre_dist(c2) == 37);
    }

    test_that("fit_tang_circle works") { // probably could benefit from more testcases
        CirclePacker nodes = CirclePacker({
            CircleNode(1, 2, 3), CircleNode(6, 7, 8), CircleNode(50, 90, 1)
        });
        nodes.init_boundary();

        expect_true(nodes.fit_tang_circle(0, 1, 2) == 2);
        expect_true(approx_equal(nodes.data[2].x, -2.47718));
        expect_true(approx_equal(nodes.data[2].y, 3.97718));
    }

    test_that("tang_circle_dist works") {
        CirclePacker nodes = CirclePacker({
            CircleNode(1, 2, 3), CircleNode(6, 7, 8), CircleNode(50, 90, 1)
        });
        nodes.init_boundary();

        expect_true(approx_equal(nodes.tang_circle_dist(0, 1, 2), 4.685548));
    }

    test_that("place_starting_three works") {
        CirclePacker nodes = CirclePacker({
            CircleNode(1, 2, 3), CircleNode(6, 7, 8), CircleNode(50, 90, 1)
        });
        nodes.init_boundary();

        nodes.place_starting_three();
        expect_true(xyr_are_equal(nodes.data[0], CircleNode(-4.984898, -1.466558, 3)));
        expect_true(xyr_are_equal(nodes.data[1], CircleNode(6.015102, 3.533442, 8)));
        expect_true(xyr_are_equal(nodes.data[2], CircleNode(-1.030204, -2.066885, 1)));
    }

    test_that("closest works") {
        CirclePacker nodes = CirclePacker({
            CircleNode(1, 2, 3), CircleNode(6, 7, 8), CircleNode(50, 90, 1), CircleNode(-1, -3, 3), CircleNode(-10, -69, 4)
        });
        nodes.init_boundary(0, 4);

        expect_true(nodes.closest(2) == 0);
        expect_true(nodes.closest(0) == nodes.closest(4));
        expect_true(nodes.closest(1) == nodes.closest(3));
    }

    test_that("closest_place works") {
        CirclePacker nodes = CirclePacker({
            CircleNode(1, 2, 3), CircleNode(6, 7, 8), CircleNode(50, 90, 1)
        });
        nodes.init_boundary();

        // maybe should also test if nodes.data got modified
        expect_true(nodes.closest_place(0, 2) == 1);
    }

    test_that("do_intersect works") {
        CirclePacker nodes = CirclePacker({
            CircleNode(1, 2, 3), CircleNode(6, 7, 8), CircleNode(50, 90, 1), CircleNode(-1, -3, 3), CircleNode(-10, -69, 4)
        });
        nodes.init_boundary(0, 4);

        expect_true(nodes.do_intersect(0, 1));
        expect_true(nodes.do_intersect(0, 3));
        expect_false(nodes.do_intersect(1, 4));
        expect_false(nodes.do_intersect(2, 3));
    }

    test_that("geod_dist works") {
        CirclePacker nodes = CirclePacker({
            CircleNode(1, 2, 3), CircleNode(6, 7, 8), CircleNode(50, 90, 1), CircleNode(-1, -3, 3), CircleNode(-10, -69, 4)
        });
        nodes.init_boundary(0, 4);

        expect_true(nodes.geod_dist(0, 1, 2) == 1);
        expect_true(nodes.geod_dist(2, 1, 0) == 2);
        expect_true(nodes.geod_dist(1, 4, 2) == 3);
        expect_true(nodes.geod_dist(3, 0, 4) == 4);
    }

    test_that("is_clear works") {
        std::pair<int, int> trial = std::make_pair(-1, -1);
        expect_true(is_clear(trial));
    }

    test_that("construct_obstruct_list works") {
        CirclePacker nodes = CirclePacker({
            CircleNode(1, 2, 3), CircleNode(6, 7, 8), CircleNode(50, 90, 1), CircleNode(-1, -3, 3), CircleNode(-10, -69, 4)
        });
        nodes.init_boundary(0, 4);

        std::vector<int> trial = nodes.construct_obstruct_list(0, 4, 0);
        std::vector<int> expected = {1, 3};
        expect_true(trial == expected);

        // NEED MORE TESTCASES!!!!!!!!!!!!!!!!!
        // in previous implementation this was leading to infinite loops due to circular references
    }

    test_that("overlap_check works for 5 connected nodes") {
        CirclePacker nodes = CirclePacker({
            CircleNode(1, 2, 3), CircleNode(6, 7, 8), CircleNode(50, 90, 1), CircleNode(-1, -3, 3), CircleNode(-10, -69, 4)
        });
        nodes.init_boundary(0, 4);
        std::pair<int ,int> trial, expected;

        trial = nodes.overlap_check(0, 1, 4);
        expected = std::make_pair(4, 1);
        expect_true(trial == expected);

        trial = nodes.overlap_check(2, 4, 3);
        expected = std::make_pair(2, 0);
        expect_true(trial == expected);

        trial = nodes.overlap_check(3, 1, 4);
        expected = std::make_pair(-1, -1);
        expect_true(trial == expected);
    }

    test_that("overlap_check works for 3 connected nodes") {
        CirclePacker nodes = CirclePacker({
            CircleNode(1, 2, 3), CircleNode(6, 7, 8), CircleNode(50, 90, 1)
        });
        nodes.init_boundary();
        std::pair<int ,int> trial, expected;

        trial = nodes.overlap_check(0, 1, 2);
        expected = std::make_pair(0, 2);
        expect_true(trial == expected);
    }

    test_that("is_degenerate_case works") {
        std::vector<double> v = {1};
        CirclePacker nodes = CirclePacker(v);
        expect_true(nodes.is_degenerate_case());

        v = {1,2};
        nodes = CirclePacker(v);
        expect_true(nodes.is_degenerate_case());

        v = {1,2,3};
        nodes = CirclePacker(v);
        expect_false(nodes.is_degenerate_case());
    }

    test_that("handle_degenerate_cases for n = 1 works") {
        Rcpp::NumericVector centroid = {0, 0};
        std::vector<double> v = {69.0};
        CirclePacker nodes = CirclePacker(v);
        Rcpp::List trial = nodes.handle_degenerate_cases(centroid, 0);

        Rcpp::NumericVector trial_x = trial[0], trial_y = trial[1];
        Rcpp::NumericVector trial_rad = trial[2], trial_centroid = trial[3];
        double trial_clRad = trial[4];

        // testing for equivalence to expected values for trial 1
        Rcpp::NumericVector x = {0}, y = {0}, rad = {69};

        expect_true(elements_are_equal(trial_x, x));
        expect_true(elements_are_equal(trial_y, y));
        expect_true(elements_are_equal(trial_rad, rad));
        expect_true(elements_are_equal(trial_centroid, centroid));
        expect_true(approx_equal(trial_clRad, 69));

        // generating test trial 2 values with new centroid and rad_scale_factor
        nodes = CirclePacker(v);
        centroid = {4, 5};
        trial = nodes.handle_degenerate_cases(centroid, 1.5);

        trial_x = trial[0], trial_y = trial[1], trial_rad = trial[2];
        trial_centroid = trial[3], trial_clRad = trial[4];

        // testing for equivalence to expected values for trial 1
        x = {4}, y = {5}, rad = {67.5};

        expect_true(elements_are_equal(trial_x, x));
        expect_true(elements_are_equal(trial_y, y));
        expect_true(elements_are_equal(trial_rad, rad));
        expect_true(elements_are_equal(trial_centroid, centroid));
        expect_true(approx_equal(trial_clRad, 69));
    }

    test_that("handle_degenerate_cases for n = 2 works") {
        // generating test trial 1 values
        std::vector<double> v = {69.0, 420.0};
        CirclePacker nodes = CirclePacker(v);
        Rcpp::NumericVector centroid = {0, 0};
        Rcpp::List trial = nodes.handle_degenerate_cases(centroid, 0);

        Rcpp::NumericVector trial_x = trial[0], trial_y = trial[1];
        Rcpp::NumericVector trial_rad = trial[2], trial_centroid = trial[3];
        double trial_clRad = trial[4];

        // testing for equivalence to expected values for trial 1
        Rcpp::NumericVector x = {-69, 420}, y = {0, 0}, rad = {69, 420};

        expect_true(elements_are_equal(trial_x, x));
        expect_true(elements_are_equal(trial_y, y));
        expect_true(elements_are_equal(trial_rad, rad));
        expect_true(elements_are_equal(trial_centroid, centroid));
        expect_true(approx_equal(trial_clRad, 244.5));

        // generating test trial 2 values with new centroid and rad_scale_factor
        v = {4, 3};
        nodes = CirclePacker(v);
        centroid = {1, 2};
        trial = nodes.handle_degenerate_cases(centroid, 0.8);

        trial_x = trial[0], trial_y = trial[1], trial_rad = trial[2];
        trial_centroid = trial[3], trial_clRad = trial[4];

        // testing for equivalence to expected values for trial 1
        x = {-3, 4}, y = {2, 2}, rad = {3.2, 2.2};

        expect_true(elements_are_equal(trial_x, x));
        expect_true(elements_are_equal(trial_y, y));
        expect_true(elements_are_equal(trial_rad, rad));
        expect_true(elements_are_equal(trial_centroid, centroid));
        expect_true(approx_equal(trial_clRad, 3.5));
    }

    test_that("process_into_clusterlist works") {
        // create trial data, 5 circles, centroid (1,1), rad_scale_factor 0.8
        CirclePacker nodes = CirclePacker({
            CircleNode(1,2,3),CircleNode(6,7,8),CircleNode(50,90,1),CircleNode(-1,-3,3),CircleNode(-10,-69,4)
        });
        Rcpp::NumericVector centroid = {1, 1};
        Rcpp::List trial = nodes.process_into_clusterlist(
            centroid, 0.5, false
        );
        Rcpp::NumericVector trial_x = trial[0], trial_y = trial[1];
        Rcpp::NumericVector trial_rad = trial[2], trial_centroid = trial[3];
        double trial_clRad = trial[4];

        // test for equivalence to expected data
        Rcpp::NumericVector x = {2, 7, 51, 0, -9};
        Rcpp::NumericVector y = {3, 8, 91, -2, -68};
        Rcpp::NumericVector rad = {2.5, 7.5, 0.5, 2.5, 3.5};

        expect_true(elements_are_equal(trial_x, x));
        expect_true(elements_are_equal(trial_y, y));
        expect_true(elements_are_equal(trial_rad, rad));
        expect_true(elements_are_equal(trial_centroid, centroid));
        expect_true(approx_equal(trial_clRad, 50.5));
    }

    test_that("fit_circle works") {
        std::vector<double> v = {3,3,2,1,1};
        CirclePacker nodes = CirclePacker(v);
        nodes.init_boundary();
        nodes.place_starting_three();
        int j = 3;
        int curr_circ = nodes.closest(2);

        j = nodes.fit_circle(curr_circ, j, false);
        expect_true(j == 4);

        CirclePacker expected_nodes = CirclePacker({
            CircleNode(-3, 1.3333333333, 3),
            CircleNode(3, 1.3333333333, 3),
            CircleNode(0, -2.66666666667, 2),
            CircleNode(-3, -2.66666666667, 1),
            CircleNode(0,0,1)
        });
        expected_nodes.init_boundary(0,3); // important that its 0,3 !

        for (int i = 0; i < 5; i++) {
            expect_true(nodes_are_equal(
                expected_nodes.data[i], nodes.data[i]
            ));
        }
    }
}
