// function for testing of two nodes' numeric parameters x, y, rad.
bool xyr_are_equal(Node& c1, Node c2, double thr = 5e-5) {
    return approx_equal(c1.x, c2.x, thr) && 
        approx_equal(c1.y, c2.y, thr) && 
        approx_equal(c1.rad, c2.rad, thr);
}

context("C++ cpp_circle_layout") {
    test_that("init_boundary works") {
        NodeVector nodes = NodeVector({1,2,3,4});
        
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
        NodeVector nodes = NodeVector({1,2,3,4});
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
        NodeVector nodes = NodeVector({1,2,3});
        nodes.init_boundary();
        
        expect_true(nodes.fwd_dist(c3, c1) == 1);
        expect_true(nodes.fwd_dist(c1, c2) == 1);
        expect_true(nodes.fwd_dist(c2, c1) == 2);
    }

    test_that("fwd_remove works") {
        int c1 = 0, c2 = 1, c5 = 4;
        NodeVector nodes = NodeVector({1,2,3,4,5});
        nodes.init_boundary(c1, c5);
    
        nodes.fwd_remove(c2, c5);
        expect_true(nodes.data[c2].nxt == c5);
        expect_true(nodes.data[c2].prv == c1);
        expect_true(nodes.data[c5].prv == c2);
        expect_true(nodes.data[c5].nxt == c1);
    }

    test_that("centre_dist works") {
        Node c1 = Node(3, 4, 1);
        expect_true(centre_dist(c1) == 5);
    
        Node c2 = Node(-12, -35, 1);
        expect_true(centre_dist(c2) == 37);
    }

    test_that("fit_tang_circle works") { // probably could benefit from more testcases
        NodeVector nodes = NodeVector({
            Node(1, 2, 3), Node(6, 7, 8), Node(50, 90, 1)
        });
        nodes.init_boundary();
        
        expect_true(nodes.fit_tang_circle(0, 1, 2) == 2);
        expect_true(approx_equal(nodes.data[2].x, -2.47718));
        expect_true(approx_equal(nodes.data[2].y, 3.97718));
    }

    test_that("tang_circle_dist works") {
        NodeVector nodes = NodeVector({
            Node(1, 2, 3), Node(6, 7, 8), Node(50, 90, 1)
        });
        nodes.init_boundary();
        
        expect_true(approx_equal(nodes.tang_circle_dist(0, 1, 2), 4.685548));
    }

    test_that("place_starting_three works") {
        NodeVector nodes = NodeVector({
            Node(1, 2, 3), Node(6, 7, 8), Node(50, 90, 1)
        });
        nodes.init_boundary();
        
        nodes.place_starting_three();
        expect_true(xyr_are_equal(nodes.data[0], Node(-4.984898, -1.466558, 3)));
        expect_true(xyr_are_equal(nodes.data[1], Node(6.015102, 3.533442, 8)));
        expect_true(xyr_are_equal(nodes.data[2], Node(-1.030204, -2.066885, 1)));
    }

    test_that("closest works") {
        NodeVector nodes = NodeVector({
            Node(1, 2, 3), Node(6, 7, 8), Node(50, 90, 1), Node(-1, -3, 3), Node(-10, -69, 4)
        });
        nodes.init_boundary(0, 4);

        expect_true(nodes.closest(2) == 0);
        expect_true(nodes.closest(0) == nodes.closest(4)); 
        expect_true(nodes.closest(1) == nodes.closest(3));
    }

    test_that("closest_place works") {
        NodeVector nodes = NodeVector({
            Node(1, 2, 3), Node(6, 7, 8), Node(50, 90, 1)
        });
        nodes.init_boundary();

        // maybe should also test if nodes.data got modified
        expect_true(nodes.closest_place(0, 2) == 1);
    }

    test_that("do_intersect works") {
        NodeVector nodes = NodeVector({
            Node(1, 2, 3), Node(6, 7, 8), Node(50, 90, 1), Node(-1, -3, 3), Node(-10, -69, 4)
        });
        nodes.init_boundary(0, 4);

        expect_true(nodes.do_intersect(0, 1));
        expect_true(nodes.do_intersect(0, 3));
        expect_false(nodes.do_intersect(1, 4));
        expect_false(nodes.do_intersect(2, 3));
    }

    test_that("geod_dist works") {
        NodeVector nodes = NodeVector({
            Node(1, 2, 3), Node(6, 7, 8), Node(50, 90, 1), Node(-1, -3, 3), Node(-10, -69, 4)
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
        NodeVector nodes = NodeVector({
            Node(1, 2, 3), Node(6, 7, 8), Node(50, 90, 1), Node(-1, -3, 3), Node(-10, -69, 4)
        });
        nodes.init_boundary(0, 4);

        std::vector<int> trial = nodes.construct_obstruct_list(0, 4, 0);
        std::vector<int> expected = {1, 3};
        expect_true(trial == expected);

        // NEED MORE TESTCASES!!!!!!!!!!!!!!!!!
        // in previous implementation this was leading to infinite loops due to circular references
    }

    test_that("overlap_check works for 5 connected nodes") {
        NodeVector nodes = NodeVector({
            Node(1, 2, 3), Node(6, 7, 8), Node(50, 90, 1), Node(-1, -3, 3), Node(-10, -69, 4)
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
        NodeVector nodes = NodeVector({
            Node(1, 2, 3), Node(6, 7, 8), Node(50, 90, 1)
        });
        nodes.init_boundary();
        std::pair<int ,int> trial, expected;

        trial = nodes.overlap_check(0, 1, 2);
        expected = std::make_pair(0, 2);
        expect_true(trial == expected);
    }
  
    test_that("is_degenerate_case works") {
        NodeVector nodes = NodeVector({1});
        expect_true(nodes.is_degenerate_case());

        nodes = NodeVector({1, 2});
        expect_true(nodes.is_degenerate_case());

        nodes = NodeVector({1, 2, 3});
        expect_false(nodes.is_degenerate_case());
    }

    // unfinished
}