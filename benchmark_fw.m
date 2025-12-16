n = 10
% eps = 1e-6 - already in generate_test_case

generate_test_case(n, "interior", "normal", 9)
generate_test_case(n, "interior", "away", 9)
generate_test_case(n, "box_boundary", "normal", 9)
generate_test_case(n, "box_boundary", "away", 9)
generate_test_case(n, "active_linear", "normal", 9)
generate_test_case(n, "active_linear", "away", 9)