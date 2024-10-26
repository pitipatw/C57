using ResourceConstrainedShortestPaths

adjlists = [0 1 1  ; 1 0 1 ;  1 1 0];
times = [0 1 2];
ub  = 10 
lb = 1 
start_times = [0, 1, 2] ; 
end_times = [1, 2, 3];
load = 1
capacity = 5;

# Initialize the problem
prob = RCSPP(
    adjlists,           # adjacency list representation of G
    [
        TimeWindowResource(times, ub, start_times, end_times),
        AdditiveResource(:load, load, 0, capacity),
        ElementaryResource(n_customers + 2, 2:n_customers + 1),
    ],
    1,                  # source node
    n_customers + 2,    # destination node
)
# Get all non-dominated shortest paths from source to destination, according to costs
# One can call `shortest_paths` with different values of `costs`
paths = shortest_paths(prob, costs)