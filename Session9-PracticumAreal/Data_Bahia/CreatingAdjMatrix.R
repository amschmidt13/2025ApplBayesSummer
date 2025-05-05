#--- Generating the neighborhood matrix

library(geobr)
library(spdep)

shape_bahia = read_micro_region(code_micro = "BA", year = 2020)

neigh_matrix = poly2nb(shape_bahia$geom, row.names = shape_bahia$code_micro) 
# We use the polygon from the geobr package inside the poly2nb function 
# to generate the neighborhood matrix. We replicate the microregions codes as
# the rows names.

plot(neigh_matrix, coords = shape_bahia$geom)
# Now we can plot the matrix with the coordinates from the polygons to 
# visually check if it was correctly generated

W = nb2listw(neigh_matrix, style = 'B')


nbinfo = nb2WB(neigh_matrix)
nbinfo$adj # show the number of neighbors of each area

num_regions = length(neigh_matrix)
# We create a vector with the number of regions

adj_matrix <- matrix(0, nrow = num_regions, ncol = num_regions)

for (i in 1:num_regions) {
  for (j in neigh_matrix[[i]]) {
    adj_matrix[i, j] = 1
  }
}

colnames(adj_matrix) = shape_bahia$code_micro
rownames(adj_matrix) = shape_bahia$code_micro

# Finally we change the names of columns and rows for better understanding and 
# then our adjacency matrix is completed