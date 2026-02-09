import bpy

def create_triangle_mesh(name, vertices_coords):
    faces = []
    
    for i in range(0, len(vertices_coords), 3):
        faces.append([i, i + 1, i + 2])
    
    create_mesh(name, vertices_coords, faces)
    
def create_mesh(name, vertices_coords, triangle_faces):
    mesh = bpy.data.meshes.new(name + "_Mesh")
    obj = bpy.data.objects.new(name, mesh)
    
    collection = bpy.context.scene.collection
    collection.objects.link(obj)
    
    mesh.from_pydata(vertices_coords, [], triangle_faces)
    
    # Update the mesh and calculate edges/normals
    mesh.update(calc_edges=True)
    
    return obj


verts = [
    (1, 1, 1),   # Vertex 0
    (-1, 1, -1), # Vertex 1
    (-1, -1, 1), # Vertex 2
    (1, -1, -1)  # Vertex 3
]


tetra_object = create_triangle_mesh("Tetrahedron", verts)


